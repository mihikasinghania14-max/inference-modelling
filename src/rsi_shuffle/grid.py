"""The (n, m, tau) Monte Carlo grid behind the variance-law figures.

This is a vectorised re-implementation of ``bayesian_simulation.py`` (preserved
unchanged in ``experiments/original-research/legacy-scripts``).  The experiment is the
same: for each combination ``(n, m, tau)``

1. start from the ordered deck, apply ``T = round(tau * t_mix)`` random transpositions
   (independent-pair convention by default);
2. draw a uniform random ``m``-subset of cards without replacement;
3. record the mean displacement of the sampled cards;

and repeat for ``R`` independent replicates, reporting mean, unbiased variance, sd,
standard error, a normal 95% interval for the mean, skewness and excess kurtosis.

Differences from the original script (none affect the model):

* replicates are shuffled in a vectorised batch;
* an explicit ``seed`` makes a run reproducible (the original used fresh OS entropy);
* the default grids are those of the original script, which differ slightly from the
  table in RSI Appendix C.3 (see ``docs/known-issues.md``, item D3).
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.stats import kurtosis, skew

from .estimator import uniform_sample_estimate
from .shuffle import PairConvention, average_displacement, displacements, mixing_time, shuffle_decks

N_VALUES = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
M_VALUES = {
    2: [1, 2], 4: [1, 2, 3, 4], 8: [1, 2, 4, 6, 8], 16: [2, 4, 6, 8, 12, 16],
    32: [2, 4, 8, 12, 16, 24, 32], 64: [4, 8, 16, 24, 32, 48, 64],
    128: [8, 16, 32, 48, 64, 96, 128], 256: [16, 32, 64, 96, 128, 192, 256],
    512: [32, 64, 96, 128, 192, 256, 384, 512],
    1024: [64, 128, 192, 256, 384, 512, 768, 1024],
    2048: [128, 256, 384, 512, 768, 1024, 1536, 2048],
    4096: [256, 512, 768, 1024, 1536, 2048, 3072, 4096],
}
TAU_VALUES = [round(i * 0.1, 1) for i in range(21)]
COLUMNS = ["N", "M", "sampling_fraction", "t_mix", "tau", "T", "mean_AD", "variance",
           "standard_deviation", "standard_error", "CI_lower", "CI_upper", "skewness", "kurtosis"]


def summarise(ad_values: np.ndarray) -> dict[str, float]:
    r = len(ad_values)
    mean = float(np.mean(ad_values))
    var = float(np.var(ad_values, ddof=1))
    sd = math.sqrt(var)
    se = sd / math.sqrt(r)
    return {"mean_AD": mean, "variance": var, "standard_deviation": sd, "standard_error": se,
            "CI_lower": mean - 1.96 * se, "CI_upper": mean + 1.96 * se,
            "skewness": float(skew(ad_values)), "kurtosis": float(kurtosis(ad_values, fisher=True))}


def run_configuration(n: int, m: int, T: int, R: int, rng: np.random.Generator,
                      pair_convention: PairConvention = "independent") -> np.ndarray:
    """``R`` independent values of the partial estimator ``AD_hat_m`` after ``T`` steps."""
    perm = shuffle_decks(n, T, R, rng, pair_convention)
    d = displacements(perm)
    return np.asarray(uniform_sample_estimate(d, m, rng), dtype=float)


def run_grid(sizes: Iterable[int] = N_VALUES, taus: Iterable[float] = TAU_VALUES,
             R: int = 1000, seed: int | None = None, out: str | Path | None = None,
             pair_convention: PairConvention = "independent",
             m_values: dict[int, list[int]] | None = None, progress: bool = False,
             share_decks: bool = False) -> list[dict]:
    """Run the grid and optionally write the CSV; returns the rows.

    With ``share_decks=True`` the R shuffled decks for a given (n, tau) are generated once and
    every m draws its own independent random subsets from them.  Each row keeps exactly the
    distribution of the original procedure, but rows with the same (n, tau) are no longer
    independent of one another.  It is several times faster for the full grid.
    """
    m_values = m_values or M_VALUES
    taus, sizes = list(taus), list(sizes)
    rng = np.random.default_rng(seed)
    rows: list[dict] = []

    def record(n: int, m: int, tau: float, t_mix: int, T: int, vals: np.ndarray) -> None:
        if R >= 2 and vals.std() > 0:
            stats = summarise(vals)
        else:                                   # zero variance: moments of order 3, 4 undefined
            stats = {k: (0.0 if k not in ("skewness", "kurtosis") else float("nan"))
                     for k in COLUMNS[6:]}
            stats["mean_AD"] = float(vals.mean())
        rows.append({"N": n, "M": m, "sampling_fraction": round(m / n, 6), "t_mix": t_mix,
                     "tau": tau, "T": T,
                     **{k: round(v, 8) if v == v else v for k, v in stats.items()}})
        if progress:
            print(f"N={n} M={m} tau={tau} T={T}", flush=True)

    for n in sizes:
        t_mix = mixing_time(n)
        if share_decks:
            for tau in taus:
                T = round(tau * t_mix)
                d = displacements(shuffle_decks(n, T, R, rng, pair_convention))
                for m in m_values[n]:
                    vals = np.asarray(uniform_sample_estimate(d, m, rng), dtype=float)
                    record(n, m, tau, t_mix, T, vals)
        else:
            for m in m_values[n]:
                for tau in taus:
                    T = round(tau * t_mix)
                    vals = run_configuration(n, m, T, R, rng, pair_convention)
                    record(n, m, tau, t_mix, T, vals)
    rows.sort(key=lambda r: (r["N"], r["M"], r["tau"]))
    if out is not None:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS)
            w.writeheader()
            w.writerows(rows)
    return rows


def main(argv: list[str] | None = None) -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sizes", type=int, nargs="+", default=N_VALUES)
    ap.add_argument("--taus", type=float, nargs="+", default=TAU_VALUES)
    ap.add_argument("-R", "--replicates", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=None, help="omit for fresh OS entropy (as in the original)")
    ap.add_argument("--pair-convention", choices=["independent", "distinct"], default="independent")
    ap.add_argument("--share-decks", action="store_true",
                    help="reuse the shuffled decks of each (n, tau) for all m (faster; see run_grid)")
    ap.add_argument("--out", default="simulation_results.csv")
    a = ap.parse_args(argv)
    rows = run_grid(a.sizes, a.taus, a.replicates, a.seed, a.out, a.pair_convention, progress=True,
                    share_decks=a.share_decks)
    print(f"wrote {len(rows)} rows to {a.out}")


if __name__ == "__main__":
    main()
