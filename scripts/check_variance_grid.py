"""Re-run a reduced part of the (n, m, tau) variance grid and compare with the stored grid.

For each selected (n, m) at tau = 1 a fresh Monte Carlo run (seeded, R = 1000 by default,
the same R as the stored grid) is compared with

  * the stored value in data/original-research/variance-grid  (variance and mean), and
  * the exact closed form theory.estimator_variance(n, m).

Both runs have Monte Carlo error: with R = 1000 the relative standard error of a variance
is about sqrt(2/(R-1)) ~ 4.5%, so two independent runs differ by ~6.3% (1 sd).  The check
passes when every |log ratio| is within 4.5 combined standard errors.

Usage:  python scripts/check_variance_grid.py [--seed 0] [-R 1000]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

from rsi_shuffle import grid, theory

STORED = Path(__file__).resolve().parents[1] / "data/original-research/variance-grid/csv/simulation_results__simulation_results.csv"
SELECTION = {64: [8, 32, 64], 256: [16, 64, 256], 512: [32, 128, 512]}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-R", "--replicates", type=int, default=1000)
    a = ap.parse_args()
    st = pd.read_csv(STORED)
    rows = grid.run_grid(sizes=list(SELECTION), taus=[1.0], R=a.replicates, seed=a.seed,
                         m_values=SELECTION)
    rel_se = math.sqrt(2.0 / (a.replicates - 1))
    ok = True
    print(f"{'n':>5} {'m':>4} {'var_new':>10} {'var_stored':>10} {'var_exact':>10} "
          f"{'new/stored':>10} {'new/exact':>9}  mean_new / AD_inf")
    for r in rows:
        s = st[(st.N == r["N"]) & (st.M == r["M"]) & (st.tau == 1.0)].iloc[0]
        exact = theory.estimator_variance(r["N"], r["M"])
        q1, q2 = r["variance"] / s["variance"], r["variance"] / exact
        pass1 = abs(math.log(q1)) < 4.5 * math.sqrt(2) * rel_se
        pass2 = abs(math.log(q2)) < 4.5 * rel_se
        ok &= pass1 and pass2
        print(f"{r['N']:>5} {r['M']:>4} {r['variance']:>10.3f} {s['variance']:>10.3f} {exact:>10.3f} "
              f"{q1:>10.3f} {q2:>9.3f}  {r['mean_AD'] / theory.ad_inf(r['N']):.4f}"
              f"{'' if pass1 and pass2 else '   <-- outside tolerance'}")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
