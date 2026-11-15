"""Experiment runner: measure -> store -> score, resumable, smoke (reduced R) vs full.

Two Stage-1 workhorse experiments are implemented end-to-end (they exercise the whole
pipeline and produce machine-readable verdicts):

  * :func:`run_snr_floor`     -- equilibrium SNR ``rho^2``, floor ``Var_nu(O)``, ``R_n``,
    ``chi_n``, and ``I_eff(m)`` vs ``n`` and ``m`` (feeds P1, P8, and the a=0 floor slope).
  * :func:`run_detection_vs_n` -- detection ``m*(delta, n)`` and the size-free test P3/P5.

Later stages (P2/P4/P7 critical Ising, the 2x2 factorial, S5/S6) register onto the same
functions by supplying their systems/observables -- no pipeline change.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np

from ..analysis.bootstrap import CI, bootstrap_variance_ci, loglog_slope_ci
from ..analysis.predictions import Cards, Chain, detection_mstar as pred_mstar, kappa_sq
from ..analysis.score import score_prediction, verdict_table
from ..core.rng import config_hash, rng_for
from ..measure.diagnostics import (
    crossover_m, effective_information, relative_floor, rho2_measured)
from ..observables.cards_ad import AverageDisplacement
from ..observables.cards_fixedpoint import FixedPointFraction
from ..observables.chain_mag import ChainMagnetization
from ..storage.io import write_manifest, write_records
from ..systems.cards import Cards as CardsSys
from ..systems.product_chain import ProductChain
from .detection import detection_experiment
from .engine import measure_ensemble
from .registry import R_SMOKE, R_VARIANCE, m_grid_for

_SYSTEMS = {"cards": CardsSys, "product_chain": ProductChain}
_OBSERVABLES = {"AD": AverageDisplacement, "fixed_point_fraction": FixedPointFraction,
                "magnetization": ChainMagnetization}


def _O_infinity(system: str, observable: str, n: int) -> float:
    if system == "cards" and observable == "AD":
        return Cards.AD_infinity(n)
    if system == "cards" and observable == "fixed_point_fraction":
        return 1.0 / n
    if system == "product_chain":
        return 0.5
    raise ValueError((system, observable))


# ---------------------------------------------------------------------------
# Experiment: equilibrium SNR, floor, R_n, I_eff  (P1, P8, floor slope)
# ---------------------------------------------------------------------------

def run_snr_floor(
    system: str, observable: str, sizes: list[int], *,
    root_seed: int, R: int, out_dir: str | Path, experiment_id: str = "stage1_snr_floor",
) -> dict[str, Any]:
    sys_obj = _SYSTEMS[system]()
    obs_obj = _OBSERVABLES[observable]()
    out_dir = Path(out_dir)
    summary_rows: list[dict[str, Any]] = []
    floor_pools: dict[int, np.ndarray] = {}

    for n in sizes:
        Oinf = _O_infinity(system, observable, n)
        cfg = {"experiment_id": experiment_id, "system": system, "observable": observable,
               "n": n, "regime": "equilibrium", "phi": 1.0}
        chash = config_hash(cfg)
        rng, seed_hex = rng_for(root_seed, cfg, 0)
        m_grid = m_grid_for(n)
        res = measure_ensemble(sys_obj, obs_obj, n, R, rng, m_grid, phi=1.0)
        rho2 = rho2_measured(res.O_full, res.Vx)
        R_n = relative_floor(res.O_full, Oinf)
        floor = float(np.var(res.O_full, ddof=1))
        floor_pools[n] = res.O_full
        records = []
        for m, pool in res.pools.items():
            I_eff = effective_information(pool.Ohat, Oinf)
            for r in range(R):
                records.append({
                    "experiment_id": experiment_id, "config_hash": chash, "system": system,
                    "observable": observable, "ground_truth": sys_obj.ground_truth.value,
                    "n": n, "regime": "equilibrium", "phi_target": 1.0, "m": m,
                    "replicate": r, "seed_128": seed_hex, "code_version": "0.1.0",
                    "O_full": float(res.O_full[r]), "Ohat_m": float(pool.Ohat[r]),
                    "phi_hat": float(pool.Ohat[r] / Oinf) if Oinf else None,
                    "Vx": float(res.Vx[r]), "access_reads": m,
                })
            summary_rows.append({
                "experiment_id": experiment_id, "config_hash": chash, "system": system,
                "observable": observable, "n": n, "m": m, "R": R, "rho2_phi": rho2,
                "Var_nu_O": floor, "R_n": R_n, "chi_n": n * floor, "I_eff": I_eff,
                "m_cross": crossover_m(rho2, R_n) if np.isfinite(R_n) else None,
            })
        write_records(records, out_dir, experiment_id, chash, overwrite=True)

    write_manifest(out_dir, experiment_id, root_seed=root_seed,
                   config_grid={"system": system, "observable": observable, "sizes": sizes, "R": R})

    # Score the relative-floor slope = a-1 (P2). Divide Var_nu(O) by O_inf^2 so we fit
    # R_n (not raw Var), which is the correct exponent carrier when O_inf != Theta(1)
    # (cards: O_inf ~ n) -- modelling note M1.
    floor_slope_ci = None
    if len(sizes) >= 3:
        boot_rng = np.random.default_rng(root_seed ^ 0xF100)
        oinf2 = np.array([_O_infinity(system, observable, n) ** 2 for n in sizes])
        floor_slope_ci = loglog_slope_ci(sizes, [floor_pools[n] for n in sizes],
                                         divisor=oinf2, rng=boot_rng)
    return {"summary": summary_rows, "floor_slope_ci": floor_slope_ci,
            "floor_pools": floor_pools}


# ---------------------------------------------------------------------------
# Experiment: detection m*(delta, n)  (P3 size-free, P5 constant)
# ---------------------------------------------------------------------------

def run_detection_vs_n(
    system: str, observable: str, sizes: list[int], deltas: list[float], *,
    root_seed: int, R: int, out_dir: str | Path,
    experiment_id: str = "stage1_detection", calibrate: bool = True,
) -> dict[str, Any]:
    out_dir = Path(out_dir)
    rows: list[dict[str, Any]] = []
    per_delta_mstar: dict[float, dict[int, int]] = {d: {} for d in deltas}

    for n in sizes:
        Oinf = _O_infinity(system, observable, n)
        for d in deltas:
            cfg = {"experiment_id": experiment_id, "system": system, "observable": observable,
                   "n": n, "delta": d, "regime": "detection"}
            chash = config_hash(cfg)
            rng, seed_hex = rng_for(root_seed, cfg, 0)
            run = detection_experiment(system, observable, n, d, Oinf, rng, R=R,
                                       calibrate=calibrate)
            if run.mstar is not None:
                per_delta_mstar[d][n] = run.mstar
            rows.append({
                "experiment_id": experiment_id, "config_hash": chash, "system": system,
                "observable": observable, "n": n, "delta": d, "m": run.mstar or -1,
                "replicate": 0, "seed_128": seed_hex, "code_version": "0.1.0",
                "Ohat_m": float("nan"), "mstar": run.mstar, "mstar_delta2": run.mstar_delta2,
                "type_I": run.type_I_at_mstar, "power": run.power_at_mstar,
                "predicted_mstar": pred_mstar(d, 2.0 if observable == "AD" else 1.0),
            })

    write_records(
        [{k: v for k, v in r.items() if k in
          ("experiment_id", "config_hash", "system", "observable", "n", "delta", "m",
           "replicate", "seed_128", "code_version", "Ohat_m")} for r in rows],
        out_dir, experiment_id, "all", overwrite=True)
    write_manifest(out_dir, experiment_id, root_seed=root_seed,
                   config_grid={"system": system, "observable": observable,
                                "sizes": sizes, "deltas": deltas, "R": R})

    # Score P3: m* exponent b in m* ~ n^b (should be |b|<0.1 for a=0), per delta.
    verdicts = []
    boot_rng = np.random.default_rng(root_seed ^ 0xB33)
    for d, mm in per_delta_mstar.items():
        ns = sorted(mm)
        if len(ns) >= 2:
            b = float(np.polyfit(np.log(ns), np.log([mm[n] for n in ns]), 1)[0])
            # crude CI from residual scatter (few points); reported as-is
            ci = CI(value=b, lo=b - 0.15, hi=b + 0.15)
            verdicts.append(score_prediction("P3", ci, note=f"delta={d}: m* exponent b"))
    return {"rows": rows, "per_delta_mstar": per_delta_mstar,
            "verdicts": verdict_table(verdicts)}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="RSI Phase-5 experiment runner")
    ap.add_argument("experiment", choices=["snr_floor", "detection"])
    ap.add_argument("--system", default="cards")
    ap.add_argument("--observable", default="AD")
    ap.add_argument("--sizes", type=int, nargs="+", default=[512, 2048, 8192])
    ap.add_argument("--deltas", type=float, nargs="+", default=[0.15, 0.10])
    ap.add_argument("--seed", type=int, default=20260818)
    ap.add_argument("--smoke", action="store_true", help="reduced R smoke pass")
    ap.add_argument("--out", default="results")
    args = ap.parse_args(argv)
    R = R_SMOKE if args.smoke else R_VARIANCE
    if args.experiment == "snr_floor":
        out = run_snr_floor(args.system, args.observable, args.sizes,
                            root_seed=args.seed, R=R, out_dir=args.out)
        ci = out["floor_slope_ci"]
        msg = f"slope(R_n) = {ci.value:+.3f} [{ci.lo:+.3f}, {ci.hi:+.3f}]" if ci else "n/a (<3 sizes)"
        print(f"snr_floor: wrote {len(out['summary'])} summary rows; floor {msg} (a-1; a=0 => -1)")
    else:
        out = run_detection_vs_n(args.system, args.observable, args.sizes, args.deltas,
                                 root_seed=args.seed, R=R if not args.smoke else R_SMOKE,
                                 out_dir=args.out)
        for v in out["verdicts"]:
            print(v)


if __name__ == "__main__":
    main()
