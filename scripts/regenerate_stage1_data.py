"""Regenerate Stage 1 datasets with a fixed seed.

The datasets in data/original-research were produced without a seed and cannot be rebuilt
bit-for-bit.  This script produces seeded counterparts (written to data/original-research/regenerated):

  calibration_pair_conventions.csv  mean/sd of AD after t_mix for n = 2..4096 and the validation
                                    sizes, R = 1000, for BOTH pair conventions  (RSI Fig. 2)
  trajectory_n512.csv               mean/sd of AD(t) for n = 512 from t = 0 to 2 t_mix  (RSI Fig. 3)
  variance_grid_seeded.csv          the full (n, m, tau) grid, R = 1000                 (RSI Figs. 6-7)

    python scripts/regenerate_stage1_data.py [--seed 20261115] [--skip-grid] [--share-decks]

The grid takes about four minutes with the original, fully independent procedure.  --share-decks reuses
the shuffled decks across m (see rsi_shuffle.grid.run_grid) and is several times faster.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from rsi_shuffle import grid, shuffle, theory

OUT = Path(__file__).resolve().parents[1] / "data" / "original-research" / "regenerated"
SIZES = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
VALIDATION = [10, 50, 100, 300, 500, 1000, 3000]
R = 1000


def calibration(seed: int) -> None:
    rng = np.random.default_rng(seed)
    rows = []
    for n in sorted(SIZES + VALIDATION):
        T = shuffle.mixing_time(n)
        row = {"n": n, "t_mix": T, "AD_inf": theory.ad_inf(n)}
        for conv in ("independent", "distinct"):
            ad = shuffle.average_displacement(shuffle.shuffle_decks(n, T, R, rng, conv))
            row[f"mean_{conv}"], row[f"sd_{conv}"] = float(ad.mean()), float(ad.std(ddof=1))
        rows.append(row)
        print("calibration", n, flush=True)
    with open(OUT / "calibration_pair_conventions.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def trajectory(seed: int) -> None:
    rng = np.random.default_rng(seed + 1)
    n = 512
    t_mix = shuffle.mixing_time(n)
    times = sorted({round(f * t_mix) for f in np.arange(0, 2.0001, 0.05)})
    deck = np.tile(np.arange(n, dtype=np.int32), (R, 1))
    rows, rows_idx, t_prev = [], np.arange(R), 0
    for t in times:
        for _ in range(t - t_prev):
            a, b = rng.integers(0, n, R), rng.integers(0, n, R)
            ca = deck[rows_idx, a].copy(); deck[rows_idx, a] = deck[rows_idx, b]; deck[rows_idx, b] = ca
        t_prev = t
        perm = np.empty_like(deck)
        np.put_along_axis(perm, deck.astype(np.int64), np.broadcast_to(np.arange(n, dtype=np.int32), (R, n)), 1)
        ad = shuffle.average_displacement(perm)
        rows.append({"t": t, "t_over_tmix": t / t_mix, "mean_AD": float(ad.mean()), "sd_AD": float(ad.std(ddof=1)),
                     "AD_over_ADinf": float(ad.mean() / theory.ad_inf(n)),
                     "phi_exact_mean_law": theory.mean_law_phi(t, n)})
    with open(OUT / "trajectory_n512.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("trajectory done", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20261115)
    ap.add_argument("--skip-grid", action="store_true")
    ap.add_argument("--share-decks", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    calibration(a.seed)
    trajectory(a.seed)
    if not a.skip_grid:
        grid.run_grid(R=R, seed=a.seed + 2, out=OUT / "variance_grid_seeded.csv",
                      share_decks=a.share_decks)
        print("grid done")


if __name__ == "__main__":
    main()
