"""Redraw summary figures of the original RSI study from the stored data.

Reads only files under ``data/original-research`` and writes PNGs to
``results/figures/original-research``.  These are *new renderings* of stored data, not
the figures of the paper (which are in ``research/original-rsi-research/figures``).

    A. ad_vs_n.png            mean AD after t_mix against AD_inf = (n^2-1)/(3n)   (figure2 data)
    B. variance_law_tau1.png  Var(AD_hat_m) / (n^2/(18 m)) against m/n, tau = 1     (variance grid)
    C. variance_vs_m_n512.png Var(AD_hat_m) against m for n = 512 at tau = 1, with the exact law

Usage:  python scripts/make_original_figures.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rsi_shuffle import theory

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "original-research"
OUT = ROOT / "results" / "figures" / "original-research"


def fig_ad_vs_n() -> None:
    """figure2 sheet: main series (n = 2..4096), a blank row, a caption row, then a second
    'validation' series at n = 10..3000 (plotted in a different colour, as in the sheet)."""
    raw = pd.read_csv(DATA / "workbooks/csv/figure2__Sheet1.csv", header=None, skiprows=1,
                      names=["n", "ad_inf", "sim", "sd"])
    num = raw.apply(pd.to_numeric, errors="coerce")
    blank = int(num["n"].isna().to_numpy().nonzero()[0][0])
    main, val = num.iloc[:blank], num.iloc[blank:].dropna(subset=["n"])
    main = main[main["n"] > 0]
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    nn = np.logspace(np.log10(main["n"].min()), np.log10(max(main["n"].max(), val["n"].max())), 200)
    ax.loglog(nn, (nn**2 - 1) / (3 * nn), color="0.3", lw=1.2, label=r"$AD_\infty=(n^2-1)/(3n)$")
    ax.errorbar(main["n"], main["sim"], yerr=main["sd"], fmt="o", ms=4, color="#0072B2", capsize=2,
                label=r"simulated mean at $t_{mix}$ ($\pm$ sd)")
    ax.errorbar(val["n"], val["sim"], yerr=val["sd"], fmt="s", ms=4, color="#D55E00", capsize=2,
                label="validation sizes")
    ax.set_xlabel("n"); ax.set_ylabel("average displacement")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "ad_vs_n.png", dpi=200); plt.close(fig)


def _grid() -> pd.DataFrame:
    return pd.read_csv(DATA / "variance-grid/csv/simulation_results__simulation_results.csv")


def fig_variance_law() -> None:
    g = _grid()
    g = g[(g.tau == 1.0) & (g.N >= 16)].copy()
    g["ratio"] = g["variance"] / (g.N**2 / (18 * g.M))
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    for n, s in g.groupby("N"):
        ax.plot(s.sampling_fraction, s.ratio, "o-", ms=3, lw=0.8, label=f"n={n}")
    f = np.linspace(0.01, 1, 200)
    ax.plot(f, 1 - f / 5, "k--", lw=1, label=r"exact law, $n\to\infty$: $1-\frac{1}{5}\,m/n$")
    ax.axhline(1, color="0.6", lw=0.8)
    ax.set_xscale("log"); ax.set_xlabel("m / n"); ax.set_ylabel(r"Var$/ (n^2/18m)$")
    ax.legend(frameon=False, fontsize=6, ncol=2)
    fig.tight_layout(); fig.savefig(OUT / "variance_law_tau1.png", dpi=200); plt.close(fig)


def fig_variance_vs_m() -> None:
    g = _grid()
    s = g[(g.tau == 1.0) & (g.N == 512)]
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    ax.loglog(s.M, s.variance, "o", color="#0072B2", label="simulation (R=1000)")
    m = np.logspace(np.log10(s.M.min()), np.log10(s.M.max()), 100)
    ax.loglog(m, [theory.estimator_variance(512, int(round(x))) for x in m], color="0.2",
              label="exact law (Lemma 1)")
    ax.loglog(m, 512**2 / (18 * m), "--", color="0.5", label=r"$n^2/(18m)$")
    ax.set_xlabel("m"); ax.set_ylabel(r"Var$(\widehat{AD}_m)$, n = 512")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "variance_vs_m_n512.png", dpi=200); plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_ad_vs_n(); fig_variance_law(); fig_variance_vs_m()
    print("wrote", *sorted(p.name for p in OUT.glob("*.png")))
