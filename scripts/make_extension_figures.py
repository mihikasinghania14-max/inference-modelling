"""Figures for the extension, drawn from exact formulas and from the data in this repository.

    A. decomposition.png     CV^2 of the estimator against m for the shuffle (n = 1024): the exact
                             law, its 1/m sampling term and the floor R_n         (rsi_shuffle.theory)
    B. floor_vs_n.png        R_n and n R_n for the shuffle, with the limit 2/5      (rsi_shuffle.theory)
    C. certification_cost.png  m*(delta) with the floor: divergence at delta = kappa sqrt(R_n)
    D. ising_floors.png      relative floors of the energy density and of |M| at T_c: exact torus energy
                             floor (Kaufman), exact small-torus |M| floor, regenerated Wolff runs
    E. curie_weiss.png       exact Var(M) against N below, at and above T_c (a - 1 = 1, 1/2, 0 slopes)
    F. field_coverage.png    coverage of the naive interval for a correlated field against xi/L

Writes PNGs to results/figures/extension.   python scripts/make_extension_figures.py
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rsi_shuffle import theory

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figures" / "extension"
NAVY, BRICK, GREY = "#1F3A5F", "#A23B2A", "0.45"
plt.rcParams.update({"font.family": "serif", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150, "savefig.dpi": 200})


def save(fig, name):
    fig.tight_layout(); fig.savefig(OUT / name); plt.close(fig)


def decomposition():
    n = 1024
    m = np.unique(np.round(np.logspace(0, np.log10(n), 120)).astype(int))
    ad2 = float(theory.ad_inf(n)) ** 2
    cv2 = np.array([theory.exact_cv2(n, int(k)) for k in m])
    floor = float(theory.relative_floor(n))
    samp = (cv2 - floor)
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    ax.loglog(m, cv2, color="k", lw=1.6, label=r"CV$^2(\hat{AD}_m)$, exact")
    ax.loglog(m, samp, color=NAVY, lw=1.2, ls="--", label="sampling term")
    ax.axhline(floor, color=BRICK, lw=1.4, label=rf"floor $R_n$ = {floor:.2e}")
    ax.loglog(m, 1 / (2 * m), color=GREY, lw=0.8, ls=":", label=r"$1/(2m)$")
    ax.set_xlabel("m"); ax.set_ylabel(r"squared coefficient of variation")
    ax.set_title("Random transposition shuffle, n = 1024", fontsize=9)
    ax.legend(frameon=False, fontsize=7.5)
    save(fig, "decomposition.png")


def floor_vs_n():
    n = np.unique(np.round(np.logspace(np.log10(3), 4, 80)).astype(int))
    R = np.array([theory.relative_floor(int(k)) for k in n])
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.0))
    ax[0].loglog(n, R, color=BRICK); ax[0].set_xlabel("n"); ax[0].set_ylabel(r"$R_n$")
    ax[0].loglog(n, 0.4 / n, color=GREY, ls=":", lw=0.9, label=r"$2/(5n)$"); ax[0].legend(frameon=False, fontsize=8)
    ax[1].semilogx(n, n * R, color=BRICK); ax[1].axhline(0.4, color=GREY, ls=":", lw=0.9)
    ax[1].set_xlabel("n"); ax[1].set_ylabel(r"$n R_n$"); ax[1].set_ylim(0, None)
    ax[1].text(4, 0.41, "limit 2/5", fontsize=8, color=GREY)
    save(fig, "floor_vs_n.png")


def certification_cost():
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    for n, c in [(256, NAVY), (1024, "#5B7DB1"), (10**5, BRICK)]:
        d = np.linspace(0.01, 0.4, 600)
        cost = np.minimum([theory.certification_cost(x, n) for x in d], n)
        ax.semilogy(d, cost, color=c, label=f"n = {n:,}")
        ax.axvline(theory.kappa() * np.sqrt(theory.relative_floor(n)), color=c, ls=":", lw=0.8)
    ax.set_xlabel(r"relative tolerance $\delta$"); ax.set_ylabel(r"$m^*(\delta)$, capped at $n$")
    ax.set_ylim(1, 3e5); ax.legend(frameon=False, fontsize=8)
    ax.set_title(r"dotted: $\delta=\kappa\sqrt{R_n}$ (certification impossible below)", fontsize=8)
    save(fig, "certification_cost.png")


def ising():
    ex = {}
    for line in (ROOT / "results/outputs/extension-verification/ising_exact.txt").read_text().splitlines():
        m = re.match(r"Kaufman L=(\d+): R\[e\]=([0-9.e+-]+)", line)
        if m:
            ex[int(m.group(1))] = float(m.group(2))
    mc = pd.read_csv(ROOT / "data/extension/regenerated/ising_critical_sizes.csv")
    fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.2), sharex=True)
    L = np.array(sorted(ex))
    ax[0].loglog(L**2, [ex[k] for k in L], "-", color=NAVY, lw=1.2, label="exact (Kaufman)")
    ax[0].errorbar(mc.n, mc.R_E, yerr=mc.R_E_se, fmt="o", ms=3.5, color=BRICK, label="Wolff, regenerated")
    ax[0].set_title(r"energy density: $R_n[E]$ vanishes", fontsize=9)
    exact_M = {16: None}
    ax[1].errorbar(mc.n, mc.R_M, yerr=mc.R_M_se, fmt="o", ms=3.5, color=BRICK, label="Wolff, regenerated")
    for Lq, v in {4: 0.06917, 6: 0.06995, 8: 0.07061}.items():
        ax[1].plot(Lq**2, v, "s", ms=4, color=NAVY, label="exact (enumeration, transfer matrix)" if Lq == 4 else None)
    ax[1].set_xscale("log"); ax[1].set_ylim(0, 0.1)
    ax[1].set_title(r"magnitude $|M|$: $R_n[|M|]$ stays flat", fontsize=9)
    for a in ax:
        a.set_xlabel(r"$n = L^2$"); a.legend(frameon=False, fontsize=7)
    ax[0].set_ylabel("relative floor"); save(fig, "ising_floors.png")


def curie_weiss():
    d = pd.read_csv(ROOT / "data/extension/regenerated/curie_weiss_exact.csv")
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    for (T, g), c in zip(d.groupby("T"), [BRICK, "#C9825B", "k", "#5B7DB1", NAVY]):
        ax.loglog(g.N, g.var_M, "o-", ms=2.5, lw=1, color=c, label=f"T = {T:g}")
    ax.set_xlabel("N"); ax.set_ylabel(r"Var$(M)$"); ax.legend(frameon=False, fontsize=7.5)
    ax.set_title("Curie-Weiss, exact: slope = a - 1", fontsize=9)
    save(fig, "curie_weiss.png")


def field_coverage():
    rows = []
    for line in (ROOT / "results/outputs/extension-verification/field.txt").read_text().splitlines():
        p = line.split()
        if len(p) == 8 and re.match(r"^\d\.\d\d$", p[0]):
            rows.append([float(x) for x in p])
    a = np.array(rows)
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    ax.plot(a[:, 0], a[:, 5], "o-", color=NAVY, label=r"interval for the field value $O(x)$")
    ax.plot(a[:, 0], a[:, 4], "o-", color=BRICK, label=r"interval for the mean $\mu$ (naive)")
    ax.plot(a[:, 0], a[:, 3], "--", color="k", lw=0.9, label="predicted from the floor")
    ax.axhline(0.95, color=GREY, lw=0.6)
    ax.set_xlabel(r"correlation length / box size, $\xi/L$"); ax.set_ylabel("coverage of nominal 95% interval")
    ax.set_ylim(0, 1.02); ax.legend(frameon=False, fontsize=7.5)
    save(fig, "field_coverage.png")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    decomposition(); floor_vs_n(); certification_cost(); curie_weiss(); field_coverage()
    if (ROOT / "data/extension/regenerated/ising_critical_sizes.csv").exists():
        ising()
    print("wrote", *sorted(p.name for p in OUT.glob("*.png")))
