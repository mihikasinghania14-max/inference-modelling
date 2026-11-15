"""Figure factory -- publication-quality, colorblind-safe, generated FROM STORED TABLES
only (never from live runs; design brief sec.8).

Each function takes tidy data (arrays or a summary frame) and returns a Matplotlib figure
following the RSI / Phase-3 visual grammar: log axes where appropriate, the rho^2(phi)
curves, precision-vs-n, m* vs delta / vs n, the Ising variance-vs-m floor, and the
information-saturation curve.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

# Okabe-Ito colorblind-safe palette
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442"]


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 130, "font.size": 10,
                         "axes.grid": True, "grid.alpha": 0.3})
    return plt


def fig_rho2_of_phi(phi: Sequence[float], measured: dict[str, Sequence[float]],
                    predicted: dict[str, Sequence[float]] | None = None):
    """rho^2(phi) trajectory curves (Phase-3 Panel A)."""
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(5, 4))
    for i, (name, ys) in enumerate(measured.items()):
        c = PALETTE[i % len(PALETTE)]
        ax.plot(phi, ys, "o", color=c, label=f"{name} (measured)")
        if predicted and name in predicted:
            ax.plot(phi, predicted[name], "-", color=c, alpha=0.7, label=f"{name} (theory)")
    ax.set_xlabel(r"mixing level $\phi$")
    ax.set_ylabel(r"$\rho^2(\phi)$")
    ax.set_title(r"Off-equilibrium signal-to-noise $\rho^2(\phi)$")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def fig_precision_vs_n(n: Sequence[int], rho2_by_n: Sequence[float], target: float,
                       label: str = "measured"):
    """Equilibrium rho^2 vs n -- flat line demonstrates size-independence (Phase-3 Panel B)."""
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.semilogx(n, rho2_by_n, "o-", color=PALETTE[0], base=2, label=label)
    ax.axhline(target, ls="--", color="k", alpha=0.6, label=f"target {target:g}")
    ax.set_xlabel("system size $n$")
    ax.set_ylabel(r"$\rho^2$")
    ax.set_title(r"Size-independence of $\rho^2$")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def fig_mstar_vs_n(n: Sequence[int], mstar_by_delta: dict[float, Sequence[float]]):
    """m* vs n at several delta -- flat curves demonstrate size-free certification (P3)."""
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(5, 4))
    for i, (d, ms) in enumerate(sorted(mstar_by_delta.items())):
        ax.loglog(n, ms, "o-", base=2, color=PALETTE[i % len(PALETTE)], label=fr"$\delta={d}$")
    ax.set_xlabel("system size $n$")
    ax.set_ylabel(r"$m^*$")
    ax.set_title(r"Certification cost $m^*(\delta,n)$")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def fig_variance_vs_m(m: Sequence[int], var_by_series: dict[str, Sequence[float]],
                      site_var: dict[str, float] | None = None):
    """Var(Ohat_m) vs m -- 1/m fall vs an m-independent floor (Phase-3 Panel D; P1)."""
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(5, 4))
    m = np.asarray(m)
    for i, (name, v) in enumerate(var_by_series.items()):
        c = PALETTE[i % len(PALETTE)]
        ax.loglog(m, v, "o-", color=c, label=name)
        if site_var and name in site_var:
            ax.loglog(m, site_var[name] / m, "--", color=c, alpha=0.6)
    ax.set_xlabel("sample size $m$")
    ax.set_ylabel(r"$\mathrm{Var}(\hat O_m)$")
    ax.set_title(r"Sampling variance vs $m$: $1/m$ decay vs floor")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def fig_information_saturation(m: Sequence[int], I_eff: Sequence[float],
                               I_max: float | None = None, m_cross: float | None = None):
    """I_eff(m) rising then saturating at 1/R_n (P8)."""
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.loglog(m, I_eff, "o-", color=PALETTE[2], label=r"$I_{\rm eff}(m)=1/CV^2$")
    if I_max:
        ax.axhline(I_max, ls="--", color="k", alpha=0.6, label=r"$I_{\max}=1/R_n$")
    if m_cross:
        ax.axvline(m_cross, ls=":", color="gray", alpha=0.7, label=r"$m_\times$")
    ax.set_xlabel("sample size $m$")
    ax.set_ylabel(r"effective information $I_{\rm eff}$")
    ax.set_title("Information saturation")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def save(fig, path: str) -> str:
    fig.savefig(path, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return path
