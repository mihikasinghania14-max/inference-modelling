"""Curie-Weiss floor-exponent experiment (CX1 / P2) -- the a=1 hard-floor result.

Measures the equilibrium magnetization floor ``Var_nu(M)`` vs system size ``N`` and fits
the exponent ``a-1`` (slope of ``log Var_nu(M)`` vs ``log N``; since ``chi_N = N Var(M) ~
N^a``, the raw-variance slope IS ``a-1`` with no ``O_inf`` needed -- appropriate here as
the symmetric ensemble has ``E[M]=0``).  Prediction (CX1): below ``T_c`` the susceptibility
is extensive (``a=1`` -> slope ``0`` -> a non-vanishing precision floor), while above
``T_c`` it self-averages (``a=0`` -> slope ``-1``).

Exactness / efficiency: mean-field spins are **exchangeable given the total magnetization**,
so we sample ``M`` from its exact 1-D law and subsample via the **hypergeometric**
distribution -- the exact exchangeable fast path (analogous to
``veil.reveal_exchangeable_first_m``).  This reaches ``N=10^5`` with no ``O(NR)`` allocation.
A test cross-checks this against the materialized-spins + veil path at small ``N``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln

from ..analysis.bootstrap import CI, loglog_slope_ci
from ..analysis.score import Verdict, score_prediction


def _logP_k(N: int, T: float, J: float = 1.0) -> np.ndarray:
    k = np.arange(N + 1)
    m = (2.0 * k - N) / N
    return gammaln(N + 1) - gammaln(k + 1) - gammaln(N - k + 1) + (J / T) * N * m * m / 2.0


def sample_up_counts(N: int, T: float, R: int, rng: np.random.Generator) -> np.ndarray:
    """Exact equilibrium draws of the up-spin count ``k`` (``R`` samples)."""
    lp = _logP_k(N, T); lp -= lp.max()
    p = np.exp(lp); p /= p.sum()
    return rng.choice(N + 1, size=R, p=p)


def magnetizations(N: int, T: float, R: int, rng: np.random.Generator) -> np.ndarray:
    k = sample_up_counts(N, T, R, rng)
    return (2.0 * k - N) / N


def subsample_Mhat(k: np.ndarray, N: int, m: int, rng: np.random.Generator) -> np.ndarray:
    """Exact ``M_hat_m`` from ``m`` uniformly-sampled spins of configs with ``k`` up-spins,
    via the hypergeometric law (exchangeable fast path)."""
    up = rng.hypergeometric(k, N - k, m)
    return (2.0 * up - m) / m


@dataclass
class FloorScaling:
    T: float
    Ns: list[int]
    floors: list[float]
    slope_ci: CI                 # slope of log Var(M) vs log N == a - 1
    a_estimate: float            # slope + 1
    verdict: Verdict


def floor_exponent(
    T: float,
    Ns: list[int],
    rng: np.random.Generator,
    *,
    R: int = 20000,
    n_boot: int = 1000,
    predicted_a: float | None = None,
) -> FloorScaling:
    """Measure ``Var_nu(M)`` vs ``N`` and fit ``a-1`` with a bootstrap CI, then score P2."""
    pools = [magnetizations(N, T, R, rng) for N in Ns]
    floors = [float(np.var(p, ddof=1)) for p in pools]
    slope_ci = loglog_slope_ci(np.array(Ns), pools, divisor=None, n_boot=n_boot,
                               rng=np.random.default_rng(rng.integers(1 << 32)))
    a_est = slope_ci.value + 1.0
    if predicted_a is None:
        predicted_a = 1.0 if T < 1.0 else 0.0
    band = (predicted_a - 1.0 - 0.1, predicted_a - 1.0 + 0.1)   # slope band around a-1
    verdict = score_prediction("P2", slope_ci, support=band,
                               note=f"Curie-Weiss T={T}: slope=a-1, predicted a={predicted_a}")
    return FloorScaling(T=T, Ns=Ns, floors=floors, slope_ci=slope_ci,
                        a_estimate=a_est, verdict=verdict)


@dataclass
class HardFloor:
    T: float
    N: int
    m_grid: list[int]
    var_Mhat: list[float]
    floor: float                 # Var_nu(M)
    ratios: list[float]          # Var(Mhat_m)/floor -- ~1 for a hard floor


def hard_floor_vs_m(
    T: float, N: int, m_grid: list[int], rng: np.random.Generator, *, R: int = 20000,
) -> HardFloor:
    """Show that subsampling cannot beat the floor below ``T_c``: ``Var(M_hat_m)`` stays at
    ``Var_nu(M)`` for all ``m`` (the a=1 hard floor)."""
    k = sample_up_counts(N, T, R, rng)
    floor = float(np.var((2.0 * k - N) / N, ddof=1))
    var_Mhat = [float(np.var(subsample_Mhat(k, N, m, rng), ddof=1)) for m in m_grid]
    return HardFloor(T=T, N=N, m_grid=m_grid, var_Mhat=var_Mhat, floor=floor,
                     ratios=[v / floor for v in var_Mhat])
