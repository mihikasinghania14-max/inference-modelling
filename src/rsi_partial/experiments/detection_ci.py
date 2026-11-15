"""Detection ``m*`` with a bootstrap CI, a delta-sweep, and an honest oracle-vs-deployable
null split (addresses critique sec.2.1-2.3).

Why this exists
---------------
`experiments/detection.py` bisects a *stochastic* power estimate and returns a single
``m*`` with no uncertainty.  Here we instead:

  * precompute, once, the per-replicate test statistic ``T = (O_inf - Ohat)/sd`` at every
    ``m`` on a geometric grid, for an equilibrium (H0) and a relaxed (H1) ensemble;
  * locate ``m*`` by interpolating the power curve's crossing of the target;
  * **bootstrap** ``m*`` by resampling replicate indices (recomputing the per-``m``
    calibrated threshold and the crossing each time) -> a genuine CI on the headline number;
  * report **two nulls side by side**: the *oracle* calibrated null (threshold = empirical
    upper-alpha quantile of the H0 statistic -- an idealized upper bound that assumes you
    know the equilibrium sampling law) and the *deployable* Gaussian null (threshold
    ``z_{1-alpha}`` -- what you actually have in the field);
  * give an **honest Type-I** by calibrating the threshold on one half of the H0 pool and
    measuring size on the held-out half (not on the pool used to set the threshold).

Sweeping ``delta`` then shows ``m* delta^2 -> kappa^2/rho^2`` as ``delta -> 0`` with the
finite-``delta`` overshoot shrinking -- the constant is *measured*, not merely bracketed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import norm

from .detection import _make_D
from ..inference.mstar import geometric_grid
from ..observe.veil import reveal_batch


def _stat_at_m(D: np.ndarray, m: int, O_infinity: float, rng) -> np.ndarray:
    """Per-replicate one-sided statistic ``T = (O_inf - Ohat_m)/sd`` (large T => below O_inf)."""
    R, n = D.shape
    mm = min(m, n)
    vals, reads = reveal_batch(D, mm, rng)
    assert reads == mm
    Ohat = vals.mean(axis=1)
    if mm >= 2:
        s2 = vals.var(axis=1, ddof=1)
        fpc = (n - mm) / (n - 1) if n > 1 else 0.0
        sd = np.sqrt(np.maximum(s2 / mm * fpc, 0.0))
    else:
        sd = np.zeros(R)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(sd > 0, (O_infinity - Ohat) / sd, 0.0)


def _mstar_from_curve(m_grid, power, target):
    """Interpolate the smallest ``m`` where ``power(m) >= target`` (linear in ``m``).

    Returns ``np.inf`` if the target is never reached, or ``m_grid[0]`` if already met."""
    for k, m in enumerate(m_grid):
        if power[k] >= target:
            if k == 0:
                return float(m)
            p0, p1 = power[k - 1], power[k]
            if p1 <= p0:
                return float(m)
            frac = (target - p0) / (p1 - p0)
            return float(m_grid[k - 1] + frac * (m - m_grid[k - 1]))
    return np.inf


@dataclass
class MStarCI:
    delta: float
    null: str                       # 'oracle' | 'deployable'
    mstar: float
    ci_lo: float
    ci_hi: float
    mstar_delta2: float
    d2_lo: float
    d2_hi: float
    type_I: float                   # held-out (oracle) or nominal-threshold (deployable)
    censored_frac: float            # fraction of bootstraps that never reached the target
    n_boot: int


def _threshold(stat0_cal: np.ndarray, alpha: float, calibrate: bool) -> float:
    if calibrate:
        return float(np.quantile(stat0_cal, 1.0 - alpha))
    return float(norm.ppf(1.0 - alpha))


def mstar_ci_for_delta(
    stat0: dict[int, np.ndarray],
    stat1: dict[int, np.ndarray],
    m_grid: list[int],
    delta: float,
    *,
    alpha: float,
    target: float,
    calibrate: bool,
    n_boot: int,
    rng: np.random.Generator,
) -> MStarCI:
    """Bootstrap CI on ``m*`` for one delta, given precomputed per-m statistics."""
    R = next(iter(stat0.values())).shape[0]
    half = R // 2
    # point estimate: per-m calibrated (or z) threshold; Type-I held out on the 2nd half
    def curve(idx0, idx1, cal_slice, test_slice):
        power = np.empty(len(m_grid)); typeI = np.empty(len(m_grid))
        for k, m in enumerate(m_grid):
            s0 = stat0[m][idx0]; s1 = stat1[m][idx1]
            thr = _threshold(s0[cal_slice], alpha, calibrate)
            power[k] = np.mean(s1 > thr)
            typeI[k] = np.mean(s0[test_slice] > thr)
        return power, typeI

    full = np.arange(R)
    cal_slice, test_slice = slice(0, half), slice(half, R)
    power, typeI = curve(full, full, cal_slice, test_slice)
    mstar = _mstar_from_curve(m_grid, power, target)
    # Type-I reported at (the grid m nearest) the point m*, held-out half
    if np.isfinite(mstar):
        kstar = int(np.argmin([abs(m - mstar) for m in m_grid]))
    else:
        kstar = len(m_grid) - 1
    typeI_at = float(typeI[kstar])

    boot = []
    censored = 0
    for _ in range(n_boot):
        idx0 = rng.integers(0, R, R); idx1 = rng.integers(0, R, R)
        pw, _ = curve(idx0, idx1, cal_slice, test_slice)
        ms = _mstar_from_curve(m_grid, pw, target)
        if np.isfinite(ms):
            boot.append(ms)
        else:
            censored += 1
    boot = np.array(boot) if boot else np.array([np.nan])
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    d2 = mstar * delta * delta
    return MStarCI(
        delta=delta, null=("oracle" if calibrate else "deployable"),
        mstar=mstar, ci_lo=float(lo), ci_hi=float(hi),
        mstar_delta2=d2, d2_lo=float(lo * delta * delta), d2_hi=float(hi * delta * delta),
        type_I=typeI_at, censored_frac=censored / n_boot, n_boot=n_boot,
    )


@dataclass
class DeltaSweep:
    system: str
    observable: str
    n: int
    oracle: list[MStarCI] = field(default_factory=list)
    deployable: list[MStarCI] = field(default_factory=list)


def delta_sweep(
    system_name: str,
    observable_name: str,
    n: int,
    deltas: list[float],
    O_infinity: float,
    rng: np.random.Generator,
    *,
    R: int = 4000,
    alpha: float = 0.05,
    target: float = 0.80,
    m_min: int = 8,
    m_max: int | None = None,
    n_boot: int = 400,
) -> DeltaSweep:
    """Sweep ``delta``; for each, report ``m*`` (+bootstrap CI) under BOTH the oracle
    (calibrated) and deployable (Gaussian) nulls.  The equilibrium (H0) ensemble is shared
    across deltas; only H1 is regenerated per delta.
    """
    if m_max is None:
        m_max = n // 2
    m_grid = geometric_grid(m_min, m_max)
    D0 = _make_D(system_name, observable_name, n, R, 1.0, rng)
    stat0 = {m: _stat_at_m(D0, m, O_infinity, rng) for m in m_grid}

    out = DeltaSweep(system=system_name, observable=observable_name, n=n)
    for d in deltas:
        D1 = _make_D(system_name, observable_name, n, R, 1.0 - d, rng)
        stat1 = {m: _stat_at_m(D1, m, O_infinity, rng) for m in m_grid}
        for calibrate, bucket in ((True, out.oracle), (False, out.deployable)):
            bucket.append(mstar_ci_for_delta(
                stat0, stat1, m_grid, d, alpha=alpha, target=target,
                calibrate=calibrate, n_boot=n_boot, rng=rng))
    return out
