"""Detection-``m*`` driver: certify a deviation ``delta`` for a (system, observable, n).

Generates an equilibrium (H0, ``phi=1``) and a relaxed (H1, ``phi=1-delta``) ensemble of
per-component ``D`` arrays, then locates the smallest ``m`` whose power reaches the target
by bisection (reusing the held arrays, so bisection is cheap).  The H1 ensemble is made
with the **physical dynamics** (the process whose convergence is certified), never a
theory shortcut.

Memory note: this driver holds the two ``R x n`` ``D`` arrays so it can subsample at any
``m`` during bisection.  That is fine at validation scale; the large-``n`` Phase-5 configs
should instead sweep a fixed geometric ``m``-grid in batches (engine.measure_ensemble) and
interpolate ``m*``.  This module is the accurate, moderate-``n`` path.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..inference.detect import power_and_type_I
from ..inference.mstar import MStarResult, detection_mstar
from ..observe.veil import reveal_batch
from ..observables.cards_ad import AverageDisplacement
from ..observables.cards_fixedpoint import FixedPointFraction
from ..observables.chain_mag import ChainMagnetization
from ..systems.cards import Cards
from ..systems.product_chain import ProductChain


def _cards_D(n, R, phi, rng, observable_name="AD"):
    if phi >= 1.0:
        P = np.argsort(rng.random((R, n)), axis=1)
    else:
        t = Cards().steps_for_phi(n, phi)
        P = np.tile(np.arange(n, dtype=np.int64), (R, 1))
        rows = np.arange(R)
        for _ in range(t):
            a = rng.integers(n, size=R); c = rng.integers(n, size=R)
            Pa = P[rows, a].copy(); P[rows, a] = P[rows, c]; P[rows, c] = Pa
    idx = np.arange(n)
    if observable_name == "fixed_point_fraction":
        return (P == idx).astype(np.float64)
    return np.abs(P - idx).astype(np.float64)


def _chain_D(n, R, phi, rng):
    return (rng.random((R, n)) < 0.5 * phi).astype(np.float64)


def _make_D(system_name, observable_name, n, R, phi, rng):
    if system_name == "cards":
        return _cards_D(n, R, phi, rng, observable_name)
    if system_name == "product_chain":
        return _chain_D(n, R, phi, rng)
    raise ValueError(f"detection driver has no fast path for system {system_name!r}")


def _subsample_stats(D, m, rng):
    R, n = D.shape
    # Selection goes through the veil (independent draw per m: no perm reuse here, so the
    # H0/H1 pools at each m are independent), so the audited path is the executed path.
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
    return Ohat, sd


@dataclass
class DetectionRun:
    system: str
    observable: str
    n: int
    delta: float
    O_infinity: float
    mstar: int | None
    mstar_delta2: float
    result: MStarResult
    type_I_at_mstar: float
    power_at_mstar: float
    calibrated: bool


def detection_experiment(
    system_name: str,
    observable_name: str,
    n: int,
    delta: float,
    O_infinity: float,
    rng: np.random.Generator,
    *,
    R: int = 4000,
    alpha: float = 0.05,
    power_target: float = 0.80,
    m_min: int = 8,
    m_max: int | None = None,
    tol: float = 0.05,
    calibrate: bool = True,
) -> DetectionRun:
    """Locate detection ``m*(delta)`` for certifying ``phi=1`` against ``phi=1-delta``.

    ``calibrate=True`` uses the exact-null calibration (control sec.5.7) so Type-I == alpha
    exactly -- the method of record for the size-free constant ``kappa^2/rho^2``.
    """
    if m_max is None:
        m_max = n // 2
    D0 = _make_D(system_name, observable_name, n, R, 1.0, rng)
    D1 = _make_D(system_name, observable_name, n, R, 1.0 - delta, rng)

    def power_of_m(m: int) -> float:
        Oh0, sd0 = _subsample_stats(D0, m, rng)
        Oh1, sd1 = _subsample_stats(D1, m, rng)
        return power_and_type_I(Oh0, sd0, Oh1, sd1, O_infinity,
                                alpha=alpha, calibrate=calibrate).power

    res = detection_mstar(power_of_m, power_target=power_target,
                          m_min=m_min, m_max=m_max, tol=tol)
    ms = res.mstar
    if ms is not None:
        Oh0, sd0 = _subsample_stats(D0, ms, rng)
        Oh1, sd1 = _subsample_stats(D1, ms, rng)
        dr = power_and_type_I(Oh0, sd0, Oh1, sd1, O_infinity, alpha=alpha, calibrate=calibrate)
        typeI, power = dr.type_I, dr.power
        md2 = ms * delta * delta
    else:
        typeI, power, md2 = float("nan"), float("nan"), float("nan")
    return DetectionRun(
        system=system_name, observable=observable_name, n=n, delta=delta,
        O_infinity=O_infinity, mstar=ms, mstar_delta2=md2, result=res,
        type_I_at_mstar=typeI, power_at_mstar=power, calibrated=calibrate,
    )
