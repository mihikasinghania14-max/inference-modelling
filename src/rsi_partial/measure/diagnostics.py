"""Per-configuration diagnostic bundle (Phase-4 sec.8.6; architecture sec.2.5).

Every quantity here is MEASURED from data -- raw moments, variances, autocorrelation.
No theoretical prediction is referenced; scoring against P1-P8 happens only in
:mod:`rsi_partial.analysis`.  Functions take plain arrays/pools so they stay decoupled
from systems and observables.

Contents
--------
* Floor / susceptibility: ``Var_nu(O)``, ``chi_n``, ``R_n`` from a full-system ``O(x)`` pool.
* SNR: ``rho2_measured`` from ``(E[O])^2 / E_x[V_x]``; per-mixing-level ``rho2(phi)``.
* Sampling-distribution shape: skewness, excess kurtosis (Gaussian-approx flags).
* Effective information ``I_eff(m)=1/CV^2`` and crossover ``m_x``.
* Autocorrelation: integrated time ``tau_int`` (windowed ACF) and ``R_eff``.
* Self-calibration: ``rho2_hat``, ``R_n_hat`` from a single observed subsample.
* Equilibration check: ordered-vs-random agreement + running-trace stationarity.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import kurtosis, skew

# ---------------------------------------------------------------------------
# Floor, susceptibility, SNR
# ---------------------------------------------------------------------------


def floor_variance(O_full_pool: np.ndarray) -> float:
    """``Var_nu(O)`` -- the m-independent floor -- from full-system population values."""
    return float(np.var(np.asarray(O_full_pool, dtype=np.float64), ddof=1))


def susceptibility(O_full_pool: np.ndarray, n: int) -> float:
    """``chi_n[O] = n * Var_nu(O)`` (per-component susceptibility)."""
    return n * floor_variance(O_full_pool)


def relative_floor(O_full_pool: np.ndarray, O_infinity: float) -> float:
    """``R_n = Var_nu(O)/O_inf^2``.  ``nan`` if ``O_inf`` ~ 0 (degenerate control)."""
    if not np.isfinite(O_infinity) or abs(O_infinity) < 1e-300:
        return float("nan")
    return floor_variance(O_full_pool) / (O_infinity * O_infinity)


def rho2_measured(O_full_pool: np.ndarray, Vx_pool: np.ndarray) -> float:
    """Measured per-component SNR ``rho^2(phi) = (E[O])^2 / E_x[V_x]``.

    ``O_full_pool``: population values ``O(x)`` over replicate states at this mixing level.
    ``Vx_pool``:     within-config population variances ``V_x`` over the same states.
    """
    EO = float(np.mean(O_full_pool))
    EVx = float(np.mean(Vx_pool))
    if EVx <= 0.0:
        return float("nan")
    return (EO * EO) / EVx


# ---------------------------------------------------------------------------
# Sampling-distribution shape (Gaussian-approximation control sec.5.7)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ShapeFlags:
    skewness: float
    excess_kurtosis: float
    gaussian_ok: bool          # |skew|<=0.3 and |exkurt|<=0.5 (Phase-4 sec.5.7)


def sampling_shape(Ohat_pool: np.ndarray, *, skew_tol=0.3, kurt_tol=0.5) -> ShapeFlags:
    """Skewness and excess kurtosis of the estimator's sampling distribution, with a
    Gaussian-adequacy flag.  Bias-corrected (Fisher) moments."""
    x = np.asarray(Ohat_pool, dtype=np.float64)
    sk = float(skew(x, bias=False))
    ku = float(kurtosis(x, fisher=True, bias=False))
    return ShapeFlags(sk, ku, abs(sk) <= skew_tol and abs(ku) <= kurt_tol)


# ---------------------------------------------------------------------------
# Effective information and crossover
# ---------------------------------------------------------------------------


def cv_squared(Ohat_pool: np.ndarray, O_infinity: float) -> float:
    """``CV^2(Ohat_m) = Var(Ohat_m)/O_inf^2`` from a replicate pool."""
    if not np.isfinite(O_infinity) or abs(O_infinity) < 1e-300:
        return float("nan")
    return float(np.var(Ohat_pool, ddof=1)) / (O_infinity * O_infinity)


def effective_information(Ohat_pool: np.ndarray, O_infinity: float) -> float:
    """``I_eff(m) = 1/CV^2`` -- rises ~ rho^2 m then saturates at ``1/R_n`` (P8)."""
    c2 = cv_squared(Ohat_pool, O_infinity)
    if not np.isfinite(c2) or c2 <= 0.0:
        return float("nan")
    return 1.0 / c2


def crossover_m(rho2: float, R_n: float) -> float:
    """Information-saturation crossover ``m_x = 1/(rho^2 R_n)`` (P8)."""
    if rho2 <= 0.0 or R_n <= 0.0:
        return float("inf")
    return 1.0 / (rho2 * R_n)


# ---------------------------------------------------------------------------
# Autocorrelation (locally-generated chains only)
# ---------------------------------------------------------------------------


def integrated_autocorr_time(series: np.ndarray, *, c: float = 5.0) -> float:
    """Integrated autocorrelation time ``tau_int`` via the automatic windowing of
    Sokal: ``tau_int = 1 + 2 sum_{k=1}^{W} rho_k``, window ``W`` the smallest index
    with ``W >= c * tau_int(W)``.

    Returns ``0.5`` for an (anti)uncorrelated series by convention (so ``R_eff ~ R``).
    """
    x = np.asarray(series, dtype=np.float64)
    x = x - x.mean()
    N = x.shape[0]
    if N < 4 or np.allclose(x, 0.0):
        return 0.5
    # ACF via FFT
    f = np.fft.fft(x, n=2 * N)
    acf = np.fft.ifft(f * np.conjugate(f))[:N].real
    acf /= acf[0]
    tau = 1.0
    for W in range(1, N):
        tau = 1.0 + 2.0 * np.sum(acf[1:W + 1])
        if W >= c * tau:
            break
    return float(max(0.5, tau))


def effective_sample_size(R: int, tau_int: float, spacing: float = 1.0) -> float:
    """``R_eff = R / (1 + 2 tau_int / spacing)`` (Phase-4 sec.4.3)."""
    return R / (1.0 + 2.0 * tau_int / spacing)


# ---------------------------------------------------------------------------
# Self-calibration (Phase-4 sec.8.7): recover rho^2, R_n from the observed subsample
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SelfCalibration:
    rho2_hat: float
    R_n_hat: float
    mstar_hat: float           # predicted detection m* at the supplied (delta, kappa^2)


def self_calibrate(
    obs_values: np.ndarray,
    O_infinity: float,
    n: int,
    *,
    delta: float = 0.10,
    kappa_sq: float = 6.182557232019766,
) -> SelfCalibration:
    """Estimate ``rho^2`` and ``R_n`` from a SINGLE observed subsample of ``m`` values,
    then predict ``m*`` -- the deployable, system-agnostic path.

    ``rho2_hat = O_inf^2 / Var_hat(D)`` (uses the known equilibrium mean and the
    subsample's within-sample variance).  ``R_n_hat`` needs a floor estimate, which a
    single subsample cannot supply directly; here we return ``nan`` for ``R_n_hat`` (it
    requires either multiple states or the identity's plateau) and predict ``m*`` in the
    size-free approximation ``kappa^2/(rho2_hat delta^2)``.  Multi-state self-calibration
    of ``R_n`` is exercised in Stage 4.
    """
    v = np.asarray(obs_values, dtype=np.float64)
    Vd = float(np.var(v, ddof=1)) if v.shape[0] >= 2 else float("nan")
    if not np.isfinite(O_infinity) or abs(O_infinity) < 1e-300 or not np.isfinite(Vd) or Vd <= 0:
        return SelfCalibration(float("nan"), float("nan"), float("nan"))
    rho2_hat = (O_infinity * O_infinity) / Vd
    mstar_hat = kappa_sq / (rho2_hat * delta * delta)
    return SelfCalibration(rho2_hat=rho2_hat, R_n_hat=float("nan"), mstar_hat=mstar_hat)


def self_calibrate_multi(
    Ohat_pools: dict[int, np.ndarray],
    O_infinity: float,
    *,
    delta: float = 0.10,
    kappa_sq: float = 6.182557232019766,
) -> SelfCalibration:
    """Recover BOTH ``rho^2`` and the floor ``R_n`` from repeated subsamples, by fitting the
    variance identity to the measured coefficient of variation (critique sec.2.4):

        ``CV^2(m) = 1/(rho^2 m) + R_n``   (linear in ``1/m``: slope ``1/rho^2``, intercept ``R_n``).

    ``Ohat_pools`` maps at least two sample sizes ``m`` to a replicate pool of ``Ohat_m``
    (independent observations at that ``m``).  This is the deployable, system-agnostic
    self-calibration: it needs no knowledge of the system, only the estimator's own scatter
    across repeats -- and unlike the single-subsample path it recovers the floor, which is
    what tells you whether certification is even feasible (``delta <= kappa sqrt(R_n)``).

    ``mstar_hat`` uses the floor-corrected prediction
    ``kappa^2/[rho2_hat (delta^2 - kappa^2 R_n_hat)]`` (``inf`` if infeasible).
    """
    if not np.isfinite(O_infinity) or abs(O_infinity) < 1e-300:
        return SelfCalibration(float("nan"), float("nan"), float("nan"))
    ms = sorted(m for m in Ohat_pools if Ohat_pools[m].shape[0] >= 2)
    if len(ms) < 2:
        return SelfCalibration(float("nan"), float("nan"), float("nan"))
    x = np.array([1.0 / m for m in ms])
    y = np.array([float(np.var(Ohat_pools[m], ddof=1)) / (O_infinity ** 2) for m in ms])
    slope, intercept = np.polyfit(x, y, 1)          # CV^2 = slope*(1/m) + intercept
    rho2_hat = 1.0 / slope if slope > 0 else float("nan")
    R_n_hat = max(intercept, 0.0)                   # floor is non-negative
    denom = rho2_hat * (delta * delta - kappa_sq * R_n_hat) if np.isfinite(rho2_hat) else 0.0
    mstar_hat = kappa_sq / denom if denom > 0 else float("inf")
    return SelfCalibration(rho2_hat=rho2_hat, R_n_hat=R_n_hat, mstar_hat=mstar_hat)


# ---------------------------------------------------------------------------
# Equilibration check (Phase-4 sec.4.3): flag, never silently use
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class EquilibrationCheck:
    ordered_mean: float
    random_mean: float
    agree: bool                # observable from the two initial conditions agrees within tol
    trace_stationary: bool     # running-trace second-half vs first-half agree within tol
    passed: bool


def equilibration_check(
    obs_from_ordered: np.ndarray,
    obs_from_random: np.ndarray,
    *,
    tol_sigma: float = 3.0,
) -> EquilibrationCheck:
    """Two-initial-condition agreement + running-trace stationarity.

    ``obs_from_ordered`` / ``obs_from_random``: pools of the observable's population value
    generated (after the same burn-in) from ordered vs random starts.  They must agree
    within ``tol_sigma`` standard errors, and each trace's two halves must agree.
    """
    a = np.asarray(obs_from_ordered, dtype=np.float64)
    b = np.asarray(obs_from_random, dtype=np.float64)
    ma, mb = float(a.mean()), float(b.mean())
    sea = a.std(ddof=1) / np.sqrt(a.shape[0]) if a.shape[0] > 1 else np.inf
    seb = b.std(ddof=1) / np.sqrt(b.shape[0]) if b.shape[0] > 1 else np.inf
    se = np.sqrt(sea ** 2 + seb ** 2)
    agree = bool(abs(ma - mb) <= tol_sigma * se) if np.isfinite(se) and se > 0 else False

    def _stationary(x: np.ndarray) -> bool:
        h = x.shape[0] // 2
        if h < 2:
            return True
        f, s = x[:h], x[h:]
        d = abs(f.mean() - s.mean())
        sef = np.sqrt(f.var(ddof=1) / f.shape[0] + s.var(ddof=1) / s.shape[0])
        return bool(d <= tol_sigma * sef) if sef > 0 else True

    trace_stationary = _stationary(a) and _stationary(b)
    return EquilibrationCheck(
        ordered_mean=ma, random_mean=mb, agree=agree,
        trace_stationary=trace_stationary, passed=agree and trace_stationary,
    )
