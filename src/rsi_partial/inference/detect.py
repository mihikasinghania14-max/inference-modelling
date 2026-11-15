"""Detection / power module (first-class, per design brief sec.4c and Phase-4 sec.8.5).

Certification test:  ``H0: phi = 1``  vs  ``H1: 1 - phi >= delta``  (one-sided; the
relaxation approaches equilibrium from below, ``phi <= 1``).

Decision rule for one observation, given the estimator ``Ohat_m``, an estimate of its
sampling standard deviation ``sd``, and the known equilibrium mean ``O_inf``:

    reject H0  <=>  T := (O_inf - Ohat_m) / sd  >  z_{1-alpha}.

``T`` is standard normal under H0 (Type-I = alpha); under ``H1: phi = 1-delta`` its mean
is ``delta / CV``, giving power ``Phi(delta/CV - z_{1-alpha})`` -- the Gaussian basis of
``m*(delta) = kappa^2/(rho^2 delta^2)``.  Power and Type-I are estimated empirically from
replicate pools; a bootstrap/exact null calibration is available for small ``m`` where the
CLT threshold is mis-sized (control Phase-4 sec.5.7).

``sd`` is the estimator's own MEASURED sampling SD -- either per-replicate
(self-calibrated / deployable) or the pooled H0 SD (oracle) -- never a theory value.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import norm


def reject_H0(
    Ohat: np.ndarray,
    sd: np.ndarray,
    O_infinity: float,
    *,
    alpha: float = 0.05,
    threshold: float | None = None,
) -> np.ndarray:
    """Vectorised one-sided decision: ``True`` where H0 (equilibrium) is rejected.

    ``threshold`` overrides the Gaussian critical value ``z_{1-alpha}`` (used by the
    calibrated small-``m`` path).
    """
    zc = norm.ppf(1.0 - alpha) if threshold is None else threshold
    Ohat = np.asarray(Ohat, dtype=np.float64)
    sd = np.asarray(sd, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.where(sd > 0.0, (O_infinity - Ohat) / sd, 0.0)
    return T > zc


@dataclass(frozen=True)
class DetectionResult:
    power: float
    type_I: float
    n_H0: int
    n_H1: int
    alpha: float
    threshold: float
    calibrated: bool


def power_and_type_I(
    Ohat_H0: np.ndarray,
    sd_H0: np.ndarray,
    Ohat_H1: np.ndarray,
    sd_H1: np.ndarray,
    O_infinity: float,
    *,
    alpha: float = 0.05,
    calibrate: bool = False,
) -> DetectionResult:
    """Estimate power and Type-I from equilibrium (H0) and deviated (H1) replicate pools.

    Parameters
    ----------
    Ohat_H0, sd_H0: estimator values and sampling-SD estimates over H0 (equilibrium)
        replicates -- give the empirical Type-I.
    Ohat_H1, sd_H1: same over H1 (relaxed to ``phi = 1-delta``) replicates -- give power.
    calibrate: if True, set the rejection threshold from the empirical H0 statistic so
        that Type-I == alpha exactly (bootstrap/exact null calibration for small ``m``),
        then measure power at that threshold.
    """
    zc = norm.ppf(1.0 - alpha)
    calibrated = False
    if calibrate:
        with np.errstate(divide="ignore", invalid="ignore"):
            T0 = np.where(sd_H0 > 0.0, (O_infinity - Ohat_H0) / sd_H0, 0.0)
        # threshold = upper-alpha empirical quantile of the H0 statistic
        zc = float(np.quantile(T0, 1.0 - alpha))
        calibrated = True
    rej_H0 = reject_H0(Ohat_H0, sd_H0, O_infinity, alpha=alpha, threshold=zc)
    rej_H1 = reject_H0(Ohat_H1, sd_H1, O_infinity, alpha=alpha, threshold=zc)
    return DetectionResult(
        power=float(np.mean(rej_H1)),
        type_I=float(np.mean(rej_H0)),
        n_H0=int(Ohat_H0.shape[0]),
        n_H1=int(Ohat_H1.shape[0]),
        alpha=alpha,
        threshold=float(zc),
        calibrated=calibrated,
    )
