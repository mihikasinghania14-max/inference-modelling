"""Frequentist estimator ``Ohat_m`` and the variance decomposition (identity star).

Everything here is a raw measurement.  No theoretical value is used.

The spine identity (Phase-4 sec.1.1; RSI App. E):

    Var(Ohat_m) = (1/m) (n-m)/(n-1) E_x[V_x]   +   Var_x(O)
                  \\_____ sampling term _____/       \\_ floor _/

:func:`estimate` gives the per-observation point estimate and its sampling-variance
estimate.  :func:`variance_split` measures both terms of the identity *independently*
from a replicate pool, so the identity can be checked rather than assumed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.types import Observation


@dataclass(frozen=True)
class PointEstimate:
    Ohat: float
    Vx_hat: float          # within-config population variance, estimated from the m values
    samp_var_raw: float    # Vx_hat / m
    samp_var_fpc: float    # Vx_hat / m * (n-m)/(n-1)   (finite-population correction)
    m: int
    n: int


def estimate(obs: Observation) -> PointEstimate:
    """Point estimate and sampling-variance estimate from a single observation.

    The sample variance of the ``m`` revealed values is an (unbiased, ddof=1) estimate
    of the within-config population variance ``V_x``; dividing by ``m`` and applying the
    finite-population correction ``(n-m)/(n-1)`` gives the sampling variance of the mean.
    """
    v = obs.values
    m, n = obs.m, obs.n
    Ohat = float(np.mean(v))
    if m >= 2:
        Vx_hat = float(np.var(v, ddof=1))
    else:
        Vx_hat = 0.0
    samp_var_raw = Vx_hat / m
    fpc = (n - m) / (n - 1) if n > 1 else 0.0
    samp_var_fpc = samp_var_raw * fpc
    return PointEstimate(
        Ohat=Ohat, Vx_hat=Vx_hat, samp_var_raw=samp_var_raw,
        samp_var_fpc=samp_var_fpc, m=m, n=n,
    )


def phi_hat(Ohat: float, O_infinity: float) -> float:
    """Relative mixing level estimate ``phi_hat = Ohat / O_inf``.

    Guarded for ``O_inf -> 0`` (the S2' degenerate control): returns ``nan`` if the
    equilibrium mean is (numerically) zero, so callers annotate rather than divide-blow-up.
    """
    if not np.isfinite(O_infinity) or abs(O_infinity) < 1e-300:
        return float("nan")
    return Ohat / O_infinity


@dataclass(frozen=True)
class VarianceSplit:
    var_Ohat: float          # measured Var(Ohat_m) across replicates
    sampling_term: float     # (1/m)(n-m)/(n-1) * mean_x[V_x]
    floor: float             # Var_x(O), measured from full-system population values
    predicted_sum: float     # sampling_term + floor
    m: int
    n: int


def variance_split(
    Ohat_pool: np.ndarray,
    Vx_pool: np.ndarray,
    O_full_pool: np.ndarray,
    m: int,
    n: int,
) -> VarianceSplit:
    """Measure both terms of identity (star) independently from a replicate pool.

    Parameters
    ----------
    Ohat_pool: estimator values ``Ohat_m`` over replicates (fresh state + fresh subset).
    Vx_pool:   within-config population variances ``V_x`` over the same states.
    O_full_pool: full-system population values ``O(x)`` over the same states.
    """
    var_Ohat = float(np.var(Ohat_pool, ddof=1))
    fpc = (n - m) / (n - 1) if n > 1 else 0.0
    sampling_term = (1.0 / m) * fpc * float(np.mean(Vx_pool))
    floor = float(np.var(O_full_pool, ddof=1))
    return VarianceSplit(
        var_Ohat=var_Ohat, sampling_term=sampling_term, floor=floor,
        predicted_sum=sampling_term + floor, m=m, n=n,
    )
