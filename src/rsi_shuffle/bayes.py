"""Conjugate Gaussian inference for the equilibrium displacement (RSI Section 5.3).

Model (RSI Eqs. 3.1, 5.1, 5.2 and Appendix F)::

    AD_hat_m | mu ~ N(mu, sigma_m^2),        sigma_m^2 = n^2 / (18 m)
    mu            ~ N(mu0, sigma0^2)

so that the posterior is ``N(mu_post, sigma_post^2)`` with precision-additive update.

Scope note
----------
``sigma_m^2 = n^2/(18 m)`` is the leading-order sampling variance for ``m << n``
(RSI Theorem 5.1).  It omits both the finite-population correction and the
equilibrium floor ``Var(AD)`` that the extension identifies; see
``docs/known-issues.md`` (items M1-M3).  The functions below implement the model
exactly as written in the paper; they do not alter it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.stats import norm


def sampling_sd(n: int, m: int) -> float:
    """Leading-order sd of the estimator, ``n / sqrt(18 m)`` (RSI Theorem 5.1)."""
    return n / math.sqrt(18.0 * m)


@dataclass(frozen=True)
class Posterior:
    mean: float
    sd: float

    def interval(self, alpha: float = 0.05) -> tuple[float, float]:
        z = norm.ppf(1.0 - alpha / 2.0)
        return self.mean - z * self.sd, self.mean + z * self.sd


def conjugate_posterior(y: float, sigma_m: float, mu0: float, sigma0: float) -> Posterior:
    """Posterior of ``mu`` given one observation ``y = AD_hat_m`` (RSI Eqs. 5.1-5.2).

    ``sigma0 = math.inf`` gives the flat-prior limit (Appendix F.1).
    """
    prec_prior = 0.0 if math.isinf(sigma0) else 1.0 / sigma0**2
    prec_data = 1.0 / sigma_m**2
    prec = prec_prior + prec_data
    mean = ((mu0 * prec_prior if prec_prior else 0.0) + y * prec_data) / prec
    return Posterior(mean=mean, sd=math.sqrt(1.0 / prec))
