"""Conjugate Gaussian Bayesian inference for the latent global observable.

Exactly the RSI paper's model (sec.5.3 / App. F):

  * likelihood      ``Ohat_m | mu ~ N(mu, sigma_m^2)``
  * conjugate prior ``mu ~ N(mu0, sigma0^2)``  (weak: ``mu0 = O_inf``, large ``sigma0``)
  * posterior       ``mu | Ohat_m ~ N(mu_post, sigma_post^2)`` with additive precision

        1/sigma_post^2 = 1/sigma0^2 + 1/sigma_m^2
        mu_post        = sigma_post^2 (mu0/sigma0^2 + Ohat_m/sigma_m^2)

  * credible interval  ``mu_post +/- z_{1-alpha/2} sigma_post``.

``sigma_m^2`` is the **measured** sampling variance of the estimator (or, only as a
labelled cross-check, the analytic ``n^2/18m`` for cards) -- never a theory value fed
in as if it were data.  Nothing here imports :mod:`rsi_partial.analysis.predictions`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.stats import norm


@dataclass(frozen=True)
class Posterior:
    mu_post: float
    sigma_post: float
    ci_lo: float
    ci_hi: float
    alpha: float

    @property
    def half_width(self) -> float:
        return 0.5 * (self.ci_hi - self.ci_lo)


def conjugate_update(
    Ohat_m: float,
    sigma_m2: float,
    *,
    mu0: float,
    sigma0: float = 100.0,
    alpha: float = 0.05,
) -> Posterior:
    """Closed-form Gaussian posterior + central ``(1-alpha)`` credible interval.

    Parameters
    ----------
    Ohat_m:   observed sample mean (the estimator).
    sigma_m2: sampling variance of the estimator (``sigma_m^2``), MEASURED.
    mu0:      prior mean (weak prior centred at ``O_inf``).
    sigma0:   prior standard deviation (weak: ``sigma0 >> sigma_m``).
    alpha:    credible level is ``1-alpha`` (two-sided).
    """
    if sigma_m2 <= 0.0:
        # Degenerate (m too small to estimate spread); fall back to prior.
        prec_data = 0.0
    else:
        prec_data = 1.0 / sigma_m2
    prec_prior = 1.0 / (sigma0 * sigma0)
    prec_post = prec_prior + prec_data
    sigma_post2 = 1.0 / prec_post
    mu_post = sigma_post2 * (mu0 * prec_prior + Ohat_m * prec_data)
    sigma_post = math.sqrt(sigma_post2)
    zc = float(norm.ppf(1.0 - alpha / 2.0))
    return Posterior(
        mu_post=mu_post,
        sigma_post=sigma_post,
        ci_lo=mu_post - zc * sigma_post,
        ci_hi=mu_post + zc * sigma_post,
        alpha=alpha,
    )


def certifies_equilibrium(
    post: Posterior,
    O_infinity: float,
    epsilon: float,
) -> bool:
    """Estimation-face certification (RSI sec.6 / App. F.2): the posterior places its
    whole credible interval within a relative tolerance ``epsilon`` of ``O_inf``.

        certified  <=>  [ci_lo, ci_hi] subset of (1 +/- epsilon) * O_inf.
    """
    lo = (1.0 - epsilon) * O_infinity
    hi = (1.0 + epsilon) * O_infinity
    return post.ci_lo >= lo and post.ci_hi <= hi
