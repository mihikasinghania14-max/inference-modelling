"""Closed-form quantities for the random transposition shuffle at equilibrium.

Every function states where its formula comes from.  Functions marked *verified* are
checked against exact enumeration of the symmetric group in
``tests/test_rsi_shuffle_theory.py``.

Notation: ``n`` cards, ``pi`` uniform on ``S_n``, ``D_i = |pi(i) - i|``,
``AD = (1/n) sum_i D_i``, ``F = n * AD`` (Spearman's footrule).
"""

from __future__ import annotations

import math
from fractions import Fraction

from scipy.stats import norm

Number = int | float | Fraction


def _F(x: Number) -> Fraction:
    return x if isinstance(x, Fraction) else Fraction(x)


def ad_inf(n: int, exact: bool = False) -> Fraction | float:
    """``AD_inf = (n^2 - 1) / (3 n)``  (RSI Proposition 2.6; *verified*)."""
    v = Fraction(n * n - 1, 3 * n)
    return v if exact else float(v)


def second_moment_displacement(n: int, exact: bool = False) -> Fraction | float:
    """``E[D^2] = (n^2 - 1)/6`` for one card at equilibrium (RSI Section 5.4)."""
    v = Fraction(n * n - 1, 6)
    return v if exact else float(v)


def var_displacement_single_card(n: int, exact: bool = False) -> Fraction | float:
    """``Var(D) = (n^2-1)/6 - AD_inf^2 = n^2/18 (1 + O(n^-2))`` (RSI Theorem 5.1)."""
    v = second_moment_displacement(n, exact=True) - ad_inf(n, exact=True) ** 2
    return v if exact else float(v)


def var_footrule(n: int, exact: bool = False) -> Fraction | float:
    """``Var(F) = (n+1)(2n^2+7)/45``  (Diaconis & Graham 1977; used by the extension,
    Section 4; *verified* by enumeration for n <= 7 in the test suite).

    The extension paper attributes this to "Singhania 2026, App. E"; RSI Appendix E
    only states ``Var = O(n^3)``.  See ``docs/known-issues.md`` (item C1).
    """
    v = Fraction((n + 1) * (2 * n * n + 7), 45)
    return v if exact else float(v)


def var_ad(n: int, exact: bool = False) -> Fraction | float:
    """``Var_nu(AD) = Var(F)/n^2``."""
    v = var_footrule(n, exact=True) / (n * n)
    return v if exact else float(v)


def relative_floor(n: int, exact: bool = False) -> Fraction | float:
    """``R_n = Var_nu(AD)/AD_inf^2 = (2n^2+7) / (5 (n-1)^2 (n+1))``  (Extension, Eq. (2);
    *verified*).  ``n R_n -> 2/5``."""
    if n < 2:
        raise ValueError("n >= 2 required")
    v = Fraction(2 * n * n + 7, 5 * (n - 1) ** 2 * (n + 1))
    return v if exact else float(v)


def expected_within_deck_variance(n: int, exact: bool = False) -> Fraction | float:
    """``E_nu[V_x]``, the mean within-configuration population variance of the
    displacements, ``(n-2)(n+1)(5n^2+n+12) / (90 n^2)``.

    Follows from RSI Eq. (E.4) with ``Var(F)`` above; *verified* by enumeration (n <= 7).
    Leading order ``n^2/18``.
    """
    v = Fraction((n - 2) * (n + 1) * (5 * n * n + n + 12), 90 * n * n)
    return v if exact else float(v)


def estimator_variance(n: int, m: int, exact: bool = False) -> Fraction | float:
    """Exact variance of the uniform-sample estimator at equilibrium (Lemma 1 of the
    extension, equivalently RSI Eq. (E.3))::

        Var(AD_hat_m) = (n - m)/(m (n - 1)) * E[V_x] + Var_nu(AD).

    *Verified* by enumeration (n <= 6).  RSI Theorem 5.1 is its leading order,
    ``n^2/(18 m)``, for ``m << n``.
    """
    if not 1 <= m <= n or n < 2:
        raise ValueError("need n >= 2 and 1 <= m <= n")
    v = Fraction(n - m, m * (n - 1)) * expected_within_deck_variance(n, exact=True) \
        + var_ad(n, exact=True)
    return v if exact else float(v)


def exact_cv2(n: int, m: int) -> float:
    """Squared coefficient of variation ``Var(AD_hat_m)/AD_inf^2`` (exact)."""
    return float(estimator_variance(n, m, exact=True) / ad_inf(n, exact=True) ** 2)


def mean_law_phi(t: int | float, n: int) -> float:
    """``phi(t) = E[AD(pi_t)] / AD_inf = 1 - (1 - 2/n)^t`` for the independent-pair chain
    started from the identity (extension v2 Proposition 2; *verified* for n = 4 by exact
    transition-matrix computation).  Not valid for the distinct-pair convention."""
    return 1.0 - (1.0 - 2.0 / n) ** t


def m_star_rsi(eps: float, alpha: float = 0.05) -> float:
    """RSI Eq. (5.7): ``m*(eps) = z_{1-alpha/2}^2 / (2 eps^2)``, the sample size for which
    the ``1-alpha`` credible half-width is ``eps * AD_inf`` (leading order, ``m* << n``).
    About 48, 192 and 768 for ``eps`` = 0.20, 0.10, 0.05 at ``alpha = 0.05``."""
    z = norm.ppf(1.0 - alpha / 2.0)
    return z * z / (2.0 * eps * eps)


def kappa(alpha: float = 0.05, beta: float = 0.20) -> float:
    """``kappa = z_{1-alpha} + z_{1-beta}`` (extension, Corollary 1); 2.486 for the
    defaults (one-sided test of size 0.05 and power 0.80)."""
    return float(norm.ppf(1.0 - alpha) + norm.ppf(1.0 - beta))


def certification_cost(delta: float, n: int, alpha: float = 0.05, beta: float = 0.20,
                       rho2: float = 2.0) -> float:
    """Extension Section 4: ``m*(delta) = kappa^2 / (rho^2 (delta^2 - kappa^2 R_n))``.

    Returns ``inf`` when ``delta <= kappa sqrt(R_n)`` (certification impossible at any
    ``m``).  Leading-order (``m << n``); the finite-population correction is not applied.
    """
    k2 = kappa(alpha, beta) ** 2
    gap = delta * delta - k2 * relative_floor(n)
    return math.inf if gap <= 0 else k2 / (rho2 * gap)
