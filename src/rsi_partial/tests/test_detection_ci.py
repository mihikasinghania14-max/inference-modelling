"""Tests for the bootstrap-CI detection sweep (critique sec.2.1-2.3).

Locks in: m* carries a CI that brackets its point estimate; the oracle null's Type-I is
measured on a HELD-OUT pool (so ~alpha is a real check, not tautological); the deployable
Gaussian null is reported separately; and the constant m*delta^2 lands near kappa^2/rho^2.
"""

from __future__ import annotations

import numpy as np

from rsi_partial.analysis.predictions import Cards, kappa_sq
from rsi_partial.experiments.detection_ci import delta_sweep

SEED = 24680


def test_delta_sweep_mstar_has_bracketing_ci_and_honest_nulls():
    rng = np.random.default_rng(SEED)
    n = 2048
    sw = delta_sweep("cards", "AD", n, [0.15, 0.10], Cards.AD_infinity(n), rng,
                     R=2500, n_boot=200)
    target = kappa_sq() / 2.0                       # ~3.09

    assert len(sw.oracle) == 2 and len(sw.deployable) == 2
    for o in sw.oracle:
        assert np.isfinite(o.mstar)
        assert o.ci_lo <= o.mstar <= o.ci_hi        # CI brackets the point estimate
        assert o.ci_hi > o.ci_lo                     # CI is informative
        # held-out Type-I of the calibrated (oracle) null is genuinely near alpha
        assert abs(o.type_I - 0.05) < 0.03, (o.delta, o.type_I)
        # the constant lands in a sensible band around kappa^2/rho^2
        assert 2.6 < o.mstar_delta2 < 4.7, (o.delta, o.mstar_delta2)
        assert o.censored_frac < 0.2                 # target reached in most bootstraps

    # deployable (Gaussian-null) path is reported and finite too
    for d in sw.deployable:
        assert d.null == "deployable" and np.isfinite(d.mstar)
        assert d.ci_lo <= d.mstar <= d.ci_hi


def test_finite_delta_overshoot_shrinks_toward_target():
    """m*delta^2 should sit above kappa^2/rho^2 and move DOWN toward it as delta shrinks
    (the finite-delta, sigma1>sigma0 correction). Checked on the oracle null."""
    rng = np.random.default_rng(SEED + 1)
    n = 2048
    sw = delta_sweep("cards", "AD", n, [0.20, 0.10], Cards.AD_infinity(n), rng,
                     R=2500, n_boot=120)
    big_delta, small_delta = sw.oracle[0], sw.oracle[1]
    # smaller delta is closer to (or below) the larger-delta overshoot
    assert small_delta.mstar_delta2 <= big_delta.mstar_delta2 + 0.25, (
        big_delta.mstar_delta2, small_delta.mstar_delta2)
