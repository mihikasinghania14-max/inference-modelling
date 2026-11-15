"""Analytical ground-truth tests -- assert against CLOSED FORMS only.

No self-referential fixtures: every assertion targets a formula from the RSI paper or the
exact product-chain algebra, not previously-saved simulator output.  Monte-Carlo
assertions use a fixed seed and a generous multiple of the standard error so they are
reproducible and non-flaky.
"""

from __future__ import annotations

import numpy as np
import pytest

from rsi_partial.analysis import predictions as P
from rsi_partial.experiments.engine import measure_ensemble
from rsi_partial.inference.estimator import estimate, phi_hat
from rsi_partial.observables.cards_ad import AverageDisplacement
from rsi_partial.observables.chain_mag import ChainMagnetization
from rsi_partial.observe.veil import reveal
from rsi_partial.systems.cards import Cards
from rsi_partial.systems.product_chain import ProductChain

SEED = 20260818


# ---------------------------------------------------------------------------
# Cards: equilibrium displacement, mean law, SNR
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", [128, 512, 3000])
def test_cards_AD_infinity(n):
    """MC mean of AD over equilibrium shuffles -> (n^2-1)/(3n)  (Prop. 2.6)."""
    rng = np.random.default_rng(SEED + n)
    R = 8000
    ad = AverageDisplacement()
    cards = Cards()
    vals = np.array([ad.population_value(cards.sample_equilibrium(n, rng)) for _ in range(R)])
    theory = P.Cards.AD_infinity(n)
    se = vals.std(ddof=1) / np.sqrt(R)
    assert abs(vals.mean() - theory) < 5 * se


def test_cards_mean_law():
    """E[AD(t)] = (1 - (1-2/n)^t) AD_inf across t (single-card kernel)."""
    n, R = 256, 4000
    rng = np.random.default_rng(SEED)
    cards = Cards()
    ad = AverageDisplacement()
    AD_inf = P.Cards.AD_infinity(n)
    for t in [n // 4, n, 2 * n]:
        # vectorised relaxation of R decks by exactly t transpositions
        Pm = np.tile(np.arange(n, dtype=np.int64), (R, 1))
        rows = np.arange(R)
        for _ in range(t):
            a = rng.integers(n, size=R); b = rng.integers(n, size=R)
            Pa = Pm[rows, a].copy(); Pm[rows, a] = Pm[rows, b]; Pm[rows, b] = Pa
        AD_t = np.abs(Pm - np.arange(n)).mean(axis=1)
        theory = P.Cards.phi_of_t(t, n) * AD_inf
        se = AD_t.std(ddof=1) / np.sqrt(R)
        assert abs(AD_t.mean() - theory) < 5 * se, (t, AD_t.mean(), theory)


def test_cards_rho2_equilibrium():
    """Measured per-card SNR at equilibrium -> 2 (RSI sec.6)."""
    n, R = 512, 8000
    rng = np.random.default_rng(SEED + 1)
    res = measure_ensemble(Cards(), AverageDisplacement(), n, R, rng, [n], phi=1.0)
    from rsi_partial.measure.diagnostics import rho2_measured
    rho2 = rho2_measured(res.O_full, res.Vx)
    assert abs(rho2 - 2.0) < 0.05


# ---------------------------------------------------------------------------
# Unbiasedness (Lemma 2.10) and the census identity
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n,m", [(512, 32), (512, 128), (1024, 64)])
def test_unbiased_fixed_state(n, m):
    """For a FIXED state, E_S[Ohat_m] = O(x) (unbiasedness of the sample mean)."""
    rng = np.random.default_rng(SEED + n + m)
    cards = Cards()
    ad = AverageDisplacement()
    state = cards.sample_equilibrium(n, rng)
    D = ad.per_component(state)
    Ox = float(D.mean())
    R = 8000
    Ohat = np.array([estimate(reveal(D, m, rng)).Ohat for _ in range(R)])
    se = Ohat.std(ddof=1) / np.sqrt(R)
    assert abs(Ohat.mean() - Ox) < 5 * se


@pytest.mark.parametrize("n", [64, 512])
def test_census_is_exact(n):
    """Ohat_n = O(x) IDENTICALLY (not merely in expectation)."""
    rng = np.random.default_rng(SEED + n)
    cards = Cards()
    ad = AverageDisplacement()
    for _ in range(20):
        state = cards.sample_equilibrium(n, rng)
        D = ad.per_component(state)
        obs = reveal(D, n, rng)
        assert np.isclose(estimate(obs).Ohat, float(D.mean()), rtol=0, atol=1e-9)
        assert obs.reads == n


# ---------------------------------------------------------------------------
# Exact estimator variance (RSI App. E.1/E.3)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("m", [32, 64, 128])
def test_cards_exact_variance(m):
    """Var(Ohat_m) matches (n-m)/(m(n-1)) E[V_x] + Var(AD)  (identity star, cards)."""
    n, R = 512, 12000
    rng = np.random.default_rng(SEED + m)
    res = measure_ensemble(Cards(), AverageDisplacement(), n, R, rng, [m], phi=1.0)
    measured = float(np.var(res.pools[m].Ohat, ddof=1))
    EVx = float(np.mean(res.Vx))
    floor = float(np.var(res.O_full, ddof=1))
    exact = P.cards_exact_var_Ohat(n, m, EVx, floor)
    # MC relative SE of a variance ~ sqrt(2/R) ~ 1.3% at R=12000; allow 6%.
    assert abs(measured - exact) / exact < 0.06, (measured, exact)


# ---------------------------------------------------------------------------
# Product chain -- EXACT laws
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("phi", [0.3, 0.6, 1.0])
def test_chain_rho2_exact(phi):
    """rho^2(phi) = phi/(2-phi) exactly (independent Bernoulli(phi/2) coordinates)."""
    n, R = 4096, 8000
    rng = np.random.default_rng(SEED + int(100 * phi))
    res = measure_ensemble(ProductChain(), ChainMagnetization(), n, R, rng, [n], phi=phi)
    from rsi_partial.measure.diagnostics import rho2_measured
    rho2 = rho2_measured(res.O_full, res.Vx)
    assert abs(rho2 - P.Chain.rho2(phi)) < 0.03


@pytest.mark.parametrize("phi,m", [(0.6, 64), (1.0, 128)])
def test_chain_var_phihat_exact(phi, m):
    """Var(phi_hat_m) = phi(2-phi)/m exactly (phi_hat = Ohat/O_inf, O_inf=1/2)."""
    n, R = 8192, 12000
    rng = np.random.default_rng(SEED + m + int(10 * phi))
    res = measure_ensemble(ProductChain(), ChainMagnetization(), n, R, rng, [m], phi=phi)
    phih = res.pools[m].Ohat / 0.5
    measured = float(np.var(phih, ddof=1))
    theory = P.Chain.var_phihat(phi, m)
    assert abs(measured - theory) / theory < 0.06, (measured, theory)


def test_phi_hat_degenerate_guard():
    """phi_hat returns nan (not a blow-up) when O_inf -> 0 (S2' degenerate control)."""
    assert np.isnan(phi_hat(0.01, 0.0))
    assert np.isclose(phi_hat(1.0, 0.5), 2.0)
