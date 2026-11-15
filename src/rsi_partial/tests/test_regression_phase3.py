"""Phase-3 regression suite -- THE GATE.

Reproduce the Phase-3 pilot findings within Monte-Carlo error before any Phase-5 run is
permitted (design brief; architecture sec.5.2).  Quantities with closed forms / robust pilot
values are asserted numerically; the near-``T_c`` Ising floor is asserted only
QUALITATIVELY (modelling note M4: the Phase-3 floor values used finite burn-in and are
labelled qualitatively-reliable; the new Wolff ground truth deliberately corrects that
bias, so reproducing the approximate numbers would signal a reintroduced bug).

Phase-3 findings reproduced here:
  1. off-equilibrium SNR laws rho^2(phi): chain phi/(2-phi) (exact), cards phi/(3/2-phi)
     (leading order);
  2. n-independence of the equilibrium SNR and of the detection m*;
  3. detection cost m* delta^2 -> kappa^2/rho^2  (cards ~3.1, chain ~6.3);
  4. Ising high-T subsample variance ~ 1/m; near-T_c an m-independent floor.
"""

from __future__ import annotations

import numpy as np
import pytest

from rsi_partial.analysis import predictions as P
from rsi_partial.experiments.detection import detection_experiment
from rsi_partial.experiments.engine import measure_ensemble
from rsi_partial.measure.diagnostics import rho2_measured
from rsi_partial.observables.cards_ad import AverageDisplacement
from rsi_partial.observables.chain_mag import ChainMagnetization
from rsi_partial.observables.ising_obs import AbsMagnetization
from rsi_partial.systems.cards import Cards
from rsi_partial.systems.product_chain import ProductChain
from rsi_partial.systems.ising2d import Ising2D, T_C

SEED = 3131


# ---------------------------------------------------------------------------
# 1. Off-equilibrium SNR laws rho^2(phi)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("phi", [0.5, 0.75, 1.0])
def test_chain_rho2_of_phi_exact(phi):
    rng = np.random.default_rng(SEED + int(100 * phi))
    res = measure_ensemble(ProductChain(), ChainMagnetization(), 4096, 6000, rng, [4096], phi=phi)
    measured = rho2_measured(res.O_full, res.Vx)
    assert abs(measured - P.Chain.rho2(phi)) < 0.03, (phi, measured, P.Chain.rho2(phi))


@pytest.mark.parametrize("phi", [0.5, 0.75, 1.0])
def test_cards_rho2_of_phi_leading_order(phi):
    """Cards rho^2(phi) ~ phi/(3/2-phi): leading order, Phase-3 matched to a few %."""
    rng = np.random.default_rng(SEED + int(100 * phi))
    res = measure_ensemble(Cards(), AverageDisplacement(), 2048, 6000, rng, [2048], phi=phi)
    measured = rho2_measured(res.O_full, res.Vx)
    predicted = P.Cards.rho2(phi)
    assert abs(measured - predicted) / predicted < 0.08, (phi, measured, predicted)


# ---------------------------------------------------------------------------
# 2. n-independence of the equilibrium SNR
# ---------------------------------------------------------------------------

def test_equilibrium_snr_size_independent():
    rng = np.random.default_rng(SEED + 7)
    for n in (512, 2048):
        res = measure_ensemble(Cards(), AverageDisplacement(), n, 6000, rng, [n], phi=1.0)
        assert abs(rho2_measured(res.O_full, res.Vx) - 2.0) < 0.06, n
    for n in (4096, 16384):
        res = measure_ensemble(ProductChain(), ChainMagnetization(), n, 6000, rng, [n], phi=1.0)
        assert abs(rho2_measured(res.O_full, res.Vx) - 1.0) < 0.04, n


# ---------------------------------------------------------------------------
# 3. Detection cost m* delta^2 -> kappa^2/rho^2, size-free (calibrated exact null)
# ---------------------------------------------------------------------------

def test_cards_detection_constant_and_size_free():
    """cards m* delta^2 ~ kappa^2/2 ~ 3.09; identical across n (size-free, a=0)."""
    target = P.Cards.detection_const()          # ~3.09
    md2 = {}
    for n in (2048, 4096):
        rng = np.random.default_rng(SEED + n)
        run = detection_experiment("cards", "AD", n, 0.10, P.Cards.AD_infinity(n), rng,
                                   R=3000, calibrate=True)
        md2[n] = run.mstar_delta2
        assert 2.4 < run.mstar_delta2 < 4.2, (n, run.mstar_delta2)
        assert abs(run.type_I_at_mstar - 0.05) < 0.02       # calibration holds
    # size-free: the two sizes agree within combined MC error (~15%)
    assert abs(md2[2048] - md2[4096]) / target < 0.20, md2


def test_chain_detection_constant():
    """chain m* delta^2 ~ kappa^2 ~ 6.18."""
    rng = np.random.default_rng(SEED + 11)
    run = detection_experiment("product_chain", "magnetization", 8192, 0.10, 0.5, rng,
                               R=3000, calibrate=True)
    assert 5.0 < run.mstar_delta2 < 7.8, run.mstar_delta2


# ---------------------------------------------------------------------------
# 4. Ising: high-T subsample variance ~ 1/m; near-T_c m-independent floor (M4: qualitative)
# ---------------------------------------------------------------------------

def _ising_var_vs_m(L, T, m_small, m_large, rng, R=1500):
    n = L * L
    res = measure_ensemble(Ising2D(cluster_flips_per_sample=6), AbsMagnetization(),
                           n, R, rng, [m_small, m_large], phi=1.0, T=T, batch=R)
    v_small = float(np.var(res.pools[m_small].Ohat, ddof=1))
    v_large = float(np.var(res.pools[m_large].Ohat, ddof=1))
    return v_small, v_large


def test_ising_high_T_variance_falls_like_one_over_m():
    """High-T: Var(M_hat_m) ~ 1/m, so Var(m=16)/Var(m=256) ~ 16 (short correlations)."""
    rng = np.random.default_rng(SEED + 21)
    v16, v256 = _ising_var_vs_m(24, 1.5 * T_C, 16, 256, rng)
    ratio = v16 / v256
    assert ratio > 4.0, ("high-T should fall ~1/m", ratio)


def test_ising_near_Tc_variance_is_floored():
    """Near T_c: an m-independent floor dominates, so Var(m=16)/Var(m=256) ~ 1 (not ~16)."""
    rng = np.random.default_rng(SEED + 22)
    v16, v256 = _ising_var_vs_m(24, T_C, 16, 256, rng)
    ratio = v16 / v256
    assert ratio < 2.5, ("near-T_c floor should NOT fall like 1/m", ratio)
