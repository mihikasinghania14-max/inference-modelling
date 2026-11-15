"""Self-calibration tests (critique sec.2.4): recover BOTH rho^2 and the floor R_n from
repeated subsamples, with no system knowledge -- the deployability claim.

Recovering R_n is what lets the deployable path know whether certification is even feasible
(delta <= kappa*sqrt(R_n)); the single-subsample path could not do this.
"""

from __future__ import annotations

import numpy as np

from rsi_partial.analysis.predictions import Cards
from rsi_partial.experiments.engine import measure_ensemble
from rsi_partial.measure.diagnostics import self_calibrate_multi
from rsi_partial.observables.cards_ad import AverageDisplacement
from rsi_partial.observables.chain_mag import ChainMagnetization
from rsi_partial.systems.cards import Cards as CardsSys
from rsi_partial.systems.curie_weiss import CurieWeiss

SEED = 555


def test_recovers_rho2_and_predicts_mstar_cards():
    """Cards a=0: rho2_hat ~ 2, floor ~ 0, and m*_hat ~ the true detection m* (~309 at 0.1)."""
    rng = np.random.default_rng(SEED)
    n = 2048
    res = measure_ensemble(CardsSys(), AverageDisplacement(), n, 6000, rng,
                           [16, 32, 64, 128, 256], phi=1.0)
    sc = self_calibrate_multi({m: p.Ohat for m, p in res.pools.items()},
                              Cards.AD_infinity(n), delta=0.10)
    assert abs(sc.rho2_hat - 2.0) < 0.2
    assert sc.R_n_hat < 1e-2                     # a=0: floor negligible at these m
    assert 250 < sc.mstar_hat < 380             # ~ kappa^2/(2*0.01) = 309


def test_recovers_hard_floor_and_flags_infeasible_curie_weiss():
    """Curie-Weiss below T_c (a=1): R_n_hat ~ O(1) (the hard floor) and m*_hat = inf."""
    rng = np.random.default_rng(SEED + 1)
    N = 4096
    res = measure_ensemble(CurieWeiss(), ChainMagnetization(), N, 4000, rng,
                           [16, 64, 256, 1024], phi=1.0, T=0.7, batch=4000)
    scale = float(np.abs(res.O_full).mean())    # |M| magnitude as the reference (E[M]=0)
    sc = self_calibrate_multi({m: p.Ohat for m, p in res.pools.items()}, scale, delta=0.30)
    assert 0.6 < sc.R_n_hat < 1.5               # extensive floor ~ 1 (in |M|^2 units)
    assert sc.mstar_hat == float("inf")         # certification infeasible -> correctly flagged
