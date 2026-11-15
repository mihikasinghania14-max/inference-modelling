"""Variance-identity tests (the spine, identity star).

Confirm EMPIRICALLY that, across (n, m) and systems,

    Var(Ohat_m) == (1/m) (n-m)/(n-1) E_x[V_x]  +  Var_x(O),

with both terms measured INDEPENDENTLY (the sampling term from within-config variances,
the floor from full-system population values).  Also checks the census limit m=n where the
sampling term vanishes and Var(Ohat_n) == Var_x(O) exactly.
"""

from __future__ import annotations

import numpy as np
import pytest

from rsi_partial.experiments.engine import measure_ensemble
from rsi_partial.inference.estimator import variance_split
from rsi_partial.observables.cards_ad import AverageDisplacement
from rsi_partial.observables.chain_mag import ChainMagnetization
from rsi_partial.systems.cards import Cards
from rsi_partial.systems.product_chain import ProductChain

SEED = 909090


@pytest.mark.parametrize("system,observable,n,phi", [
    (Cards(), AverageDisplacement(), 512, 1.0),
    (ProductChain(), ChainMagnetization(), 4096, 0.7),
])
def test_variance_identity(system, observable, n, phi):
    rng = np.random.default_rng(SEED)
    R = 12000
    m_grid = [16, 32, 64, 128, 256]
    res = measure_ensemble(system, observable, n, R, rng, m_grid, phi=phi)
    for m in m_grid:
        vs = variance_split(res.pools[m].Ohat, res.Vx, res.O_full, m, n)
        ratio = vs.var_Ohat / vs.predicted_sum
        assert 0.90 < ratio < 1.10, (system.name, m, ratio, vs)


def test_census_variance_equals_floor():
    """At m=n the sampling term is 0, so Var(Ohat_n) == Var_x(O) exactly (same draws)."""
    rng = np.random.default_rng(SEED + 1)
    n, R = 256, 4000
    res = measure_ensemble(Cards(), AverageDisplacement(), n, R, rng, [n], phi=1.0)
    vs = variance_split(res.pools[n].Ohat, res.Vx, res.O_full, n, n)
    assert vs.sampling_term == 0.0
    # Ohat_n == O_full identically, so their variances coincide to machine precision.
    assert np.isclose(vs.var_Ohat, vs.floor, rtol=1e-9)


def test_floor_is_m_independent():
    """The measured floor (Var_x(O)) does not depend on m -- it is a state-level quantity."""
    rng = np.random.default_rng(SEED + 2)
    n, R = 512, 8000
    res = measure_ensemble(Cards(), AverageDisplacement(), n, R, rng, [16, 64, 256], phi=1.0)
    floor = float(np.var(res.O_full, ddof=1))
    for m in (16, 64, 256):
        vs = variance_split(res.pools[m].Ohat, res.Vx, res.O_full, m, n)
        assert vs.floor == pytest.approx(floor)
