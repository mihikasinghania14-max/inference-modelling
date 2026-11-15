"""Equilibration-gate tests (critique sec.1.3-1.4; Phase-4 negative control sec.5.2).

The certificate must DETECT under-burn-in, not be fooled by it: a two-initial-condition
check flags a chain that has not relaxed, and the generation gate refuses to record it.
"""

from __future__ import annotations

import numpy as np
import pytest

from rsi_partial.experiments.equilibration import (
    gated_equilibrium, two_start_check)
from rsi_partial.observables.ising_obs import EnergyDensity
from rsi_partial.systems.ising2d import Ising2D, T_C

SEED = 8642
L = 24
N = L * L
T_HIGH = 1.6 * T_C          # fast mixing -> adequate burn reaches stationarity quickly


def test_gate_flags_under_burn():
    rng = np.random.default_rng(SEED)
    chk = two_start_check(Ising2D(), EnergyDensity(), N, 250, rng, T=T_HIGH, sweeps=2)
    assert not chk.passed                       # ordered vs random disagree -> flagged


def test_gate_passes_adequate_burn():
    rng = np.random.default_rng(SEED + 1)
    chk = two_start_check(Ising2D(), EnergyDensity(), N, 250, rng, T=T_HIGH, sweeps=40)
    assert chk.passed                           # both starts relaxed to the same energy


def test_gated_equilibrium_raises_on_under_burn():
    rng = np.random.default_rng(SEED + 2)
    with pytest.raises(RuntimeError):
        gated_equilibrium(Ising2D(), EnergyDensity(), N, 200, rng,
                          T=T_HIGH, burn_sweeps=2, check_R=250)


def test_gated_equilibrium_returns_on_adequate_burn():
    rng = np.random.default_rng(SEED + 3)
    pool = gated_equilibrium(Ising2D(), EnergyDensity(), N, 200, rng,
                             T=T_HIGH, burn_sweeps=40, check_R=250)
    assert pool.shape == (200,)
    assert -1.0 < pool.mean() < -0.3            # a sane high-T energy density
