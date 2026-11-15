"""Ising ground-truth tests: Wolff cluster sampler vs Onsager; Glauber stationarity.

The Wolff sampler is validated where it is unambiguous:
  * energy density below and above ``T_c`` (finite-size corrections small there);
  * spontaneous magnetization below ``T_c`` (Onsager's exact ``m(T)``).
At ``T_c`` itself the infinite-lattice Onsager values are singular, so finite-``L`` samples
are expected to differ -- that regime is validated only qualitatively (energy near the
critical value), consistent with modelling note M4.

Glauber (the certified physical dynamics) is checked for detailed balance via
stationarity: started from a Wolff-equilibrium ensemble, Glauber sweeps must leave the
energy unchanged within Monte-Carlo error.
"""

from __future__ import annotations

import numpy as np
import pytest

from rsi_partial.analysis.predictions import (
    ising_onsager_energy, ising_onsager_magnetization)
from rsi_partial.observables.ising_obs import AbsMagnetization, EnergyDensity
from rsi_partial.systems.ising2d import Ising2D, T_C

SEED = 77


def _wolff_means(L, T, count, flips, warmup, obs, rng):
    ise = Ising2D(cluster_flips_per_sample=flips)
    n = L * L
    return np.array([obs.population_value(s)
                     for s in ise.wolff_stream(n, rng, T=T, count=count, warmup=warmup)])


@pytest.mark.parametrize("L,T,tol", [
    (24, 2.0, 0.02),            # below T_c: ordered phase, Wolff efficient
    (24, 1.15 * T_C, 0.02),     # above T_c but moderate: still efficient
])
def test_wolff_energy_matches_onsager(L, T, tol):
    rng = np.random.default_rng(SEED + int(1000 * T))
    E = _wolff_means(L, T, count=800, flips=6, warmup=150, obs=EnergyDensity(), rng=rng)
    onsager = ising_onsager_energy(T)
    assert abs(E.mean() - onsager) < abs(onsager) * tol, (E.mean(), onsager)


def test_wolff_magnetization_below_Tc():
    """<|M|> at T=2.0 matches Onsager spontaneous magnetization within finite-size error."""
    rng = np.random.default_rng(SEED + 5)
    M = _wolff_means(32, 2.0, count=800, flips=5, warmup=150, obs=AbsMagnetization(), rng=rng)
    onsager = ising_onsager_magnetization(2.0)        # ~0.9113
    assert abs(M.mean() - onsager) < 0.02, (M.mean(), onsager)


def test_wolff_high_T_energy():
    """High-T energy (Wolff needs more decorrelation there; clusters shrink)."""
    rng = np.random.default_rng(SEED + 9)
    T = 1.5 * T_C
    E = _wolff_means(24, T, count=1000, flips=15, warmup=200, obs=EnergyDensity(), rng=rng)
    onsager = ising_onsager_energy(T)
    assert abs(E.mean() - onsager) < abs(onsager) * 0.03, (E.mean(), onsager)


def test_glauber_preserves_equilibrium():
    """Detailed-balance sanity: checkerboard-Glauber sweeps from Wolff equilibrium leave
    <E> stationary (warm one chain, cheap)."""
    rng = np.random.default_rng(SEED + 3)
    ise = Ising2D(cluster_flips_per_sample=6)
    en = EnergyDensity()
    T, L, n = 2.0, 24, 24 * 24
    eq = [c.copy() for c in ise.wolff_stream(n, rng, T=T, count=300, warmup=150)]
    E0 = np.array([en.population_value(c) for c in eq])
    E1 = np.empty(len(eq))
    for r, c in enumerate(eq):
        s = c.copy()
        for _ in range(15):
            s = ise.step(s, rng, T=T)
        E1[r] = en.population_value(s)
    se = np.sqrt(E0.var(ddof=1) / len(eq) + E1.var(ddof=1) / len(eq))
    assert abs(E0.mean() - E1.mean()) < 4 * se, (E0.mean(), E1.mean(), se)


def test_checkerboard_step_batched_matches_single():
    """The vectorised step accepts a batch (b,L,L) and each lattice is a valid update;
    a batch of identical lattices with a shared rng stream stays self-consistent in shape."""
    rng = np.random.default_rng(SEED + 8)
    ise = Ising2D()
    L, n = 16, 256
    s1 = ise.initial_state(n, rng, kind="random")
    batch = np.stack([s1, s1, s1])
    out = ise.step(batch, rng, T=2.5)
    assert out.shape == (3, L, L)
    assert out.dtype == np.int8 and set(np.unique(out)).issubset({-1, 1})


def test_odd_L_checkerboard_rejected():
    """Checkerboard Glauber requires even L (periodic bipartite); odd L must raise."""
    rng = np.random.default_rng(0)
    ise = Ising2D()
    s = np.ones((15, 15), dtype=np.int8)
    with pytest.raises(ValueError):
        ise.step(s, rng, T=2.5)


def test_offequilibrium_relaxes_toward_equilibrium():
    """Hot-start Glauber relaxation: energy magnitude grows monotonically toward the
    equilibrium value as sweeps increase (a genuine off-equilibrium trajectory, phi<1)."""
    rng = np.random.default_rng(SEED + 12)
    ise = Ising2D(cluster_flips_per_sample=6)
    en = EnergyDensity()
    T, L, n = 2.0, 24, 24 * 24
    eq = np.mean([en.population_value(c) for c in ise.wolff_stream(n, rng, T=T, count=200, warmup=150)])
    means = []
    for sw in (2, 10, 40):
        st = ise.glauber_relax_batch(n, 200, rng, T=T, sweeps=sw, start="random")
        means.append(np.mean([en.population_value(st[r]) for r in range(st.shape[0])]))
    # monotone toward equilibrium, and still short of it at few sweeps (phi<1)
    assert means[0] > means[1] > means[2], means      # energy decreasing (more negative)
    assert means[2] > eq                              # 40 sweeps not yet fully equilibrated
    assert means[2] / eq > 0.85                       # but well on the way (phi_E > 0.85)


def test_ising_sample_equilibrium_rejects_offeq():
    """phi<1 to the equilibrium sampler raises instead of silently returning equilibrium."""
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError):
        Ising2D().sample_equilibrium(256, rng, phi=0.7, T=2.0)


def test_wolff_seed_reproducible():
    """Same seed -> identical Wolff sample (reproducibility)."""
    n = 16 * 16
    a = Ising2D().sample_equilibrium(n, np.random.default_rng(123), T=2.3)
    b = Ising2D().sample_equilibrium(n, np.random.default_rng(123), T=2.3)
    assert np.array_equal(a, b)
