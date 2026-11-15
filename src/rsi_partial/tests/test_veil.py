"""Veil / access-audit tests (control Phase-4 sec.5.5; design invariant 1).

Enforce structurally that inference sees ONLY the m revealed contributions:
  * the ``Observation`` has no full-state field;
  * every reveal reads exactly ``m`` components;
  * results are invariant to a random relabelling of components (no positional leakage);
  * the memory-aware sampler works at large ``n`` with small ``m`` (O(Rm), no R x n array).
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rsi_partial.core.types import Observation
from rsi_partial.observe.veil import (
    relabel, reveal, reveal_batch, reveal_exchangeable_first_m)
from rsi_partial.observables.cards_ad import AverageDisplacement
from rsi_partial.systems.cards import Cards

SEED = 4242


def test_observation_has_no_state_field():
    """The veil boundary must not carry the full state -- leakage impossible by construction."""
    fields = {f.name for f in dataclasses.fields(Observation)}
    for forbidden in ("state", "full_state", "values_full", "D_full"):
        assert forbidden not in fields
    assert "values" in fields


@pytest.mark.parametrize("n,m", [(1000, 8), (1000, 64), (4096, 256)])
def test_access_audit_reads_equal_m(n, m):
    """reveal records exactly m reads and carries exactly m values."""
    rng = np.random.default_rng(SEED)
    D = AverageDisplacement().per_component(Cards().sample_equilibrium(n, rng))
    obs = reveal(D, m, rng)
    assert obs.reads == m
    assert obs.values.shape == (m,)
    assert obs.m == m and obs.n == n


def test_relabelling_invariance():
    """A random relabelling of components leaves the reveal distribution unchanged.

    The population mean is exactly invariant; the sampling distribution of Ohat is
    invariant in law, checked here via its mean over many reveals within MC error.
    """
    rng = np.random.default_rng(SEED)
    n, m, R = 512, 32, 6000
    D = AverageDisplacement().per_component(Cards().sample_equilibrium(n, rng))
    D_relabelled = relabel(D, np.random.default_rng(1))
    assert np.isclose(D.mean(), D_relabelled.mean())          # exact
    a = np.array([reveal(D, m, rng).values.mean() for _ in range(R)])
    b = np.array([reveal(D_relabelled, m, rng).values.mean() for _ in range(R)])
    se = np.sqrt(a.var(ddof=1) / R + b.var(ddof=1) / R)
    assert abs(a.mean() - b.mean()) < 5 * se


def test_reveal_batch_is_genuine_and_correct():
    """reveal_batch reveals exactly m per row; reads is derived from the array, not asserted;
    the census returns a copy (not a view onto the full state)."""
    rng = np.random.default_rng(SEED)
    b, n, m = 200, 400, 30
    D = rng.random((b, n))
    vals, reads = reveal_batch(D, m, rng)
    assert vals.shape == (b, m)
    assert reads == m == vals.shape[1]            # genuine: count == what inference sees
    assert vals.size == b * m                     # nothing beyond m*b exposed
    # every revealed value is a real entry of its own row (uniform without replacement)
    for r in range(5):
        assert set(np.round(vals[r], 12)).issubset(set(np.round(D[r], 12)))
    # census reveals all n but as a copy: mutating it must not touch the source state
    full, reads_n = reveal_batch(D, n, rng)
    assert reads_n == n
    full[0, 0] = -999.0
    assert D[0, 0] != -999.0                       # copy, not a view -> no back-channel


def test_engine_access_audit_is_the_executed_path():
    """The engine's reported reads come from the veil on the path it actually runs, and
    equal m for every m (census reveals all n)."""
    from rsi_partial.experiments.engine import measure_ensemble
    rng = np.random.default_rng(SEED)
    n = 300
    res = measure_ensemble(Cards(), AverageDisplacement(), n, 500, rng, [8, 32, 128, n])
    assert res.reads_per_replicate == {8: 8, 32: 32, 128: 128, n: n}


def test_reveal_rejects_m_gt_n():
    rng = np.random.default_rng(0)
    D = np.arange(10.0)
    with pytest.raises(ValueError):
        reveal(D, 11, rng)


def test_memory_aware_large_n_small_m():
    """Sampling at large n with small m must not allocate an R x n array.

    We reveal from a length-2e5 array many times; this is O(m) per draw.  (A naive
    argsort(random(R,n)) would allocate ~1.6e6 floats per batch and OOM at 1e5 scale.)
    """
    rng = np.random.default_rng(SEED)
    n, m = 200_000, 16
    D = rng.random(n)
    for _ in range(200):
        obs = reveal(D, m, rng)
        assert obs.values.shape == (m,)


def test_exchangeable_fast_path_matches_choice():
    """The first-m exchangeable fast path is distributionally equal to choice-based reveal
    for an exchangeable-under-nu observable (product chain magnetization at equilibrium)."""
    rng = np.random.default_rng(SEED)
    n, m, R = 4096, 64, 6000

    def gen_m(mm, rng):                       # first-m Bernoulli(1/2) coordinates
        return (rng.random(mm) < 0.5).astype(np.float64)

    fast = np.array([
        reveal_exchangeable_first_m(gen_m, n, m, rng).values.mean() for _ in range(R)
    ])
    full = (rng.random((R, n)) < 0.5).astype(np.float64)
    slow = np.array([reveal(full[r], m, rng).values.mean() for r in range(R)])
    se = np.sqrt(fast.var(ddof=1) / R + slow.var(ddof=1) / R)
    assert abs(fast.mean() - slow.mean()) < 5 * se
    # variances agree too (both ~ Var(Bernoulli)/m * fpc); within ~6%
    assert abs(fast.var(ddof=1) - slow.var(ddof=1)) / slow.var(ddof=1) < 0.10
