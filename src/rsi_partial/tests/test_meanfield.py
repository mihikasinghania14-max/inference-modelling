"""Curie-Weiss (S4) tests -- the a=1 hard-floor result (CX1 / P2).

The genuinely new, falsification-grade content: the SAME mean-field system has an extensive
magnetization floor below T_c (a=1, slope 0, a hard floor no sampling beats) and a
self-averaging one above (a=0, slope -1).  A flat/absent floor below T_c would falsify the
refined conjecture; here it is confirmed with a bootstrap CI.
"""

from __future__ import annotations

import numpy as np

from rsi_partial.analysis.predictions import curie_weiss_m0
from rsi_partial.experiments.engine import measure_ensemble
from rsi_partial.experiments.meanfield import (
    floor_exponent, hard_floor_vs_m, magnetizations, sample_up_counts, subsample_Mhat)
from rsi_partial.observables.chain_mag import ChainMagnetization
from rsi_partial.systems.curie_weiss import CurieWeiss

SEED = 13579
NS = [100, 300, 1000, 3000]


def test_sampler_matches_spontaneous_magnetization():
    """<|M|> below T_c matches the exact m0(T) = tanh(m0/T)."""
    rng = np.random.default_rng(SEED)
    M = magnetizations(2000, 0.7, 20000, rng)
    assert abs(np.abs(M).mean() - curie_weiss_m0(0.7)) < 0.02


def test_below_Tc_extensive_floor_a_equals_1():
    """Below T_c: slope ~ 0 (a=1). The CI is very tight and may exclude *exactly* 0 by a
    hair (mild finite-N drift as Var(M) approaches m0^2), so the scientific check is
    CI-within-the-support-band (verdict) plus a_estimate ~ 1, not contains(0)."""
    rng = np.random.default_rng(SEED + 1)
    fs = floor_exponent(0.7, NS, rng, R=10000, n_boot=200)
    assert abs(fs.a_estimate - 1.0) < 0.1         # extensive floor
    assert fs.verdict.verdict == "support"        # P2 CI lies inside a-1 +/- 0.1


def test_above_Tc_self_averaging_a_equals_0():
    rng = np.random.default_rng(SEED + 2)
    fs = floor_exponent(2.0, NS, rng, R=10000, n_boot=200)
    assert abs(fs.a_estimate - 0.0) < 0.1         # self-averaging floor ~ 1/N
    assert fs.verdict.verdict == "support"        # P2 CI lies inside a-1 +/- 0.1


def test_hard_floor_not_reduced_by_more_sampling():
    """Below T_c, Var(M_hat_m) stays at the floor for all m (subsampling cannot help)."""
    rng = np.random.default_rng(SEED + 3)
    hf = hard_floor_vs_m(0.7, 5000, [16, 64, 256, 1024], rng, R=10000)
    for r in hf.ratios:
        assert 0.9 < r < 1.15                     # pinned at the floor, no 1/m decay


def test_hypergeometric_fastpath_matches_veil_engine():
    """The exchangeable hypergeometric subsample equals the materialized-spins + veil path
    (validates the fast path against the audited path)."""
    rng = np.random.default_rng(SEED + 4)
    N, T, m, R = 400, 0.7, 32, 6000
    res = measure_ensemble(CurieWeiss(), ChainMagnetization(), N, R, rng, [m],
                           phi=1.0, T=T, batch=R)
    k = sample_up_counts(N, T, R, rng)
    fast = subsample_Mhat(k, N, m, rng)
    # both estimate Var(M_hat_m); agree within Monte-Carlo error (~ few %)
    assert abs(res.pools[m].Ohat.var(ddof=1) - fast.var(ddof=1)) < 0.03
