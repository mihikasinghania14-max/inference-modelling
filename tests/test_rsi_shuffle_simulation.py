"""Statistical and structural tests of the simulation code (seeded, so deterministic)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from rsi_shuffle import bayes, estimator, grid, shuffle, theory


def test_mixing_time_values():
    assert shuffle.mixing_time(512) == 1597                        # RSI Fig. 3 / Table
    assert [shuffle.mixing_time(n) for n in (2, 4, 8, 16, 32)] == [1, 3, 8, 22, 55]
    assert shuffle.mixing_time(512, "ceil") == 1598


def test_deck_is_a_permutation_and_ad_formula():
    rng = np.random.default_rng(0)
    for conv in ("independent", "distinct"):
        perm = shuffle.shuffle_deck(32, 100, rng, conv)
        assert sorted(perm.tolist()) == list(range(32))
        assert shuffle.average_displacement(perm) == pytest.approx(np.abs(perm - np.arange(32)).mean())
    perms = shuffle.shuffle_decks(16, 40, 50, rng)
    assert all(sorted(row.tolist()) == list(range(16)) for row in perms)


def test_pair_conventions_differ_at_n2():
    """With n = 2 and one step: independent pairs swap w.p. 1/2 (mean AD 0.5 = AD_inf);
    distinct pairs always swap (AD = 1).  This is the discrepancy documented as D1."""
    rng = np.random.default_rng(1)
    ind = shuffle.average_displacement(shuffle.shuffle_decks(2, 1, 20000, rng, "independent"))
    dis = shuffle.average_displacement(shuffle.shuffle_decks(2, 1, 2000, rng, "distinct"))
    assert abs(ind.mean() - 0.5) < 4 * ind.std() / math.sqrt(len(ind))
    assert np.all(dis == 1.0)


def test_batch_and_single_agree_in_distribution():
    n, t, R = 16, 22, 3000
    rng = np.random.default_rng(2)
    a = shuffle.average_displacement(shuffle.shuffle_decks(n, t, R, rng))
    b = np.array([shuffle.average_displacement(shuffle.shuffle_deck(n, t, rng)) for _ in range(R)])
    se = math.sqrt(a.var() / R + b.var() / R)
    assert abs(a.mean() - b.mean()) < 4.5 * se


def test_equilibrium_mean_at_mixing_time_large_n():
    """RSI Fig. 2: mean AD after t_mix matches AD_inf (here n = 256, 400 decks)."""
    n, R = 256, 400
    rng = np.random.default_rng(3)
    ad = shuffle.average_displacement(shuffle.shuffle_decks(n, 2 * shuffle.mixing_time(n), R, rng))
    se = ad.std(ddof=1) / math.sqrt(R)
    assert abs(ad.mean() - theory.ad_inf(n)) < 4.5 * se


def test_uniform_estimator_unbiased_and_first_m_biased():
    n, m, R = 64, 8, 4000
    rng = np.random.default_rng(4)
    perms = np.argsort(rng.random((R, n)), axis=1)               # uniform permutations
    d = shuffle.displacements(perms)
    u = estimator.uniform_sample_estimate(d, m, rng)
    f = estimator.first_m_window_estimate(d, m)
    ad_inf = theory.ad_inf(n)
    assert abs(u.mean() - ad_inf) < 4.5 * u.std(ddof=1) / math.sqrt(R)
    # RSI Corollary 2.7: boundary cards move further, so the first-m window is biased upward.
    assert f.mean() - ad_inf > 8 * f.std(ddof=1) / math.sqrt(R)


def test_estimator_variance_matches_exact_law():
    n, m, R = 64, 16, 6000
    rng = np.random.default_rng(5)
    perms = np.argsort(rng.random((R, n)), axis=1)
    est = estimator.uniform_sample_estimate(shuffle.displacements(perms), m, rng)
    assert est.var(ddof=1) == pytest.approx(theory.estimator_variance(n, m), rel=0.08)


def test_conjugate_posterior_algebra():
    post = bayes.conjugate_posterior(y=168.2, sigma_m=10.67, mu0=170.7, sigma0=100.0)
    # RSI Fig. 8 (n = 512, m = 128): posterior ~ N(168.23, 10.61^2)
    assert post.mean == pytest.approx(168.23, abs=0.02)
    assert post.sd == pytest.approx(10.61, abs=0.01)
    flat = bayes.conjugate_posterior(y=3.0, sigma_m=2.0, mu0=0.0, sigma0=math.inf)
    assert (flat.mean, flat.sd) == pytest.approx((3.0, 2.0))
    assert bayes.sampling_sd(512, 128) == pytest.approx(10.667, abs=1e-3)


def test_grid_small_run_columns_and_values(tmp_path):
    out = tmp_path / "grid.csv"
    rows = grid.run_grid(sizes=[64], taus=[1.0], R=400, seed=7, out=out,
                         m_values={64: [8, 64]})
    assert [r["M"] for r in rows] == [8, 64]
    assert out.read_text().splitlines()[0].split(",") == grid.COLUMNS
    for r in rows:
        assert r["mean_AD"] == pytest.approx(theory.ad_inf(64), rel=0.05)
        assert r["variance"] == pytest.approx(theory.estimator_variance(64, r["M"]), rel=0.3)
