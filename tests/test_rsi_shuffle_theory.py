"""Exact checks of the closed forms in ``rsi_shuffle.theory`` by enumerating S_n."""

from __future__ import annotations

import itertools
import math
from fractions import Fraction

import numpy as np
import pytest

from rsi_shuffle import theory


def _all_displacements(n):
    return [tuple(abs(p[i] - i) for i in range(n)) for p in itertools.permutations(range(n))]


@pytest.mark.parametrize("n", [2, 3, 4, 5, 6, 7])
def test_equilibrium_moments_by_enumeration(n):
    D = _all_displacements(n)
    N = len(D)
    F = [sum(d) for d in D]
    mean_F = Fraction(sum(F), N)
    var_F = Fraction(sum(f * f for f in F), N) - mean_F**2
    assert mean_F / n == theory.ad_inf(n, exact=True)                       # RSI Prop. 2.6
    assert var_F == theory.var_footrule(n, exact=True)                      # (n+1)(2n^2+7)/45
    # E_nu[V_x]: mean within-deck population variance of the displacements.
    EV = sum(Fraction(sum(x * x for x in d), n) - Fraction(sum(d), n) ** 2 for d in D) / N
    assert EV == theory.expected_within_deck_variance(n, exact=True)
    # single-card second moment (RSI Section 5.4)
    second = Fraction(sum(sum(x * x for x in d) for d in D), N * n)
    assert second == theory.second_moment_displacement(n, exact=True)


@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_exact_estimator_variance_lemma1(n):
    """Var(AD_hat_m) = (n-m)/(m(n-1)) E[V_x] + Var(AD), for every m (extension Lemma 1)."""
    D = _all_displacements(n)
    for m in range(1, n + 1):
        vals = []
        for d in D:
            for S in itertools.combinations(range(n), m):
                vals.append(Fraction(sum(d[i] for i in S), m))
        mean = sum(vals) / len(vals)
        var = sum(v * v for v in vals) / len(vals) - mean**2
        assert mean == theory.ad_inf(n, exact=True)                          # Lemma 2.10
        assert var == theory.estimator_variance(n, m, exact=True)


def test_relative_floor_limit_and_closed_form():
    for n in (3, 10, 100):
        want = theory.var_ad(n, exact=True) / theory.ad_inf(n, exact=True) ** 2
        assert theory.relative_floor(n, exact=True) == want
    assert abs(1000 * theory.relative_floor(1000) - 0.4) < 1e-3          # n R_n -> 2/5


def test_leading_order_variance_law():
    """RSI Theorem 5.1: Var(AD_hat_m) ~ n^2/(18 m) for m << n."""
    n, m = 4096, 64
    ratio = theory.estimator_variance(n, m) / (n * n / (18 * m))
    assert abs(ratio - 1) < 0.02


def test_sample_size_law_values():
    # RSI Section 5.4: about 48, 192, 768 for 20%, 10%, 5%.
    assert [round(theory.m_star_rsi(e)) for e in (0.20, 0.10, 0.05)] == [48, 192, 768]
    assert abs(theory.kappa() - 2.486) < 1e-3
    # certification impossible below kappa*sqrt(R_n)
    n = 1000
    floor_delta = theory.kappa() * math.sqrt(theory.relative_floor(n))
    assert math.isinf(theory.certification_cost(0.9 * floor_delta, n))
    assert math.isfinite(theory.certification_cost(1.5 * floor_delta, n))


def _exact_mean_ad_path(n, tmax):
    """Exact E[AD(pi_t)] for the independent-pair chain from the identity, n small."""
    states = list(itertools.permutations(range(n)))
    index = {s: k for k, s in enumerate(states)}
    # one-step kernel on `deck` (deck[pos] = card): swap positions a, b, each pair w.p. 1/n^2
    P = [[Fraction(0)] * len(states) for _ in states]
    for s in states:
        for a in range(n):
            for b in range(n):
                t = list(s)
                t[a], t[b] = t[b], t[a]
                P[index[s]][index[tuple(t)]] += Fraction(1, n * n)
    ad = []
    for s in states:  # AD of a deck: |position - card|
        ad.append(Fraction(sum(abs(p - c) for p, c in enumerate(s)), n))
    dist = [Fraction(0)] * len(states)
    dist[index[tuple(range(n))]] = Fraction(1)
    out = [sum(d * a for d, a in zip(dist, ad))]
    for _ in range(tmax):
        dist = [sum(dist[i] * P[i][j] for i in range(len(states))) for j in range(len(states))]
        out.append(sum(d * a for d, a in zip(dist, ad)))
    return out


def test_mean_law_exact_n4():
    n, tmax = 4, 6
    exact = _exact_mean_ad_path(n, tmax)
    ad_inf = theory.ad_inf(n, exact=True)
    for t, v in enumerate(exact):
        assert float(v / ad_inf) == pytest.approx(theory.mean_law_phi(t, n), abs=1e-12)
