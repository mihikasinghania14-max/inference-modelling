"""S2 -- Random-transposition shuffle (cards); the RSI theory anchor.

State is a 0-indexed permutation ``pi`` of ``{0,...,n-1}`` (``pi[i]`` = image of
position ``i``).  Displacement per card ``D_i = |pi[i] - i|``; the average-displacement
observable ``AD = (1/n) sum_i D_i`` lives in :mod:`rsi_partial.observables.cards_ad`.

Exact results (RSI paper; all **[EST]**):
  * ``AD_inf = (n^2 - 1) / (3 n)``                                   (Prop. 2.6)
  * mean law ``E[AD(t)] = (1 - (1 - 2/n)^t) AD_inf``, i.e. ``phi(t) = 1 - (1-2/n)^t``
  * equilibrium per-card SNR ``rho^2(1) = 2``.

Dynamics :meth:`step` draws two positions ``a, b`` independently and uniformly (with
replacement, so ``a == b`` self-loops occur with prob 1/n, making the chain aperiodic --
RSI Def. 2.1 / App. C.2) and swaps them.  Off-equilibrium states at a target mixing
level are produced by running this physical dynamics from the identity for the number
of steps implied by the exact mean law (:meth:`relax_to_phi`).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from ..core.interfaces import System
from ..core.types import GroundTruth


class Cards(System):
    name = "cards"
    ground_truth = GroundTruth.EXACT_MARGINAL

    def sample_equilibrium(
        self, n: int, rng: np.random.Generator, *, phi: float = 1.0, **params: Any
    ) -> np.ndarray:
        if phi >= 1.0:
            return rng.permutation(n).astype(np.int64)
        # Off-equilibrium: exact marginal is not a product form, so use the physical
        # dynamics from the identity.  (The mean law is exact; the joint law is not
        # analytically constructible, hence we run the actual certified process.)
        return self.relax_to_phi(n, rng, phi=phi, **params)

    def initial_state(
        self, n: int, rng: np.random.Generator, *, kind: str = "ordered", **params: Any
    ) -> np.ndarray:
        if kind == "ordered":
            return np.arange(n, dtype=np.int64)     # identity permutation (phi = 0)
        if kind == "random":
            return rng.permutation(n).astype(np.int64)
        raise ValueError(kind)

    def step(self, state: np.ndarray, rng: np.random.Generator, **params: Any) -> np.ndarray:
        n = state.shape[0]
        a = int(rng.integers(n))
        b = int(rng.integers(n))
        state[a], state[b] = state[b], state[a]
        return state

    def steps_for_phi(self, n: int, phi: float) -> int:
        """Number of transpositions ``t`` s.t. ``1 - (1-2/n)^t`` is closest to ``phi``."""
        if phi <= 0.0:
            return 0
        if phi >= 1.0:
            # 'equilibrium via dynamics' is not used; caller should draw uniformly.
            phi = 1.0 - 1e-12
        base = 1.0 - 2.0 / n
        # solve (1-2/n)^t = 1 - phi  ->  t = log(1-phi)/log(1-2/n)
        t = math.log(1.0 - phi) / math.log(base)
        return max(0, int(round(t)))

    def relax_to_phi(
        self, n: int, rng: np.random.Generator, *, phi: float, **params: Any
    ) -> np.ndarray:
        t = self.steps_for_phi(n, phi)
        state = self.initial_state(n, rng, kind="ordered")
        # Vectorised block of t transpositions (with-replacement pairs, incl. self-loops).
        if t > 0:
            a = rng.integers(n, size=t)
            b = rng.integers(n, size=t)
            for k in range(t):
                ak, bk = int(a[k]), int(b[k])
                state[ak], state[bk] = state[bk], state[ak]
        return state

    def true_phi(self, t: int, n: int, **params: Any) -> float | None:
        return 1.0 - (1.0 - 2.0 / n) ** t

    def t_mix(self, n: int, **params: Any) -> int | None:
        # RSI Eq. 2.1: t_mix = round( (1/2) n ln n ).
        return int(round(0.5 * n * math.log(n)))

    def O_infinity(self, n: int, **params: Any) -> float | None:
        # AD observable equilibrium mean.
        return (n * n - 1) / (3.0 * n)

    def rho2(self, phi: float, **params: Any) -> float | None:
        # Equilibrium value is EXACT (=2).  The off-equilibrium form phi/(3/2 - phi)
        # is empirical/leading-order and therefore lives in analysis.predictions, not
        # here.  Return the exact equilibrium value only; None off equilibrium.
        if phi >= 1.0:
            return 2.0
        return None

    def susceptibility_exponent_a(self, **params: Any) -> float | None:
        # a = 0 via R_n ~ 1/n (empirical/strongly-suggested).  Returned as metadata
        # for scoring convenience; the estimator never uses it.
        return 0.0
