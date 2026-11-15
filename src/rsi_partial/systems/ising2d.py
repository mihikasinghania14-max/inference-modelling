"""S3 -- 2-D Ising model on an ``L x L`` periodic lattice (``n = L^2``), ``J = 1``.

Two *distinct* algorithms, per the design brief's ground-truth/dynamics separation:

  * :meth:`sample_equilibrium` -- **Wolff single-cluster** sampler.  Bypasses critical
    slowing, so near-``T_c`` ground truth is never corrupted by finite burn-in (the fix
    for the Phase-3 caveat).  A Swendsen-Wang sweep is available as a cross-check.
  * :meth:`step` -- **single-site Glauber (heat-bath)** dynamics: the physical process
    whose convergence to equilibrium is what the partial-observation certificate tests.
    One ``step`` is one full sweep of ``n`` heat-bath site updates.

``T_c = 2 / ln(1 + sqrt(2)) ~ 2.269``.  Equilibrium observable values are analytic only
via Onsager (see :mod:`rsi_partial.analysis.predictions`); this class exposes no faked
metadata.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from ..core.interfaces import System
from ..core.types import GroundTruth

T_C = 2.0 / math.log(1.0 + math.sqrt(2.0))


class Ising2D(System):
    name = "ising2d"
    ground_truth = GroundTruth.MCMC_CLUSTER
    T_c = T_C                      # 2/ln(1+sqrt2) ~ 2.269 (exposed for convenience)

    def __init__(self, *, cluster_flips_per_sample: int = 4) -> None:
        # Number of Wolff cluster flips between recorded equilibrium samples (decorrelation).
        self.cluster_flips_per_sample = cluster_flips_per_sample
        self._mask_cache: dict[int, np.ndarray] = {}

    def _checkerboard(self, L: int) -> np.ndarray:
        """Cached boolean mask, True on one sublattice ((i+j) even).  Requires even ``L``:
        with periodic boundaries the square lattice is bipartite only for even ``L`` (an
        odd wrap creates same-colour neighbours and breaks the two-sublattice update)."""
        m = self._mask_cache.get(L)
        if m is None:
            if L % 2 != 0:
                raise ValueError(
                    f"checkerboard Glauber needs even L (periodic bipartite); got L={L}")
            ii, jj = np.indices((L, L))
            m = ((ii + jj) % 2 == 0)
            self._mask_cache[L] = m
        return m

    # ---- initial conditions ----------------------------------------------
    def initial_state(
        self, n: int, rng: np.random.Generator, *, kind: str = "ordered", **params: Any
    ) -> np.ndarray:
        L = self._L(n)
        if kind == "ordered":
            return np.ones((L, L), dtype=np.int8)
        if kind == "random":
            return np.where(rng.random((L, L)) < 0.5, 1, -1).astype(np.int8)
        raise ValueError(kind)

    @staticmethod
    def _L(n: int) -> int:
        L = int(round(math.sqrt(n)))
        if L * L != n:
            raise ValueError(f"n={n} is not a perfect square (Ising needs n=L^2)")
        return L

    # ---- Wolff cluster equilibrium sampler -------------------------------
    def _wolff_flip(self, s: np.ndarray, T: float, rng: np.random.Generator) -> None:
        """One Wolff single-cluster flip, in place."""
        L = s.shape[0]
        p_add = 1.0 - math.exp(-2.0 / T)  # J=1
        i0 = int(rng.integers(L)); j0 = int(rng.integers(L))
        seed_spin = s[i0, j0]
        stack = [(i0, j0)]
        s[i0, j0] = -seed_spin            # flip seed immediately; marks it visited
        while stack:
            i, j = stack.pop()
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ni = (i + di) % L; nj = (j + dj) % L
                if s[ni, nj] == seed_spin and rng.random() < p_add:
                    s[ni, nj] = -seed_spin
                    stack.append((ni, nj))

    #: Default Wolff warmup (cluster flips from a random start).  Measured to equilibrate
    #: the energy to <0.5% of Onsager for T >= T_c at L up to 32 (below T_c needs ~100
    #: flips to build the ordered domains; 20 was too short and left samples under-mixed).
    DEFAULT_WARMUP = 100

    def sample_equilibrium(
        self, n: int, rng: np.random.Generator, *, phi: float = 1.0, T: float = T_C,
        _warmup: int | None = None, **params: Any
    ) -> np.ndarray:
        """Draw one equilibrium configuration at temperature ``T`` via Wolff.

        A Wolff warmup from a random start is applied, then several cluster flips before
        the recorded sample.  Wolff decorrelates in ~O(1) cluster flips near ``T_c`` (no
        critical slowing), but from a random start below ``T_c`` it needs ~100 flips to
        build the ordered domains -- hence :data:`DEFAULT_WARMUP`.  For many samples at one
        ``(n,T)`` prefer :meth:`wolff_stream` (warm once, then decorrelate) -- far cheaper.

        This method draws **equilibrium only**.  ``phi<1`` is rejected rather than silently
        ignored (critique sec.1.1): Ising off-equilibrium states come from the physical
        dynamics via :meth:`glauber_relax_batch`, not from the equilibrium sampler.
        """
        if phi < 1.0:
            raise ValueError(
                "Ising sample_equilibrium draws equilibrium only; for off-equilibrium use "
                "glauber_relax_batch(..., sweeps=k) (checkerboard-Glauber relaxation)")
        L = self._L(n)
        s = np.where(rng.random((L, L)) < 0.5, 1, -1).astype(np.int8)
        warmup = _warmup if _warmup is not None else max(self.DEFAULT_WARMUP,
                                                         4 * self.cluster_flips_per_sample)
        for _ in range(warmup):
            self._wolff_flip(s, T, rng)
        for _ in range(self.cluster_flips_per_sample):
            self._wolff_flip(s, T, rng)
        return s

    def wolff_stream(
        self, n: int, rng: np.random.Generator, *, T: float = T_C, count: int,
        warmup: int = 100,
    ):
        """Yield ``count`` decorrelated equilibrium configs sharing one warmed chain.

        Cheaper than independent :meth:`sample_equilibrium` calls when many samples are
        needed at the same ``(n, T)``; still Wolff, so still rapid ground truth.
        """
        L = self._L(n)
        s = np.where(rng.random((L, L)) < 0.5, 1, -1).astype(np.int8)
        for _ in range(warmup):
            self._wolff_flip(s, T, rng)
        for _ in range(count):
            for _ in range(self.cluster_flips_per_sample):
                self._wolff_flip(s, T, rng)
            yield s.copy()

    # ---- single-site Glauber dynamics (the certified physical process) ---
    def step(self, state: np.ndarray, rng: np.random.Generator, *, T: float = T_C,
             **params: Any) -> np.ndarray:
        """One full **checkerboard heat-bath (Glauber) sweep**, fully vectorised.

        The square lattice with periodic boundaries is bipartite (even ``L``): every site's
        four neighbours lie on the opposite sublattice.  So all sites of one colour can be
        heat-bath-updated **simultaneously** with their (fixed) opposite-colour neighbours,
        then the other colour -- two NumPy operations instead of ``n`` Python iterations.
        This is exactly single-site Glauber (heat-bath, detailed balance preserved), just
        applied in a valid parallel order, and is what makes the ``L=64..128`` critical
        lattices reachable (critique sec.4.1).

        Accepts a single ``(L, L)`` lattice or a batch ``(b, L, L)`` (the mask broadcasts
        over the batch axis).  Returns a new array.
        """
        beta = 1.0 / T
        s = state.astype(np.int8, copy=True)
        L = s.shape[-1]
        if s.shape[-2] != L:
            raise ValueError("Ising lattice must be square")
        even = self._checkerboard(L)
        for color in (even, ~even):
            nb = (np.roll(s, 1, axis=-2) + np.roll(s, -1, axis=-2)
                  + np.roll(s, 1, axis=-1) + np.roll(s, -1, axis=-1))
            p_up = 1.0 / (1.0 + np.exp(-2.0 * beta * nb))
            flip = np.where(rng.random(s.shape) < p_up, 1, -1).astype(np.int8)
            s = np.where(color, flip, s).astype(np.int8)     # mask (L,L) broadcasts over batch
        return s

    def glauber_relax_batch(
        self, n: int, b: int, rng: np.random.Generator, *, T: float = T_C,
        sweeps: int, start: str = "random",
    ) -> np.ndarray:
        """Return ``(b, L, L)``: ``b`` lattices after ``sweeps`` checkerboard-Glauber sweeps
        from a common start type -- the **off-equilibrium** generator (the physical process
        whose convergence to equilibrium is being certified).

        ``start='random'`` is a hot (``T=inf``) start; ``start='ordered'`` a cold (all-up)
        start.  The resulting mixing level is a measured output (compare the observable to
        its equilibrium value), not an input -- for Ising ``phi`` is defined against the
        **measured** finite-``L`` equilibrium mean, and the well-posed detection observable
        is the **energy density** (nonzero ``O_inf = u(T)``), not ``|M|`` whose mean
        vanishes as ``L -> inf`` above ``T_c`` (critique sec.1.2).
        """
        L = self._L(n)
        if start == "random":
            s = np.where(rng.random((b, L, L)) < 0.5, 1, -1).astype(np.int8)
        elif start == "ordered":
            s = np.ones((b, L, L), dtype=np.int8)
        else:
            raise ValueError(f"start must be 'random' or 'ordered', got {start!r}")
        for _ in range(sweeps):
            s = self.step(s, rng, T=T)
        return s

    def t_mix(self, n: int, *, T: float = T_C, **params: Any) -> int | None:
        return None  # not analytic; measured via tau_int

    def O_infinity(self, n: int, **params: Any) -> float | None:
        return None  # measured (Onsager references live in analysis.predictions)
