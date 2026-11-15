"""S4 -- Curie-Weiss (mean-field) Ising on the complete graph, ``N`` spins.

Hamiltonian ``H = -(J/2N) (sum_i s_i)^2`` (each spin couples to the mean field), so the
energy per spin is ``-J m^2 / 2`` with ``m = (1/N) sum_i s_i``.  With ``J=1`` the critical
temperature is ``T_c = J = 1``.  This is the CX1 counterexample: below ``T_c`` the
magnetization susceptibility is **extensive**, ``chi_N[M] = N Var(M) = Theta(N)`` (``a=1``),
because the equilibrium magnetization distribution is bimodal at ``+/- m0(T)`` and
``Var_nu(M) -> m0^2 = Theta(1)`` does not self-average.  The prediction is a **hard,
non-vanishing precision floor** ``delta_floor = Theta(1)`` -- no sample of any size can
certify the magnetization below it.

Exact rapid equilibrium sampler (no MCMC, no critical slowing): with all-to-all coupling,
the energy depends only on the total magnetization, so given ``M`` the spins are uniform
over all configurations with that sum.  Hence we (1) sample the number of up-spins
``k = N_+`` from its exact 1-D law ``P(k) proportional to C(N,k) exp(beta J N m(k)^2 / 2)``,
then (2) place ``k`` up-spins at uniformly random positions.  This draws exact equilibrium
configurations (both symmetry sectors) in ``O(N)``.

The physical dynamics :meth:`step` is mean-field heat-bath (Glauber): each spin relaxes in
the instantaneous mean field ``h = J m``.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.special import gammaln

from ..core.interfaces import System
from ..core.types import GroundTruth

T_C_CW = 1.0                    # J = 1


class CurieWeiss(System):
    name = "curie_weiss"
    ground_truth = GroundTruth.EXACT_MARGINAL      # exact 1-D magnetization law
    T_c = T_C_CW

    def __init__(self, *, J: float = 1.0) -> None:
        self.J = J

    # ---- exact equilibrium magnetization law -----------------------------
    def _logP_k(self, N: int, T: float) -> np.ndarray:
        """Unnormalised log-probability of ``k`` up-spins, ``k = 0..N`` (exact)."""
        k = np.arange(N + 1)
        m = (2.0 * k - N) / N
        logbinom = gammaln(N + 1) - gammaln(k + 1) - gammaln(N - k + 1)
        beta = 1.0 / T
        return logbinom + beta * self.J * N * m * m / 2.0

    def sample_equilibrium(
        self, n: int, rng: np.random.Generator, *, phi: float = 1.0, T: float = T_C_CW,
        **params: Any
    ) -> np.ndarray:
        if phi < 1.0:
            raise ValueError(
                "Curie-Weiss sample_equilibrium draws equilibrium only; use glauber_relax "
                "for off-equilibrium")
        N = int(n)
        logp = self._logP_k(N, T)
        logp -= logp.max()
        p = np.exp(logp); p /= p.sum()
        k = int(rng.choice(N + 1, p=p))            # exact number of up-spins
        s = -np.ones(N, dtype=np.int8)
        if k > 0:
            up = rng.choice(N, size=k, replace=False)
            s[up] = 1
        return s

    def initial_state(
        self, n: int, rng: np.random.Generator, *, kind: str = "ordered", **params: Any
    ) -> np.ndarray:
        if kind == "ordered":
            return np.ones(int(n), dtype=np.int8)          # cold start (all up)
        if kind == "random":
            return np.where(rng.random(int(n)) < 0.5, 1, -1).astype(np.int8)
        raise ValueError(kind)

    # ---- mean-field Glauber dynamics -------------------------------------
    def step(self, state: np.ndarray, rng: np.random.Generator, *, T: float = T_C_CW,
             **params: Any) -> np.ndarray:
        """One heat-bath sweep in the instantaneous mean field ``h = J m`` (vectorised).

        Each spin is resampled to ``+1`` with probability ``1/(1+exp(-2 beta h))``.  The
        field is refreshed once per sweep from the current magnetization -- a standard
        mean-field Glauber update.
        """
        beta = 1.0 / T
        s = state
        m = float(s.mean())
        h = self.J * m
        p_up = 1.0 / (1.0 + np.exp(-2.0 * beta * h))
        return np.where(rng.random(s.shape) < p_up, 1, -1).astype(np.int8)

    def glauber_relax_batch(
        self, n: int, b: int, rng: np.random.Generator, *, T: float = T_C_CW,
        sweeps: int, start: str = "ordered",
    ) -> np.ndarray:
        """``b`` independent chains after ``sweeps`` mean-field Glauber sweeps (off-eq)."""
        s = (np.ones((b, n), dtype=np.int8) if start == "ordered"
             else np.where(rng.random((b, n)) < 0.5, 1, -1).astype(np.int8))
        beta = 1.0 / T
        for _ in range(sweeps):
            m = s.mean(axis=1, keepdims=True)                 # per-chain mean field
            p_up = 1.0 / (1.0 + np.exp(-2.0 * beta * self.J * m))
            s = np.where(rng.random(s.shape) < p_up, 1, -1).astype(np.int8)
        return s

    def t_mix(self, n: int, **params: Any) -> int | None:
        return None

    def O_infinity(self, n: int, **params: Any) -> float | None:
        # Below T_c the symmetric ensemble has E[M]=0 (bimodal); the floor is measured
        # from Var_nu(M) directly (its exponent is a-1), not via a ratio to O_inf.
        return None
