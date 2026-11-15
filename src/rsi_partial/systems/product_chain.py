"""S1 -- Product / hypercube chain on ``{0,1}^n``.

Cleanest applied case: independent coordinates, so the equilibrium (and every
off-equilibrium mixing level) has an EXACT product marginal.  Following Phase-4
sec.2(a): at mixing level ``phi`` the coordinates are independent ``Bernoulli(phi/2)``.

Exact results (all **[EST]**):
  * magnetization ``O = (1/n) sum_i x_i`` has ``E[O] = phi/2``; equilibrium ``O_inf = 1/2``.
  * ``rho2(phi) = phi / (2 - phi)``  (so ``rho2(1) = 1``).
  * ``Var(phi_hat_m) = phi(2 - phi)/m`` and susceptibility exponent ``a = 0``.

The physical dynamics :meth:`step` is independent-site (zero-coupling) Glauber, which
relaxes each coordinate to ``Bernoulli(1/2)``; it is provided for completeness and for
the autocorrelation machinery.  Because the off-equilibrium marginal is exact, detection
experiments use :meth:`sample_equilibrium` with the target ``phi`` directly.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..core.interfaces import System
from ..core.types import GroundTruth


class ProductChain(System):
    name = "product_chain"
    ground_truth = GroundTruth.EXACT_MARGINAL

    def sample_equilibrium(
        self, n: int, rng: np.random.Generator, *, phi: float = 1.0, **params: Any
    ) -> np.ndarray:
        p = 0.5 * phi
        return (rng.random(n) < p).astype(np.int8)

    def initial_state(
        self, n: int, rng: np.random.Generator, *, kind: str = "ordered", **params: Any
    ) -> np.ndarray:
        if kind == "ordered":
            return np.zeros(n, dtype=np.int8)      # all-zero ordered state (phi = 0)
        if kind == "random":
            return self.sample_equilibrium(n, rng, phi=1.0)
        raise ValueError(kind)

    def step(
        self,
        state: np.ndarray,
        rng: np.random.Generator,
        *,
        flip_rate: float = 1.0,
        **params: Any,
    ) -> np.ndarray:
        """Independent-site Glauber: each site is resampled to Bernoulli(1/2) w.p.
        ``flip_rate`` (heat-bath update with zero coupling)."""
        out = state.copy()
        resample = rng.random(state.shape[0]) < flip_rate
        out[resample] = (rng.random(int(resample.sum())) < 0.5).astype(np.int8)
        return out

    def true_phi(self, t: int, n: int, *, flip_rate: float = 1.0, **params: Any) -> float | None:
        # Starting from all-zeros, P(x_i=1 at time t) = 0.5 (1 - (1-flip_rate)^t);
        # phi = E[O]/O_inf = 2 * P(x_i=1) = 1 - (1-flip_rate)^t.
        return 1.0 - (1.0 - flip_rate) ** t

    def O_infinity(self, n: int, **params: Any) -> float | None:
        return 0.5

    def rho2(self, phi: float, **params: Any) -> float | None:
        return phi / (2.0 - phi)

    def susceptibility_exponent_a(self, **params: Any) -> float | None:
        return 0.0
