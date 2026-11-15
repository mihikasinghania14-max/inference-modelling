"""S2' -- Fixed-point-fraction observable on the cards system.

``D_i = 1{pi[i] == i}`` (indicator that card ``i`` is a fixed point).
``O = F = (1/n) sum_i D_i`` = fraction of fixed points.

This is the **degenerate-SNR control** (Phase-4 S2'): the number of fixed points of a
uniform permutation is asymptotically Poisson(1), so ``E[F] = 1/n -> 0`` and the
per-component SNR ``rho^2 = E[D]^2 / Var(D) ~ 1/n -> 0`` (condition (i) fails).  It
exercises the ``rho^2``/``R_n`` machinery on an ``O_inf -> 0`` case from day one
(design brief Stage-0 requirement; modelling note M3).

Nothing here divides by ``O_inf``; the estimator works from raw moments and the analysis
layer annotates the degenerate regime.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..core.interfaces import Observable


class FixedPointFraction(Observable):
    name = "fixed_point_fraction"

    def per_component(self, state: np.ndarray) -> np.ndarray:
        idx = np.arange(state.shape[0], dtype=np.int64)
        return (state == idx).astype(np.float64)

    def E_infinity(self, n: int, **params: Any) -> float | None:
        # E[F] = (1/n) * n * P(pi[i]=i) = 1/n  (exact for a uniform permutation).
        return 1.0 / n

    def rho2(self, phi: float, **params: Any) -> float | None:
        # Degenerate: rho^2 = E[D]^2/Var(D) with D ~ Bernoulli(1/n)
        #   = (1/n)^2 / [ (1/n)(1 - 1/n) ] = 1/(n-1)  ->  ~ 1/n.
        # Return None: this is the control whose SNR the simulator MEASURES.
        return None

    def a(self, **params: Any) -> float | None:
        # Correlations are summable (condition (ii) not binding); size-limitation here
        # is from rho^2 -> 0, not from a correlation floor.
        return 0.0
