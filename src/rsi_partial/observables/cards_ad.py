"""Average-displacement observable ``AD`` on the cards system (RSI's observable).

``D_i = |pi[i] - i|`` (0-indexed; identical in law to the paper's 1-indexed
``|sigma(i) - i|``).  ``O = AD = (1/n) sum_i D_i``.

Note (modelling note M1): ``D_i`` grows with ``n`` and ``AD_inf = Theta(n)``, so ``AD``
is not a bounded / ``Theta(1)`` observable.  The framework's scale-free claims concern
only the dimensionless ratios (``rho^2``, ``R_n``, ``CV``); the code works in those and
never assumes bounded-by-constant contributions.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..core.interfaces import Observable


class AverageDisplacement(Observable):
    name = "AD"

    def per_component(self, state: np.ndarray) -> np.ndarray:
        idx = np.arange(state.shape[0], dtype=np.int64)
        return np.abs(state - idx).astype(np.float64)

    def E_infinity(self, n: int, **params: Any) -> float | None:
        return (n * n - 1) / (3.0 * n)

    def rho2(self, phi: float, **params: Any) -> float | None:
        # Exact equilibrium value only (=2); off-equilibrium is empirical -> analysis.
        return 2.0 if phi >= 1.0 else None

    def a(self, **params: Any) -> float | None:
        return 0.0
