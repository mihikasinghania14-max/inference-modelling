"""Magnetization observable on the product chain (S1).

``D_i = x_i in {0,1}``; ``O = (1/n) sum_i x_i``.  Exact: ``E_inf = 1/2``,
``rho^2(phi) = phi/(2 - phi)``, ``a = 0``.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..core.interfaces import Observable


class ChainMagnetization(Observable):
    name = "magnetization"

    def per_component(self, state: np.ndarray) -> np.ndarray:
        return state.astype(np.float64)

    def E_infinity(self, n: int, **params: Any) -> float | None:
        return 0.5

    def rho2(self, phi: float, **params: Any) -> float | None:
        return phi / (2.0 - phi)

    def a(self, **params: Any) -> float | None:
        return 0.0
