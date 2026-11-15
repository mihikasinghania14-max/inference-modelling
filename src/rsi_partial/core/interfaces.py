"""Abstract interfaces: ``System`` and ``Observable``.

New systems/observables register by implementing these ABCs; no other layer needs
to change (the decoupling the design brief requires).  Metadata methods return ``None``
where a quantity is not analytically known -- **never a faked value** (epistemic
guardrail: the simulator measures; theory lives in the analysis layer).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

import numpy as np

from .types import GroundTruth


class System(ABC):
    """A finite Markov chain with a distinct rapid equilibrium sampler and dynamics."""

    name: str
    ground_truth: GroundTruth

    # ---- generation -------------------------------------------------------
    @abstractmethod
    def sample_equilibrium(
        self, n: int, rng: np.random.Generator, *, phi: float = 1.0, **params: Any
    ) -> np.ndarray:
        """Return one equilibrium (or, if ``phi<1`` and exact, off-equilibrium) state.

        Uses a *correct rapid* sampler, kept distinct from :meth:`step`.  For systems
        with an exact off-equilibrium marginal (S1), ``phi`` selects the mixing level;
        otherwise off-equilibrium states are produced by running :meth:`step` from an
        ordered start (see :meth:`relax_to_phi`).
        """

    @abstractmethod
    def step(self, state: np.ndarray, rng: np.random.Generator, **params: Any) -> np.ndarray:
        """Apply one step of the *physical dynamics* whose convergence is certified."""

    def initial_state(
        self, n: int, rng: np.random.Generator, *, kind: str = "ordered", **params: Any
    ) -> np.ndarray:
        """Return an initial condition for the dynamics (``'ordered'`` or ``'random'``)."""
        raise NotImplementedError

    def relax_to_phi(
        self,
        n: int,
        rng: np.random.Generator,
        *,
        phi: float,
        **params: Any,
    ) -> np.ndarray:
        """Off-equilibrium state at mixing level ``phi`` via the physical dynamics.

        Default implementation runs :meth:`step` from an ordered start for the number
        of steps implied by :meth:`true_phi`.  Systems with an exact off-equilibrium
        marginal override :meth:`sample_equilibrium` with ``phi`` instead.
        """
        raise NotImplementedError

    # ---- analytic metadata (None where unknown) --------------------------
    def true_phi(self, t: int, n: int, **params: Any) -> float | None:
        return None

    def t_mix(self, n: int, **params: Any) -> int | None:
        return None

    def O_infinity(self, n: int, **params: Any) -> float | None:
        return None

    def rho2(self, phi: float, **params: Any) -> float | None:
        return None

    def susceptibility_exponent_a(self, **params: Any) -> float | None:
        return None


class Observable(ABC):
    """An intensive additive observable ``O(x) = (1/n) sum_i D_i(x)``."""

    name: str

    @abstractmethod
    def per_component(self, state: np.ndarray) -> np.ndarray:
        """Return the length-``n`` array of per-component contributions ``D_i``.

        This is the ONLY function the veil samples from.
        """

    def population_value(self, state: np.ndarray) -> float:
        """Full-system observable ``O(x)`` (ground-truth diagnostic; NOT fed to inference)."""
        return float(np.mean(self.per_component(state)))

    # ---- analytic metadata (None where unknown) --------------------------
    def E_infinity(self, n: int, **params: Any) -> float | None:
        return None

    def rho2(self, phi: float, **params: Any) -> float | None:
        return None

    def a(self, **params: Any) -> float | None:
        return None
