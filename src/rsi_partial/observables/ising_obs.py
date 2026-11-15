"""Observables on the 2-D Ising system (S3).

State is an ``(L, L)`` array of spins ``+/-1`` (periodic BC), ``n = L^2``.

  * ``abs_magnetization``: per-component readout is the **signed** spin ``s_i``; the
    population value is ``|mean(s)|`` (modelling note M5 -- the intensive additive
    observable is the signed mean spin, with ``|.|`` applied at the population level and
    the +/-M sign symmetry folded in analysis).
  * ``energy_density``: ``D_i = -s_i (s_right + s_down)`` so ``sum_i D_i`` equals the
    total bond energy (each bond counted once on a periodic lattice); population value is
    the energy per site.
  * ``block_mean``: ``D_i`` = mean spin of the ``b x b`` block containing site ``i``
    (coarse-grained observable for the CX4 comparison).

Equilibrium values are analytic only in special cases (Onsager); the observables return
``None`` for metadata and the simulator measures.  Onsager references live in
:mod:`rsi_partial.analysis.predictions`.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..core.interfaces import Observable


class AbsMagnetization(Observable):
    name = "abs_magnetization"

    def per_component(self, state: np.ndarray) -> np.ndarray:
        return state.astype(np.float64).ravel()

    def population_value(self, state: np.ndarray) -> float:
        return float(abs(np.mean(state)))


class EnergyDensity(Observable):
    name = "energy_density"

    def per_component(self, state: np.ndarray) -> np.ndarray:
        s = state.astype(np.float64)
        right = np.roll(s, -1, axis=1)
        down = np.roll(s, -1, axis=0)
        d = -s * (right + down)          # sum over sites = total bond energy (J=1)
        return d.ravel()

    def population_value(self, state: np.ndarray) -> float:
        return float(np.mean(self.per_component(state)))


class BlockMean(Observable):
    name = "block_mean"

    def __init__(self, block: int = 4) -> None:
        self.block = block
        self.name = f"block_mean_b{block}"

    def per_component(self, state: np.ndarray) -> np.ndarray:
        s = state.astype(np.float64)
        L = s.shape[0]
        b = self.block
        if L % b != 0:
            raise ValueError(f"block size {b} does not divide L={L}")
        # block-average each site's value: reshape into (L/b, b, L/b, b), mean over blocks
        nb = L // b
        blocks = s.reshape(nb, b, nb, b).mean(axis=(1, 3))          # (nb, nb) block means
        expanded = np.repeat(np.repeat(blocks, b, axis=0), b, axis=1)  # back to (L, L)
        return expanded.ravel()

    def population_value(self, state: np.ndarray) -> float:
        return float(np.mean(state))
