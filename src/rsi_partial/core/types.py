"""Core data types, including the physical veil boundary.

The single most important object here is :class:`Observation`.  It is the *only*
thing the inference layer is allowed to see.  It carries exactly ``m`` revealed
per-component contributions and metadata -- **it deliberately has no field for the
full state**, so full-system leakage into inference is impossible by construction
(design brief Design Invariant 1; Phase-4 sec.5.5).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import numpy as np


class GroundTruth(str, Enum):
    """How an equilibrium ensemble for a config is generated."""

    EXACT_MARGINAL = "exact_marginal"     # product-form / direct draw (S1, S2, S6 density)
    MCMC_CLUSTER = "mcmc_cluster"          # Wolff / Swendsen-Wang (S3)
    TRANSFER_MATRIX = "transfer_matrix"    # exact 1-D chain (S5)
    METROPOLIS = "metropolis"              # mean-field exact/Metropolis (S4)


class Regime(str, Enum):
    EQUILIBRIUM = "equilibrium"
    OFF_EQUILIBRIUM = "off_equilibrium"


@dataclass(frozen=True)
class Observation:
    """The veil.  Inference sees this and nothing else.

    Parameters
    ----------
    values:
        Exactly ``m`` revealed per-component contributions ``D_i`` for ``i`` in the
        uniform random subset ``S``.  This is the whole payload.
    n, m:
        System size and sample size.
    seed_hex:
        128-bit child seed (hex) that produced this observation, for provenance.
    config_hash:
        Hash of the originating configuration.
    reads:
        Audited number of component reads used to build ``values`` (must equal ``m``);
        see the access audit in :mod:`rsi_partial.observe.veil`.
    """

    values: np.ndarray
    n: int
    m: int
    seed_hex: str = ""
    config_hash: str = ""
    reads: int = -1

    # NOTE: there is intentionally no `state` field.  Do not add one.

    def __post_init__(self) -> None:
        if self.values.shape[0] != self.m:
            raise ValueError(
                f"Observation carries {self.values.shape[0]} values but m={self.m}"
            )
        if self.m > self.n:
            raise ValueError(f"m={self.m} exceeds n={self.n}")


@dataclass(frozen=True)
class ConfigKey:
    """Canonical identity of a single experiment configuration (pre-replicate).

    Hashing this (see :func:`rsi_partial.core.rng.config_hash`) yields deterministic
    filenames and the per-replicate seed-spawn key.
    """

    experiment_id: str
    system: str
    observable: str
    ground_truth: str
    n: int
    regime: str = Regime.EQUILIBRIUM.value
    phi_target: float = 1.0
    # optional / system-specific (kept in a sorted dict for a stable hash)
    params: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        d = {
            "experiment_id": self.experiment_id,
            "system": self.system,
            "observable": self.observable,
            "ground_truth": self.ground_truth,
            "n": self.n,
            "regime": self.regime,
            "phi_target": self.phi_target,
        }
        # flatten params with a prefix so the hash is stable and collision-free
        for k in sorted(self.params):
            d[f"param.{k}"] = self.params[k]
        return d
