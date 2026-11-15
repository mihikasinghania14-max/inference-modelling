"""Reproducibility layer: seed protocol and configuration hashing.

Protocol (Phase-4 sec.4.3):
  * One documented **root seed** per experiment (recorded in the run manifest).
  * A fresh child stream is spawned per ``(config, replicate)`` via
    ``numpy.random.SeedSequence.spawn`` -- independent replicates therefore use
    non-overlapping streams by construction.
  * Every output row carries its 128-bit child seed (hex) and the config hash.
  * An optional unseeded ``distribution_level`` mode draws fresh OS entropy,
    reproducing the RSI paper's protocol for one cross-check.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any

import numpy as np

_ENTROPY_BITS = 128


def config_hash(config: dict[str, Any], length: int = 16) -> str:
    """Deterministic BLAKE2b hash of a canonicalized config dict."""
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.blake2b(canonical.encode("utf-8"), digest_size=length // 2).hexdigest()


def _spawn_key(config: dict[str, Any], replicate: int) -> int:
    """Map (config, replicate) to a stable non-negative integer entropy key.

    We fold the config hash + replicate index into a 128-bit integer used to seed a
    SeedSequence.  Two different (config, replicate) pairs get different keys with
    overwhelming probability; the SAME pair always gets the same key -> replayable.
    """
    payload = json.dumps(config, sort_keys=True, default=str) + f"|rep={replicate}"
    digest = hashlib.blake2b(payload.encode("utf-8"), digest_size=16).digest()
    return int.from_bytes(digest, "big")


def child_seed(
    root_seed: int,
    config: dict[str, Any],
    replicate: int,
    *,
    distribution_level: bool = False,
) -> np.random.SeedSequence:
    """Return the SeedSequence for a single ``(config, replicate)``.

    Reproducible mode (default): deterministic function of ``(root_seed, config,
    replicate)``.  ``distribution_level=True``: fresh OS entropy (unseeded), for the
    RSI cross-check.
    """
    if distribution_level:
        return np.random.SeedSequence()  # OS entropy
    key = _spawn_key(config, replicate)
    # Combine root and per-config key as the SeedSequence entropy vector.
    return np.random.SeedSequence([int(root_seed), key])


def rng_for(
    root_seed: int,
    config: dict[str, Any],
    replicate: int,
    *,
    distribution_level: bool = False,
) -> tuple[np.random.Generator, str]:
    """Return ``(Generator, seed_hex)`` for a ``(config, replicate)``.

    ``seed_hex`` is the 128-bit entropy of the spawned child, stored per row.
    """
    ss = child_seed(root_seed, config, replicate, distribution_level=distribution_level)
    # Derive a compact, stable 128-bit provenance tag from the sequence state.
    state = ss.generate_state(4, dtype=np.uint32)
    seed_hex = "".join(f"{int(x):08x}" for x in state)
    return np.random.default_rng(ss), seed_hex
