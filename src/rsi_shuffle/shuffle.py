"""The random transposition shuffle and the average displacement observable.

Conventions (as in the RSI paper)
---------------------------------
* Cards are labelled ``0..n-1`` here (the paper uses ``1..n``); card ``i`` starts at
  position ``i``.
* ``perm[i]`` is the position of card ``i`` after the shuffle, i.e. ``pi_t(i)``.
* ``D_i = |pi(i) - i|`` is the displacement of card ``i`` and
  ``AD(pi) = (1/n) * sum_i D_i``.

Pair convention
---------------
Definition 2.1 and Appendix C.2 of the RSI paper draw the two positions
``a_t, b_t`` independently and uniformly from ``{1..n}``, so that ``a_t = b_t``
(a self-loop) occurs with probability ``1/n``.  This is ``pair_convention="independent"``
and is the default.

Some of the early scripts (``transposition_shuffle.py``, ``card1_position.py`` and
``transposition_sim.py`` in the original script folder) instead draw two *distinct*
positions with ``random.sample(range(n), 2)``.  That is ``pair_convention="distinct"``.
The two conventions define different chains (the distinct-pair chain has no self-loops
and is periodic for n = 2); they agree to leading order in ``1/n``.  See
``docs/known-issues.md`` (item D1) for which datasets were generated with which
convention.
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np

PairConvention = Literal["independent", "distinct"]


def mixing_time(n: int, rounding: Literal["nearest", "ceil"] = "nearest") -> int:
    """Leading-order mixing time ``(1/2) n ln n`` of RSI Eq. (2.1).

    ``rounding="nearest"`` is the convention stated in the paper and used by
    ``bayesian_simulation.py``; ``"ceil"`` is used by ``random_transposition_shuffle.py``.
    The two differ by at most one step.
    """
    if n < 2:
        raise ValueError("n must be at least 2")
    x = 0.5 * n * math.log(n)
    return int(round(x)) if rounding == "nearest" else int(math.ceil(x))


def _draw_pairs(rng: np.random.Generator, n: int, size: int | tuple[int, ...],
                pair_convention: PairConvention) -> tuple[np.ndarray, np.ndarray]:
    a = rng.integers(0, n, size=size)
    if pair_convention == "independent":
        b = rng.integers(0, n, size=size)
    elif pair_convention == "distinct":
        if n < 2:
            raise ValueError("distinct pairs need n >= 2")
        # b is uniform on the n-1 positions different from a.
        b = (a + 1 + rng.integers(0, n - 1, size=size)) % n
    else:  # pragma: no cover
        raise ValueError(f"unknown pair_convention {pair_convention!r}")
    return a, b


def shuffle_deck(n: int, t: int, rng: np.random.Generator,
                 pair_convention: PairConvention = "independent") -> np.ndarray:
    """Apply ``t`` random transpositions to the ordered deck; return ``perm`` with
    ``perm[card] = final position``.

    Single-deck reference implementation (pure Python loop over the ``t`` swaps).
    """
    deck = np.arange(n, dtype=np.int64)          # deck[position] = card
    if t > 0:
        a, b = _draw_pairs(rng, n, t, pair_convention)
        for i, j in zip(a.tolist(), b.tolist()):
            deck[i], deck[j] = deck[j], deck[i]
    perm = np.empty(n, dtype=np.int64)
    perm[deck] = np.arange(n, dtype=np.int64)
    return perm


def shuffle_decks(n: int, t: int, R: int, rng: np.random.Generator,
                  pair_convention: PairConvention = "independent") -> np.ndarray:
    """``R`` independent shuffles at once, vectorised across replicates.

    Returns an ``(R, n)`` array ``perm`` with ``perm[r, card] = final position`` in
    replicate ``r``.  Distributionally identical to calling :func:`shuffle_deck` ``R``
    times; memory is ``O(R n)``.
    """
    deck = np.tile(np.arange(n, dtype=np.int32), (R, 1))     # deck[r, position] = card
    rows = np.arange(R)
    for _ in range(t):
        a, b = _draw_pairs(rng, n, R, pair_convention)
        ca = deck[rows, a].copy()
        deck[rows, a] = deck[rows, b]
        deck[rows, b] = ca
    perm = np.empty_like(deck)
    np.put_along_axis(perm, deck.astype(np.int64), np.broadcast_to(
        np.arange(n, dtype=np.int32), (R, n)), axis=1)
    return perm


def displacements(perm: np.ndarray) -> np.ndarray:
    """Per-card displacements ``|pi(i) - i|`` (works on ``(n,)`` or ``(R, n)`` arrays)."""
    perm = np.asarray(perm)
    return np.abs(perm - np.arange(perm.shape[-1]))


def average_displacement(perm: np.ndarray) -> np.ndarray | float:
    """``AD(pi) = (1/n) sum_i |pi(i) - i|`` (RSI Definition 2.4); vectorised over rows."""
    return displacements(perm).mean(axis=-1)
