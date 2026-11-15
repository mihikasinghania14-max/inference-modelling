"""Partial-observation estimators of the average displacement.

``uniform_sample_estimate`` is the estimator of RSI Definition 2.9: a uniform random
``m``-subset of the cards, drawn without replacement, and the mean of their
displacements.  It is unbiased for ``AD(pi)`` (Lemma 2.10).

``first_m_window_estimate`` is the deterministic window ("the first ``m`` cards") that the
paper discusses and *discards* (Section 3.2 and Appendix C.4) because it is biased by the
position dependence of ``E|sigma(i) - i|`` (Corollary 2.7).  It is provided so that the
bias can be demonstrated; it must not be used for inference.
"""

from __future__ import annotations

import numpy as np


def uniform_sample_estimate(d: np.ndarray, m: int, rng: np.random.Generator) -> np.ndarray | float:
    """Mean of ``m`` displacements drawn uniformly without replacement.

    ``d`` is ``(n,)`` (one deck) or ``(R, n)`` (``R`` decks; an independent subset is drawn
    for each row).
    """
    d = np.asarray(d)
    n = d.shape[-1]
    if not 1 <= m <= n:
        raise ValueError(f"need 1 <= m <= n, got m={m}, n={n}")
    if d.ndim == 1:
        idx = rng.choice(n, size=m, replace=False)
        return float(d[idx].mean())
    # Rows: the m smallest of n i.i.d. uniforms give a uniform m-subset (O(R n) memory).
    idx = np.argpartition(rng.random(d.shape), m - 1, axis=1)[:, :m] if m < n else \
        np.tile(np.arange(n), (d.shape[0], 1))
    return np.take_along_axis(d, idx, axis=1).mean(axis=1)


def first_m_window_estimate(d: np.ndarray, m: int) -> np.ndarray | float:
    """Mean displacement of cards ``0..m-1`` (biased; see module docstring)."""
    d = np.asarray(d)
    n = d.shape[-1]
    if not 1 <= m <= n:
        raise ValueError(f"need 1 <= m <= n, got m={m}, n={n}")
    return d[..., :m].mean(axis=-1)
