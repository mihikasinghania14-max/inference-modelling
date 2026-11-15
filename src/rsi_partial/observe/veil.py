"""The veil: reveal a uniform random ``m``-subset (without replacement) of the
per-component contributions, and nothing else.

Memory strategy (Phase-4 sec.8.4; design brief sec.3):
  * The naive ``argsort(random(R, n))`` subset draw is banned -- it is ``O(Rn)`` and
    OOM-ed at ``n ~ 1e5`` in Phase 3.
  * :func:`reveal` draws ``m`` indices without replacement in ``O(m)`` and gathers
    ``values[idx]``.  Callers loop one replicate at a time, so no ``R x n`` array is
    ever materialised: the sampling hot path is ``O(Rm)``.
  * :func:`reveal_exchangeable_first_m` is the ``O(m)`` fast path for observables that
    are exchangeable under ``nu`` (the first ``m`` freshly generated components are a
    valid uniform subset).  Guarded by an explicit flag and cross-checked in tests.

Access audit: every :class:`~rsi_partial.core.types.Observation` records ``reads`` =
number of component reads; the runner asserts ``reads == m`` per replicate
(control Phase-4 sec.5.5).
"""

from __future__ import annotations

import numpy as np

from ..core.types import Observation


def reveal(
    values_full: np.ndarray,
    m: int,
    rng: np.random.Generator,
    *,
    seed_hex: str = "",
    config_hash: str = "",
) -> Observation:
    """Reveal ``m`` of the ``n`` precomputed contributions, uniformly without replacement.

    ``values_full`` is the observable's full per-component array (length ``n``).  It is
    computed once at generation time; the veil selects ``m`` of its entries.  The full
    array is not stored on the returned :class:`Observation`.
    """
    n = int(values_full.shape[0])
    if m > n:
        raise ValueError(f"m={m} exceeds n={n}")
    if m == n:
        idx = np.arange(n)
    else:
        idx = rng.choice(n, size=m, replace=False)
    picked = np.asarray(values_full[idx], dtype=np.float64)
    return Observation(
        values=picked, n=n, m=m, seed_hex=seed_hex, config_hash=config_hash, reads=int(m)
    )


def reveal_batch(
    D_batch: np.ndarray,
    m: int,
    rng: np.random.Generator,
    *,
    perm: np.ndarray | None = None,
) -> tuple[np.ndarray, int]:
    """Vectorised veil for a whole batch: reveal exactly ``m`` components per replicate
    (row of ``D_batch``), uniformly without replacement.

    This is **the** selection path the engine uses, so the audited path is the executed
    path (critique sec.3.1).  Returns ``(values, reads_per_replicate)`` where ``values`` is
    ``(b, m)`` and ``reads_per_replicate`` is **derived from the returned array width**
    (``values.shape[1]``) -- a genuine count of components exposed to inference, not a
    hard-coded assertion.  Inference must consume only ``values``; the caller may keep the
    full ``D_batch`` for ground-truth diagnostics but must never feed it to an estimator.

    ``perm`` (an optional ``(b, n)`` array of per-row permutations) lets the caller reuse
    one draw across an ``m``-grid via nested prefixes ``perm[:, :m]`` (common random
    numbers).  Note (critique sec.2.5): nested prefixes make different ``m`` correlated --
    fine for a monotone ``m*`` crossing, but variance-vs-``m`` fits should treat the points
    as dependent.  Pass ``perm=None`` for an independent draw at this ``m``.
    """
    b, n = D_batch.shape
    if m > n:
        raise ValueError(f"m={m} exceeds n={n}")
    if m == n:
        vals = np.array(D_batch, copy=True)       # census: reveals all n (the Ô_n=O check)
    else:
        if perm is None:
            perm = np.argsort(rng.random((b, n)), axis=1)
        vals = np.take_along_axis(D_batch, perm[:, :m], axis=1)
    reads_per_replicate = int(vals.shape[1])       # genuine: what inference actually sees
    return vals, reads_per_replicate


def reveal_exchangeable_first_m(
    generate_m,
    n: int,
    m: int,
    rng: np.random.Generator,
    *,
    seed_hex: str = "",
    config_hash: str = "",
) -> Observation:
    """Exchangeable fast path: ``generate_m(m, rng)`` returns exactly ``m`` per-component
    contributions drawn from the equilibrium marginal.

    Valid ONLY when the observable's components are exchangeable under ``nu`` (every
    component has common inclusion probability ``m/n``; e.g. product chain, cards ``AD``
    at equilibrium).  Never use for non-exchangeable systems (Ising) or off-equilibrium
    states where position matters.
    """
    picked = np.asarray(generate_m(m, rng), dtype=np.float64)
    if picked.shape[0] != m:
        raise ValueError(f"generate_m returned {picked.shape[0]} values, expected {m}")
    return Observation(
        values=picked, n=n, m=m, seed_hex=seed_hex, config_hash=config_hash, reads=int(m)
    )


def relabel(values_full: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Randomly permute component labels (control Phase-4 sec.5.5): results must be
    invariant, precluding positional leakage."""
    return values_full[rng.permutation(values_full.shape[0])]
