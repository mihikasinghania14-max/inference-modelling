"""Experiment engine: turn (system, observable, grid) into replicate-pool measurements.

Responsibilities:
  * Generate ``R`` replicate states (batched, so peak memory is ``O(batch * n)``, never
    ``O(R n)``), through the ``System``/``Observable`` interface -- with vectorised fast
    paths for the exact systems (chain, cards) keyed by name.
  * Apply the veil at each ``m`` (a prefix of a fresh per-replicate random permutation is
    a valid uniform ``m``-subset) and record the estimator pool.
  * Emit per-``m`` pools: ``Ohat``, sampling-SD estimate, plus the state-level ``O(x)`` and
    within-config variance ``V_x`` (ddof=0, matching RSI App. E) needed for the identity,
    floor, and SNR diagnostics.

The engine never sees a theoretical value.  It reads only ``m`` components per replicate
per ``m`` (access audit holds); reusing a batch's states across the ``m``-grid is a
common-random-numbers optimisation that leaves each ``m``'s within-pool replicates
independent.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.interfaces import Observable, System
from ..observe.veil import reveal_batch

DEFAULT_BATCH = 512


# ---------------------------------------------------------------------------
# Vectorised D-array generators (fast paths) -- return a (b, n) float array of D_i
# ---------------------------------------------------------------------------


def _cards_D_equilibrium(n: int, b: int, rng: np.random.Generator) -> np.ndarray:
    perms = np.argsort(rng.random((b, n)), axis=1)        # b uniform permutations
    idx = np.arange(n)
    return np.abs(perms - idx).astype(np.float64)


def _cards_fixedpoint_equilibrium(n: int, b: int, rng: np.random.Generator) -> np.ndarray:
    perms = np.argsort(rng.random((b, n)), axis=1)
    idx = np.arange(n)
    return (perms == idx).astype(np.float64)


def _cards_D_relaxed(n: int, b: int, phi: float, rng: np.random.Generator,
                     *, observable_name: str) -> np.ndarray:
    """Relax ``b`` decks to mixing level ``phi`` via vectorised random transpositions
    (physical dynamics), then read the observable's per-component contribution."""
    from ..systems.cards import Cards
    t = Cards().steps_for_phi(n, phi)
    P = np.tile(np.arange(n, dtype=np.int64), (b, 1))
    rows = np.arange(b)
    for _ in range(t):
        a = rng.integers(n, size=b)
        c = rng.integers(n, size=b)
        Pa = P[rows, a].copy()
        P[rows, a] = P[rows, c]
        P[rows, c] = Pa
    idx = np.arange(n)
    if observable_name == "fixed_point_fraction":
        return (P == idx).astype(np.float64)
    return np.abs(P - idx).astype(np.float64)


def _chain_D(n: int, b: int, phi: float, rng: np.random.Generator) -> np.ndarray:
    p = 0.5 * phi
    return (rng.random((b, n)) < p).astype(np.float64)


def _ising_D_from_states(states, observable, b, n) -> np.ndarray:
    D = np.empty((b, n), dtype=np.float64)
    for r in range(b):
        D[r] = observable.per_component(states[r])
    return D


def _vectorised_D_batch(
    system: System, observable: Observable, n: int, b: int, phi: float,
    rng: np.random.Generator, *, T: float | None,
    sweeps: int | None = None, start: str = "random",
) -> np.ndarray | None:
    """Return a (b, n) D-array via a fast path, or ``None`` if none applies."""
    sname, oname = system.name, observable.name
    if sname == "cards":
        if phi >= 1.0:
            return _cards_fixedpoint_equilibrium(n, b, rng) if oname == "fixed_point_fraction" \
                else _cards_D_equilibrium(n, b, rng)
        return _cards_D_relaxed(n, b, phi, rng, observable_name=oname)
    if sname == "product_chain":
        return _chain_D(n, b, phi, rng)
    if sname == "ising2d":
        TT = T if T is not None else getattr(system, "T_c", None)
        if sweeps is not None:
            # OFF-EQUILIBRIUM: physical checkerboard-Glauber relaxation (the certified process).
            states = system.glauber_relax_batch(n, b, rng, T=TT, sweeps=sweeps, start=start)  # type: ignore[attr-defined]
            return _ising_D_from_states(states, observable, b, n)
        if phi >= 1.0:
            # EQUILIBRIUM: warm one Wolff chain per batch, read b decorrelated configs.
            states = system.wolff_stream(n, rng, T=TT, count=b)  # type: ignore[attr-defined]
            return _ising_D_from_states(list(states), observable, b, n)
        raise ValueError(
            "off-equilibrium Ising needs sweeps=... (Glauber relaxation); phi<1 alone is "
            "undefined for Ising -- pass sweeps and start to measure_ensemble")
    return None


def _generic_D_batch(
    system: System, observable: Observable, n: int, b: int, phi: float,
    rng: np.random.Generator, *, T: float | None,
) -> np.ndarray:
    """Fallback: generate each state via the interface (works for any system, incl. Ising)."""
    D = np.empty((b, n), dtype=np.float64)
    kw = {} if T is None else {"T": T}
    for r in range(b):
        state = system.sample_equilibrium(n, rng, phi=phi, **kw)
        D[r] = observable.per_component(state)
    return D


# ---------------------------------------------------------------------------
# Pool measurement
# ---------------------------------------------------------------------------


@dataclass
class MPool:
    """Estimator pool at a single sample size ``m``."""
    m: int
    n: int
    Ohat: np.ndarray            # (R,) estimator values
    sd: np.ndarray              # (R,) per-replicate sampling-SD estimates (ddof=1, FPC)


@dataclass
class EnsembleResult:
    """Everything measured for one (system, observable, n, phi) ensemble."""
    n: int
    phi: float
    O_full: np.ndarray          # (R,) population values O(x)
    Vx: np.ndarray              # (R,) within-config population variances V_x (ddof=0)
    pools: dict[int, MPool] = field(default_factory=dict)
    reads_per_replicate: dict[int, int] = field(default_factory=dict)  # access audit


def measure_ensemble(
    system: System,
    observable: Observable,
    n: int,
    R: int,
    rng: np.random.Generator,
    m_grid: list[int],
    *,
    phi: float = 1.0,
    T: float | None = None,
    sweeps: int | None = None,
    start: str = "random",
    batch: int = DEFAULT_BATCH,
    force_generic: bool = False,
) -> EnsembleResult:
    """Generate ``R`` replicate states and measure estimator pools at every ``m`` in
    ``m_grid`` (each ``m <= n``).  Batched to keep memory ``O(batch * n)``.

    Equilibrium is the default (``phi=1``).  For **off-equilibrium Ising**, pass
    ``sweeps=k`` (and optionally ``start='random'|'ordered'``): the ensemble is produced by
    ``k`` checkerboard-Glauber sweeps from that start (the physical dynamics), and the
    resulting mixing level is measured, not assumed.  For cards, off-equilibrium is selected
    by ``phi<1`` (relaxation by transpositions).

    ``force_generic`` bypasses the vectorised fast paths (used by tests to cross-check the
    fast path against the interface path).
    """
    m_grid = sorted({int(m) for m in m_grid if 1 <= m <= n})
    O_full = np.empty(R, dtype=np.float64)
    Vx = np.empty(R, dtype=np.float64)
    Ohat = {m: np.empty(R, dtype=np.float64) for m in m_grid}
    sd = {m: np.empty(R, dtype=np.float64) for m in m_grid}
    reads: dict[int, int] = {}              # genuine per-replicate read counts from the veil

    done = 0
    while done < R:
        b = min(batch, R - done)
        if force_generic:
            D = _generic_D_batch(system, observable, n, b, phi, rng, T=T)
        else:
            D = _vectorised_D_batch(system, observable, n, b, phi, rng, T=T,
                                    sweeps=sweeps, start=start)
            if D is None:
                D = _generic_D_batch(system, observable, n, b, phi, rng, T=T)
        # state-level quantities (ground-truth diagnostics -- NOT fed to inference)
        O_full[done:done + b] = D.mean(axis=1)
        Vx[done:done + b] = D.var(axis=1, ddof=0)          # V_x = (1/n) sum (D_i - O)^2
        # ALL inference values come through the veil: the engine never indexes D itself
        # for an estimator (the audited path is the executed path, critique sec.3.1).
        perm = np.argsort(rng.random((b, n)), axis=1)      # one draw reused as nested prefixes
        for m in m_grid:
            vals, reads_m = reveal_batch(D, m, rng, perm=perm)
            reads[m] = reads_m                             # genuine count = vals.shape[1]
            Ohat[m][done:done + b] = vals.mean(axis=1)
            if m >= 2:
                s2 = vals.var(axis=1, ddof=1)
                fpc = (n - m) / (n - 1) if n > 1 else 0.0
                sd[m][done:done + b] = np.sqrt(np.maximum(s2 / m * fpc, 0.0))
            else:
                sd[m][done:done + b] = 0.0
        done += b

    # The audit is now a real post-condition, not an assertion by fiat: the veil reported
    # exactly `m` reads for every m (census m=n reveals all n, the Ô_n=O consistency check).
    for m in m_grid:
        assert reads[m] == m, f"veil exposed {reads[m]} != {m} components at m={m}"

    pools = {m: MPool(m=m, n=n, Ohat=Ohat[m], sd=sd[m]) for m in m_grid}
    return EnsembleResult(n=n, phi=phi, O_full=O_full, Vx=Vx, pools=pools,
                          reads_per_replicate=reads)
