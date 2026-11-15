"""Equilibration gate for locally-generated (MCMC) ensembles (Phase-4 sec.4.3, sec.5.2).

The diagnostics existed (`measure.diagnostics.equilibration_check`, `integrated_autocorr_time`)
but nothing forced a generator to prove stationarity before recording (critique sec.1.3-1.4).
This module wires them into a usable gate:

  * :func:`two_start_check` -- generate the observable from an ORDERED and a RANDOM start
    under the physical dynamics and require agreement (the two-initial-condition test);
    a deliberately under-burned run must FAIL this (the negative control).
  * :func:`tau_int_of_dynamics` -- integrated autocorrelation time of the observable along a
    single Glauber chain, to *set* replicate spacing rather than assert it.
  * :func:`gated_equilibrium` -- draw an equilibrium ensemble only if the two-start check
    passes at the requested burn-in; otherwise raise so a biased ensemble is never recorded.
"""

from __future__ import annotations

import numpy as np

from ..core.interfaces import Observable, System
from ..measure.diagnostics import (
    EquilibrationCheck, equilibration_check, integrated_autocorr_time)


def _relax_pool(system: System, observable: Observable, n: int, R: int,
                rng: np.random.Generator, *, T: float, sweeps: int, start: str) -> np.ndarray:
    """Observable population values over ``R`` chains after ``sweeps`` sweeps from ``start``."""
    states = system.glauber_relax_batch(n, R, rng, T=T, sweeps=sweeps, start=start)  # type: ignore[attr-defined]
    return np.array([observable.population_value(states[r]) for r in range(states.shape[0])])


def two_start_check(
    system: System, observable: Observable, n: int, R: int, rng: np.random.Generator,
    *, T: float, sweeps: int, tol_sigma: float = 3.0,
) -> EquilibrationCheck:
    """Two-initial-condition agreement after ``sweeps`` of the physical dynamics.

    Ordered (cold) and random (hot) starts must relax to the same observable mean within
    ``tol_sigma`` standard errors; if ``sweeps`` is too small they disagree (the certificate
    detects under-burn-in instead of being fooled by it)."""
    a = _relax_pool(system, observable, n, R, rng, T=T, sweeps=sweeps, start="ordered")
    b = _relax_pool(system, observable, n, R, rng, T=T, sweeps=sweeps, start="random")
    return equilibration_check(a, b, tol_sigma=tol_sigma)


def tau_int_of_dynamics(
    system: System, observable: Observable, n: int, rng: np.random.Generator,
    *, T: float, n_sweeps: int = 2000, burn: int = 200,
) -> float:
    """Integrated autocorrelation time of the observable along one Glauber chain (in sweeps),
    used to choose replicate spacing (>= 2 tau_int)."""
    s = system.initial_state(n, rng, kind="random")
    for _ in range(burn):
        s = system.step(s, rng, T=T)
    trace = np.empty(n_sweeps)
    for t in range(n_sweeps):
        s = system.step(s, rng, T=T)
        trace[t] = observable.population_value(s)
    return integrated_autocorr_time(trace)


def gated_equilibrium(
    system: System, observable: Observable, n: int, R: int, rng: np.random.Generator,
    *, T: float, burn_sweeps: int, check_R: int = 200, tol_sigma: float = 3.0,
) -> np.ndarray:
    """Return an equilibrium observable pool ONLY if the two-start check passes at
    ``burn_sweeps``; otherwise raise ``RuntimeError`` (never silently record a biased set).

    This is the campaign guard for locally-generated Ising ground truth when the rapid
    Wolff sampler is not used: it refuses to proceed until stationarity is demonstrated.
    """
    chk = two_start_check(system, observable, n, check_R, rng, T=T,
                          sweeps=burn_sweeps, tol_sigma=tol_sigma)
    if not chk.passed:
        raise RuntimeError(
            f"equilibration gate FAILED at burn={burn_sweeps} sweeps "
            f"(ordered {chk.ordered_mean:.4f} vs random {chk.random_mean:.4f}); "
            f"increase burn-in or use the Wolff sampler")
    return _relax_pool(system, observable, n, R, rng, T=T, sweeps=burn_sweeps, start="random")
