r"""``m*`` finders -- two distinct certification-cost notions (modelling note M2).

  * **detection-m\***: smallest ``m`` giving power ``>= 1-beta`` at size ``alpha`` for a
    deviation ``delta``.  Phase-5 headline (P3-P5).  Target ``kappa^2/(rho^2 delta^2)``.
  * **estimation-m\***: smallest ``m`` whose posterior credible half-width is
    ``<= epsilon * O_inf`` at level ``1-alpha``.  Target ``z_{1-alpha/2}^2/(rho^2 eps^2)``
    (RSI Eq. 5.7).

Both reduce to locating the crossing of a MONOTONE (in ``m``) metric against a target,
which :func:`find_threshold_m` does on a geometric grid + integer bisection to a relative
tolerance (default ``+/-5%``).  The metric callable is supplied by the experiment layer,
so this module stays agnostic to how pools are generated (decoupling).  No theory value
is used here; predicted targets live in :mod:`rsi_partial.analysis.predictions`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Literal

import numpy as np


def geometric_grid(m_min: int, m_max: int, *, base: int = 2) -> list[int]:
    """Geometric grid ``{m_min, base*m_min, ...} <= m_max`` (Phase-4 sec.4.2)."""
    grid: list[int] = []
    m = m_min
    while m <= m_max:
        grid.append(int(m))
        m *= base
    if grid and grid[-1] != m_max:
        grid.append(int(m_max))
    return grid


@dataclass
class MStarResult:
    mstar: int | None                 # None if no m in range meets the target
    metric_at_mstar: float
    target: float
    direction: str                    # 'increasing' | 'decreasing'
    tol: float
    evaluations: dict[int, float] = field(default_factory=dict)
    bracket: tuple[int, int] | None = None
    infeasible: bool = False          # target unreachable within [m_min, m_max]


def find_threshold_m(
    evaluate: Callable[[int], float],
    target: float,
    *,
    direction: Literal["increasing", "decreasing"],
    m_min: int,
    m_max: int,
    tol: float = 0.05,
    base: int = 2,
) -> MStarResult:
    """Smallest integer ``m`` in ``[m_min, m_max]`` at which ``evaluate(m)`` meets ``target``.

    ``direction='increasing'``: metric rises with ``m`` (e.g. power); we seek the first
    ``m`` with ``metric >= target``.  ``direction='decreasing'``: metric falls with ``m``
    (e.g. relative credible half-width); we seek the first ``m`` with ``metric <= target``.

    Strategy: scan a geometric grid to bracket the crossing, then integer-bisect the
    bracket until its width is within a relative tolerance ``tol``.  ``evaluate`` may be
    stochastic (Monte-Carlo power); callers should use enough replicates that its noise
    near the target is small compared with ``tol``.
    """
    evals: dict[int, float] = {}

    def _meets(m: int) -> bool:
        if m not in evals:
            evals[m] = float(evaluate(m))
        v = evals[m]
        return v >= target if direction == "increasing" else v <= target

    grid = geometric_grid(m_min, m_max, base=base)

    # Bracket on the grid.
    lo: int | None = None    # largest grid m that does NOT meet target
    hi: int | None = None    # smallest grid m that DOES meet target
    for m in grid:
        if _meets(m):
            hi = m
            break
        lo = m

    if hi is None:
        # Target never met within range -> infeasible / needs larger m.
        best_m = grid[-1]
        return MStarResult(
            mstar=None, metric_at_mstar=evals[best_m], target=target,
            direction=direction, tol=tol, evaluations=evals,
            bracket=None, infeasible=True,
        )

    if lo is None:
        # Even the smallest m meets the target: m* <= m_min (report m_min).
        return MStarResult(
            mstar=hi, metric_at_mstar=evals[hi], target=target, direction=direction,
            tol=tol, evaluations=evals, bracket=(m_min, hi), infeasible=False,
        )

    # Integer bisection between lo (fails) and hi (meets).
    while hi - lo > 1 and (hi - lo) / lo > tol:
        mid = (lo + hi) // 2
        if _meets(mid):
            hi = mid
        else:
            lo = mid

    return MStarResult(
        mstar=hi, metric_at_mstar=evals[hi], target=target, direction=direction,
        tol=tol, evaluations=evals, bracket=(lo, hi), infeasible=False,
    )


def detection_mstar(
    power_of_m: Callable[[int], float],
    *,
    power_target: float = 0.80,
    m_min: int = 4,
    m_max: int = 100_000,
    tol: float = 0.05,
) -> MStarResult:
    """Smallest ``m`` with detection power ``>= power_target`` (default ``1-beta=0.80``)."""
    return find_threshold_m(
        power_of_m, power_target, direction="increasing",
        m_min=m_min, m_max=m_max, tol=tol,
    )


def estimation_mstar(
    rel_halfwidth_of_m: Callable[[int], float],
    epsilon: float,
    *,
    m_min: int = 4,
    m_max: int = 100_000,
    tol: float = 0.05,
) -> MStarResult:
    """Smallest ``m`` with posterior relative credible half-width ``<= epsilon``."""
    return find_threshold_m(
        rel_halfwidth_of_m, epsilon, direction="decreasing",
        m_min=m_min, m_max=m_max, tol=tol,
    )
