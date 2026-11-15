"""Bootstrap confidence intervals for headline quantities (>=1000 resamples).

Every Phase-5 headline number is reported with a bootstrap 95% CI (Phase-4 sec.4.3).
Generic :func:`bootstrap_ci` plus specialisations for the quantities that are not simple
means: a variance, a ratio of pooled quantities, and a log-log slope (for the floor
exponent ``a-1`` and the ``m*`` exponent ``b``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class CI:
    value: float
    lo: float
    hi: float
    level: float = 0.95

    def excludes(self, x: float) -> bool:
        return x < self.lo or x > self.hi

    def contains(self, x: float) -> bool:
        return self.lo <= x <= self.hi


def bootstrap_ci(
    data: np.ndarray,
    statistic: Callable[[np.ndarray], float],
    *,
    n_boot: int = 1000,
    level: float = 0.95,
    rng: np.random.Generator | None = None,
) -> CI:
    """Percentile bootstrap CI for ``statistic(data)`` by resampling rows with replacement."""
    rng = rng or np.random.default_rng()
    data = np.asarray(data)
    N = data.shape[0]
    point = float(statistic(data))
    boot = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, N, size=N)
        boot[b] = statistic(data[idx])
    a = (1.0 - level) / 2.0
    lo, hi = np.quantile(boot, [a, 1.0 - a])
    return CI(value=point, lo=float(lo), hi=float(hi), level=level)


def bootstrap_variance_ci(pool: np.ndarray, *, n_boot=1000, level=0.95,
                          rng=None) -> CI:
    """CI on ``Var(pool)`` (e.g. the floor ``Var_nu(O)`` or ``Var(Ohat_m)``)."""
    return bootstrap_ci(pool, lambda x: float(np.var(x, ddof=1)),
                        n_boot=n_boot, level=level, rng=rng)


def bootstrap_ratio_ci(num_pool: np.ndarray, den_pool: np.ndarray,
                       stat=np.mean, *, n_boot=1000, level=0.95, rng=None) -> CI:
    """CI on ``stat(num)/stat(den)`` resampling the two pools jointly (paired if equal len)."""
    rng = rng or np.random.default_rng()
    num = np.asarray(num_pool); den = np.asarray(den_pool)
    point = float(stat(num) / stat(den))
    Nn, Nd = num.shape[0], den.shape[0]
    boot = np.empty(n_boot)
    for b in range(n_boot):
        boot[b] = float(stat(num[rng.integers(0, Nn, Nn)]) / stat(den[rng.integers(0, Nd, Nd)]))
    a = (1 - level) / 2
    lo, hi = np.quantile(boot, [a, 1 - a])
    return CI(point, float(lo), float(hi), level)


def loglog_slope_ci(
    x: np.ndarray,
    y_samples: list[np.ndarray],
    *,
    divisor: np.ndarray | None = None,
    n_boot: int = 1000,
    level: float = 0.95,
    rng: np.random.Generator | None = None,
) -> CI:
    """CI on the slope of ``log(Var(pool)/divisor) vs log(x)``.

    ``x``: the size grid (e.g. ``n`` values).  ``y_samples``: for each ``x`` a replicate
    pool whose variance is the numerator (e.g. per-``n`` pools of ``O_full`` -> Var_nu(O)).
    ``divisor``: per-``x`` scale applied to the variance before the log (pass ``O_inf**2``
    to fit the *relative floor* ``R_n`` and hence the exponent ``a-1`` -- essential when
    ``O_inf`` is not ``Theta(1)``, e.g. cards where ``O_inf ~ n``; modelling note M1).
    The slope is refitted on each bootstrap resample of every pool.
    """
    rng = rng or np.random.default_rng()
    logx = np.log(np.asarray(x, dtype=np.float64))
    div = np.ones(len(x)) if divisor is None else np.asarray(divisor, dtype=np.float64)

    def _slope(pools: list[np.ndarray]) -> float:
        y = np.array([np.var(p, ddof=1) for p in pools]) / div
        return float(np.polyfit(logx, np.log(y), 1)[0])

    point = _slope(y_samples)
    boot = np.empty(n_boot)
    for b in range(n_boot):
        res = [p[rng.integers(0, p.shape[0], p.shape[0])] for p in y_samples]
        boot[b] = _slope(res)
    a = (1 - level) / 2
    lo, hi = np.quantile(boot, [a, 1 - a])
    return CI(point, float(lo), float(hi), level)
