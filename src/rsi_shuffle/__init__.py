"""rsi_shuffle: code for the original RSI study (random transposition shuffle).

This package re-implements, in reusable form, the computations described in
*Bayesian Inference in Random Transposition Shuffles* (RSI India, 2026):

* :mod:`rsi_shuffle.shuffle`    -- the random transposition chain and the average
  displacement observable;
* :mod:`rsi_shuffle.estimator`  -- the uniform random-sample estimator and the
  (biased, discarded) first-m window;
* :mod:`rsi_shuffle.bayes`      -- the conjugate Gaussian posterior of Section 5.3;
* :mod:`rsi_shuffle.theory`     -- closed-form quantities (equilibrium moments,
  variance law, sample-size law, and the exact footrule variance used by the
  extension);
* :mod:`rsi_shuffle.grid`       -- the (n, m, tau) Monte Carlo grid behind the
  variance-law figures.

The original scripts are preserved unchanged in
``experiments/original-research/legacy-scripts``.
"""

from . import bayes, estimator, grid, shuffle, theory

__all__ = ["bayes", "estimator", "grid", "shuffle", "theory"]
