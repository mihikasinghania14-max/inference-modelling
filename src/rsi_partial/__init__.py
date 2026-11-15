"""rsi_partial — Partial-observation equilibrium-certification research simulator.

Phase-5 instrument for the RSI research programme "Partial Observation of Random
Systems".  Layered, reproducible, memory-efficient.  See ``Phase5_architecture.md``.

Epistemic guardrails (enforced structurally):
  * The inference layer receives only the ``m`` revealed contributions (the "veil").
  * No theoretical prediction is hard-coded into an estimator; predictions live in
    ``rsi_partial.analysis.predictions`` and are used only for scoring.
  * Ground-truth generation (rapid sampler) is distinct from the physical dynamics.
"""

__version__ = "0.1.0"
