# rsi_partial: partial-observation equilibrium-certification simulator

A modular, reproducible, memory-efficient research instrument for the programme
*"Partial Observation of Random Systems."* (Phase 5 of the project plan.) It certifies the equilibrium **mean of a
chosen intensive additive observable** from a uniform random `m`-subset of a system's `n`
components, and measures every quantity the Phase-4 predictions (P1–P8) require.

The science is specified in `research/research-extension/theoretical-framework/Phase4_conjecture_and_design.md`, building on the RSI paper
(`research/original-rsi-research/paper/`) and the Phase 3 pilots
(`Phase3_findings.md`, `Phase3_interpretation.md` in the same folder). This package is *execution*, not redesign. The design of
the package is in `research/research-extension/supplementary-materials/Phase5_architecture.md` and its running, epistemically
tagged results log is `RESULTS.md` in the same folder. Where a docstring refers to the "design brief", it means the build specification that
accompanied the start of Phase 5 (a task description for writing the instrument, not a research document; it is not included in this
repository); its scientific content is `Phase4_conjecture_and_design.md`. "Phase-4 sec." and "modelling note Mk" refer to
`Phase4_conjecture_and_design.md` and `Phase5_architecture.md`.

## The spine

For any intensive additive observable estimated from a uniform `m`-subset without
replacement, the law of total variance gives, exactly,

```
Var(Ô_m) = (1/m)·(n−m)/(n−1)·E_x[V_x]  +  Var_x(O)      (identity ★)
CV²(Ô_m) ≈ 1/(ρ²·m)  +  R_n            (m ≪ n)
```

a `1/m` sampling term plus an `m`-independent floor `R_n = Var_ν(O)/O_∞² = χ_n/(nO_∞²)`.
Size-free certification holds iff the per-component SNR `ρ²` is bounded **and** the
observable susceptibility is summable (`a=0`, where `χ_n ≍ n^a`).

## Layout

```
rsi_partial/
  core/         interfaces (System, Observable), rng/seed protocol, the Observation veil
  systems/      S1 product_chain, S2 cards, S3 ising2d (Wolff GT + Glauber dynamics)
  observables/  AD, fixed_point_fraction (S2′), magnetization, |M|/energy/block (Ising)
  observe/      veil.py — memory-aware O(Rm) subsampler + access audit
  inference/    estimator, bayes (conjugate Gaussian), detect (power/Type-I), mstar
  measure/      diagnostics — floor, χ_n, R_n, ρ²(ϕ), skew/kurtosis, τ_int, I_eff, self-calib
  analysis/     predictions (theory, ISOLATED), bootstrap CIs, score (verdicts), figures
  experiments/  engine (pool measurement), detection, registry (frozen grids), run (CLI)
  storage/      schema + Parquet/CSV io + JSON manifests
  tests/        analytical · veil · identity · ising_gt · regression_phase3 (the gate) · guardrails
```

## Guardrails (enforced structurally)

- **The veil is physical.** Inference sees only the `m` revealed values; `Observation` has
  no state field. An access audit asserts exactly `m` reads (`tests/test_veil.py`).
- **No prediction is hard-coded into an estimator.** All theory lives in
  `analysis/predictions.py`; a test asserts no measuring layer imports it.
- **Ground truth ≠ dynamics.** `sample_equilibrium` (rapid: Wolff/exact) is distinct from
  `step` (the certified physical dynamics: Glauber/transposition).
- **Measure → score → plot.** Figures come from stored tables, never live runs.

## Quickstart

From the repository root:

```bash
pip install -e ".[test]"
pytest src/rsi_partial/tests                # 74 tests: analytical, veil, identity, Ising ground truth, regression gate, guardrails
python -m rsi_partial.experiments.run snr_floor  --system cards --sizes 512 2048 8192 --smoke --out /tmp/res
python -m rsi_partial.experiments.run detection --system cards --sizes 2048 8192 --deltas 0.15 0.10 --smoke --out /tmp/res
bash scripts/reproduce_extension_cards_smoke.sh   # the two runs above, compared with the stored results
```

`--smoke` uses `R = 300` replicates per configuration (full runs: `R = 2000` for variance configurations and `R = 5000` for
detection, see `experiments/registry.py`). The default root seed is `20260818`.

**Gate rule:** the analytical and Phase-3 regression suites must be green before any run. They are (74 tests when this
repository was assembled).

## What is implemented and what is not

* Implemented end-to-end: `snr_floor` and `detection` experiments for the card shuffle (`cards`) and the product chain
  (`product_chain`), observables `AD`, `fixed_point_fraction`, `magnetization`.
* Implemented as components with tests but with **no campaign runner**: the 2D Ising systems (`systems/ising2d.py`: Wolff
  ground-truth sampler, checkerboard Glauber dynamics) and observables (`observables/ising_obs.py`), and the Curie–Weiss system
  (`systems/curie_weiss.py`). The Ising critical-floor measurements reported in the extension paper were not produced by a
  runner in this package that could be found (`docs/known-issues.md` X4).
* `registry.py` lists the planned grids for all stages; only Stage 1 was run (smoke scale).
