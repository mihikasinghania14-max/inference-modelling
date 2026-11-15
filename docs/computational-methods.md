# Computational methods

## Packages

| Package | Role |
|---|---|
| `src/rsi_shuffle` | Stage 1. `shuffle` (chain, both pair conventions, vectorised over replicates), `estimator` (uniform and first-`m` estimators), `bayes` (conjugate posterior), `theory` (closed forms with exact rational arithmetic), `grid` (the `(n, m, τ)` grid and its CLI). |
| `src/rsi_partial` | Stage 2 measurement instrument (see its README). Layers: `core` (interfaces, seed protocol), `systems` (cards, product chain, 2D Ising, Curie–Weiss), `observables`, `observe` (the veil), `inference` (estimator, Bayes, detection, `m*`), `measure` (floor, `χ_n`, `ρ²`, diagnostics), `analysis` (predictions, bootstrap, scoring, figures), `experiments` (engine, registry, runner), `storage`. |

## Algorithms

* **Shuffle.** `shuffle_decks` stores a deck array of shape `(R, n)` and swaps `R` pairs at once per step with fancy indexing; the final permutation is obtained by inverting the deck. Cost `O(t · R)` plus `O(R n)` memory. `shuffle_deck` is a single-deck pure-Python reference used in tests to confirm the batch version has the same distribution.
* **Uniform `m`-subset per replicate.** The indices of the `m` smallest of `n` i.i.d. uniforms (`argpartition`) form a uniform subset; each row gets its own.
* **Pair convention.** `independent`: `a, b` uniform on `[n]`. `distinct`: `b = (a + 1 + U) mod n` with `U` uniform on `[n−1]`, which is uniform on the positions different from `a`.
* **Exact theory.** `theory` returns `fractions.Fraction` values on request so tests can compare closed forms with enumeration without rounding.
* **Wolff cluster sampler** (`wolff.c`, `rsi_partial.systems.ising2d`): single-cluster updates with `p_add = 1 − e^{−2K}` at `K = ½ ln(1+√2)`. In `wolff.c` the number of cluster flips between measurements is a fixed integer `k = ⌈n/⟨C⟩⌉` estimated from a warm-up; the schedule therefore does not depend on the cluster sizes realised in the run.
* **Exact Ising references.** Enumeration of all `2^16` configurations (`L = 4`), transfer matrix (`L = 6, 8`), Kaufman's partition function with mpmath (`L` up to 128).

## Randomness and seeds

* `rsi_shuffle` functions take a `numpy.random.Generator`; the CLI takes `--seed` (omit for OS entropy, as in the original scripts). Tests use fixed seeds.
* `rsi_partial` derives all randomness from a root seed (default `20260818`) as described in `docs/methodology.md`. Runs are bit-reproducible given the seed and library versions; the stored smoke results were reproduced exactly on a different platform and library versions (see `docs/reproducibility.md`).

## Dependencies

NumPy ≥ 1.24, SciPy ≥ 1.10, pandas ≥ 2.0, PyArrow ≥ 12, Matplotlib ≥ 3.7, openpyxl ≥ 3.1 (spreadsheet export); pytest ≥ 7 for tests; SymPy ≥ 1.12 and mpmath ≥ 1.3 for `shuffle_algebra.py` and `ising_exact.py`; a C99 compiler for `wolff.c`. Python ≥ 3.11 (the stored runs used Python 3.14; this repository was assembled and tested with Python 3.13).

## Performance (assembled on a 1-core cloud container)

| Task | Time |
|---|---|
| `pytest` (103 tests) | about 80 s |
| `scripts/check_variance_grid.py` | about 2 s |
| `scripts/reproduce_extension_cards_smoke.sh` | about 7 s |
| `experiments/extension/verification/run_all.sh` | about 25 s |
| `scripts/regenerate_stage1_data.py` | about 5 min (2 cores unused: single process) |
| `experiments/extension/campaign/run_ising_campaign.sh` | hours (`L = 96, 128` dominate); `analyse_ising_campaign.py` about 30 s |
| `curie_weiss_exact.py` | seconds |

## Critical-Ising campaign (added with the regenerated data)

`experiments/extension/campaign/ising_wolff.c` is a single-cluster Wolff sampler on the `L × L` torus at `T_c`. Between two measurements a fixed number `k = ⌈n/⟨C⟩⌉` of clusters is flipped, `⟨C⟩` being estimated in a warm-up phase, so the measurement times do not depend on the cluster sizes realised afterwards. `analyse_ising_campaign.py` computes the floors from block means with a jackknife, applies the two-start gate and fits slopes against `log n` for `L ≥ 16` with a block bootstrap. Protocol and results: `data/extension/regenerated/README.md`.
