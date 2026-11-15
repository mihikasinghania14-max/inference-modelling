# Extension data: what exists and what is missing

## Exists (stored in `results/outputs`)

| Dataset | Location | How it was generated | Reproduce |
|---|---|---|---|
| Stage 1 smoke results of the measurement instrument, random transposition shuffle (`cards`), observable `AD`: equilibrium SNR `ρ²`, floor `Var_ν(O)`, `R_n`, `χ_n` and `I_eff(m)` at `n = 512, 2048, 8192` (`stage1_snr_floor__*.parquet`, 25 summary rows) and detection cost `m*(δ, n)` at `n = 2048, 8192`, `δ = 0.15, 0.10` (`stage1_detection__all.parquet`) | `results/outputs/extension-cards-smoke/` with JSON manifests (root seed `20260818`, `R = 300`, library versions, configuration hash) | `python -m rsi_partial.experiments.run … --smoke` | `bash scripts/reproduce_extension_cards_smoke.sh` reproduces the numeric columns bit-for-bit (verified when this repository was assembled) |
| Independent verification outputs | `results/outputs/extension-verification/` | `experiments/extension/verification` | `bash experiments/extension/verification/run_all.sh` (identical output, verified) |

## Not found

The following datasets are described in the extension paper, in the working notes or in the integrated manuscript, but **no raw output, summary table or generating script was found** in any location accessible when this repository was assembled (searched: all connected folders of the author's machine, the attached project, uploaded files, and the unpacked source archives of the manuscripts). They have **not** been reconstructed or replaced by re-runs presented as the originals.

| Missing item | Described in |
|---|---|
| 2D Ising critical runs, `L ≤ 128` (Wolff sampler, two-start equilibration gate, bootstrap-CI slopes): energy and `|M|` floors, slopes `−0.82` and `+0.02`, `a[M] = 0.889`, floor `0.066`; and the earlier run to `L = 80` summarised in `Phase5_criticality_confirmed.md` | extension paper, Section 5 and Figure 2; `Phase5_criticality_confirmed.md` |
| Curie–Weiss runs (slopes `−0.996` above and `+0.005` below `T_c`), CX1 | extension paper, Section 5; `RESULTS.md` (the instrument's log records the result but no stored result file was found) |
| Correlated-field experiments (reflecting boundaries; coverage from 91% to about 12%) | project notes; extension draft |
| Algorithm 1 validation runs (single-realisation floor estimator) | extension paper, Sections 6–7 (proposed, unvalidated); integrated manuscript, Section 8 |
| Full `snr_floor` and `detection` runs at full `R` (`R = 2000` and `5000`) and the full Phase-5 P1–P8 campaign | `src/rsi_partial/experiments/registry.py`; `RESULTS.md` lists them as unrun |
| Scripts producing Figures 1–2 of the extension paper (`fig_decomposition.png`, `fig_criticality.png`) | `research/research-extension/paper/source` holds the figures only |

Independent re-computations of several of the corresponding quantities (exact small-torus Ising floors, Kaufman energy floors, Wolff `L = 16, 32`, Curie–Weiss limits, a torus field experiment) are in `results/outputs/extension-verification`. They are checks, not the original data.

## Regenerated (new runs, not the originals)

`data/extension/regenerated` holds new seeded computations of two of the missing items: the critical 2D Ising floors for `L = 8 … 128` and the Curie–Weiss moments. They are documented, with protocol, seeds and results, in `regenerated/README.md`. The Ising regeneration gives an `|M|` floor of about 0.072 and an energy slope of −0.863, which agree with the exact values and with the integrated manuscript and differ from the extension paper's 0.066 and −0.82 (`docs/known-issues.md`, X1, X2, R1). The correlated-field, Algorithm 1 and full-campaign items above are still missing.
