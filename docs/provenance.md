# Provenance of the files in this repository

This note says which material is the author's original work product, which was reconstructed or regenerated when the repository was assembled, and which parts were prepared with software assistance. It exists so that a reader can tell stored results from re-runs.

## 1. Original research material (as found)

`research/` holds the RSI paper (PDF), the extension paper (PDF and its 15 September LaTeX source), the integrated manuscript, and the working notes of the extension project, copied without edits except where a stray file name was redacted. `data/original-research/workbooks` and `data/original-research/variance-grid` hold the spreadsheets behind the RSI figures. `experiments/original-research/legacy-scripts` holds the author's scripts as written. `results/outputs/extension-cards-smoke` is the stored output of the measurement instrument's smoke run.

## 2. Computations added when the repository was assembled

| Location | What | Status |
|---|---|---|
| `experiments/extension/verification`, `results/outputs/extension-verification` | Independent checks of the extension's algebra and several numerical claims | Checks, with stored output |
| `experiments/extension/campaign`, `data/extension/regenerated` | New seeded Wolff runs for the critical Ising model, exact Curie–Weiss moments | Regenerated, not the original runs |
| `scripts/regenerate_stage1_data.py`, `data/original-research/regenerated` | Seeded Stage 1 datasets (calibration for both pair conventions, `n = 512` trajectory, full grid) | Regenerated, not the original runs |
| `src/rsi_shuffle` | A clean, tested implementation of the Stage 1 simulation and theory | Code written for this repository, following the papers and the author's scripts |
| `tests/`, `notebooks/`, `scripts/make_*_figures.py` | Tests, walkthroughs and figures | Written for this repository |

None of the regenerated files replaces an original. Where a regenerated value differs from a number in a paper, the difference is recorded in `docs/known-issues.md` and the paper is not changed.

## 3. Disclosure of AI assistance

An AI assistant, Claude (Anthropic), was used in this project under the author's direction. Its role, reconstructed from the working notes, file dates and audit documents available when the repository was assembled (the session transcripts themselves are not part of the repository, so the author should check each item):

* **Extension phase (August–October 2026).** Helped to design and write the `rsi_partial` measurement instrument and its tests, to run its computations, and to draft the extension paper and working notes (`research/research-extension`), from the author's specification.
* **Verification and audit (October 2026).** Wrote the independent verification scripts in `experiments/extension/verification` and an audit of the integrated manuscript.
* **Integrated manuscript (November 2026).** Helped to prepare the manuscript draft in the *Bayesian Analysis* format.
* **This repository (November 2026).** Wrote `src/rsi_shuffle`, the tests, notebooks, figure and regeneration scripts, the Ising campaign code and the documentation; ran the regenerated computations; and compiled the inventory and `docs/known-issues.md`.

No AI assistance for the Stage 1 (RSI) paper, the RSI data or the scripts under `experiments/original-research` is recorded in the available material; the author should state whether any was used. The mathematical content, the interpretation of the results and the decision to publish are the author's responsibility. The assistant's computations were re-run and checked against exact results where these exist (`docs/known-issues.md` records what did not agree).

*For the author to confirm, correct or extend before the repository is made public. A journal's own policy on AI-use disclosure may require different wording in the manuscript.*

## 4. Not included

Files that were found but deliberately not copied (private working documents, task descriptions, application material) are listed in `docs/inventory.md`.
