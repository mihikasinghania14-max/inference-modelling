# Data

Everything under `data/original-research` belongs to the Stage 1 (RSI) study. The spreadsheets are the authoritative copies; `csv/` folders hold exports of every sheet written by `scripts/export_xlsx_to_csv.py` (cell values as cached in the workbook; formulas are not re-evaluated). No data file contains personal information beyond the author's name in the workbook properties.

All data were produced by Monte Carlo simulation of the random transposition shuffle with `R = 1000` independent replicates per configuration unless stated otherwise. None of the generating scripts fixed a random seed (as the RSI paper states), so these files cannot be regenerated bit-for-bit; see `docs/reproducibility.md`. Seeded re-runs of the same experiments are in `original-research/regenerated`; they are new runs, not substitutes for the stored files.

## Stage 1 data

### `original-research/workbooks/`

| File | Sheets | What it holds | Used for |
|---|---|---|---|
| `RSI Data.xlsx` | `G1 AD vs n (theory+simulated)` | For `n = 2…4096` (and validation sizes 10…3000): `AD_∞ = (n²−1)/(3n)`, the simulated mean of `AD` after `t_mix` steps, its sd over 1000 replicates, and the relative error (a formula `(C−B)/B`) | Figure 2 |
| | `Sample trajectory shuffle n=512` | Mean `AD(t)` and sd at selected `t` up to `2 t_mix` for `n = 512` | Figure 3 |
| | `Normalised graph` | `AD(t)/AD_∞` against `t/t_mix` for four sizes, with working columns | Figure 4 (working sheet) |
| | `Theoretical values` | `n`, `AD_∞` (4 s.f.), `t_mix` | tables and figures |
| | `1000 x 12 data` | The 1000 raw replicates of `AD` at `t_mix` for each of 12 sizes (`n = 2…4096`) with `AVERAGE` formulas in rows 3 | raw data of Figure 2 |
| | `ideas` | Scratch notes (a six-card example and a to-do remark); not data | — |
| `figure2.xlsx` | `Sheet1` | The Figure 2 series, followed by a second block of "validation" sizes (`n = 10…3000`) | Figure 2 |
| `figure3.xlsx` | `Sheet1` | `t`, mean `AD(t)`, sd, `n = 512` | Figure 3 |
| `figure4.xlsx` | `n=8` … `n=4096` | `AD(t)/AD_∞` against `t/t_mix` for `θ = 0.1, 0.2, …, 2.0` | Figure 4 and Appendix Figure 1 (`appendix fig1.xlsx` in the source folder is byte-identical and not duplicated) |

**Caveat (D1).** The `n = 2` and `n = 4` entries of the Figure 2 data are consistent with a chain that swaps two *distinct* positions, not with the independent-pair chain the paper defines. See `docs/known-issues.md`.

### `original-research/variance-grid/`

| File | What it holds |
|---|---|
| `simulation_results.xlsx` | Output of `bayesian_simulation.py`: one row per `(N, M, τ)`, 1617 rows (12 values of `N`, 2 to 8 values of `M` each, 21 values of `τ = 0, 0.1, …, 2.0`), `R = 1000`. Columns: `N, M, sampling_fraction, t_mix, tau, T, mean_AD, variance, standard_deviation, standard_error, CI_lower, CI_upper, skewness, kurtosis` (`T = round(τ · t_mix)`; `CI` is the normal 95% interval for the mean; `kurtosis` is excess kurtosis). `Sheet1` and `Sheet2` hold a hand-extracted `n = 512, m = 64` trajectory over `τ`. |
| `csv/simulation_results__simulation_results.csv` | CSV export of the main sheet; used by `scripts/check_variance_grid.py` and `scripts/make_original_figures.py`. 182 cells are empty: in the 91 rows with zero variance (all 77 rows with `τ = 0`, and 14 small-`n` rows where `T = round(τ · t_mix)` is 0) skewness and kurtosis are undefined. |

Sampling scheme: uniform `m`-subset without replacement; independent pair draws. Verified: at `τ = 1` the mean agrees with `AD_∞` (mean relative bias −0.4% for `n ≥ 64`). This file was kept in the author's folder `INVALIDATED 512`, but it is **not** invalid (see `docs/known-issues.md` D2, D3).

### `original-research/first-m-window-invalidated/`

`first_m_window_by_n.xlsx` (one sheet per `n`, `m` across, rows mean, variance, sd, se, CI, coefficient of variation, skewness, kurtosis) and `first_m_window_n512.xlsx` (`n = 512`). Generated with the deterministic window "cards `0…m−1`", which the paper discards as biased (Section 3.2). Kept as evidence of the bias; not used elsewhere. The `n = 512` file reports mean 238.9 at `m = 32` against `AD_∞ = 170.7`.

## Stage 2 data

See `extension/README.md`: only the stored smoke-run results (in `results/outputs/extension-cards-smoke`) exist; the critical-system datasets were not found.

## File size and version control

The largest data file is 204 KB. Git LFS is not needed for this repository (largest tracked file about 5.8 MB, a figure PNG).

## Regenerated data

| Folder | Contents | README |
|---|---|---|
| `original-research/regenerated` | Seeded Stage 1 calibration (both pair conventions), `n = 512` trajectory and the full variance grid | `original-research/regenerated/README.md` |
| `extension/regenerated` | Wolff runs for the critical Ising model (`L = 8 … 128`, two starts each), exact Curie–Weiss moments | `extension/regenerated/README.md` |

These were generated when the repository was assembled (see `docs/provenance.md`). They are never presented as the original runs behind the papers.
