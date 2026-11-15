# Reproducibility

Status of each result class when this repository was assembled (Python 3.13, NumPy/SciPy as pinned in `requirements.txt` ranges, Linux).

| Result | Command | Outcome |
|---|---|---|
| Closed forms for the shuffle (`AD_∞`, `Var F`, `E[V]`, exact estimator variance, mean law at `n = 4`) | `pytest tests/test_rsi_shuffle_theory.py` | pass, by exact enumeration of `S_n` (`n ≤ 7`) |
| Simulation code of Stage 1 | `pytest tests/test_rsi_shuffle_simulation.py` | pass; statistical tests use fixed seeds and 4.5-sigma tolerances |
| Instrument test suite (`rsi_partial`) | `pytest src/rsi_partial/tests` | 74 tests pass |
| Stored Stage-1 smoke results of the instrument | `bash scripts/reproduce_extension_cards_smoke.sh` | numeric columns **identical bit-for-bit** to the stored copy, although the stored run used Windows, Python 3.14, NumPy 2.5.0 and the reproduction used Linux, Python 3.13 |
| Verification scripts and the two Wolff runs | `bash experiments/extension/verification/run_all.sh` | all seven outputs identical to the stored outputs |
| Stage 1 variance grid | `python scripts/check_variance_grid.py` | statistically consistent with the stored grid and with the exact law (9 configurations, passes) |
| Stage 1 figures | `python scripts/make_original_figures.py` | regenerated from stored data |
| Stage 1 data (`RSI Data.xlsx`, figure sheets, the grid) | not regenerable bit-for-bit | generated without a fixed seed, as stated in the RSI paper |
| Extension: Ising `L ≤ 128`, Curie–Weiss, field, Algorithm 1 runs behind the paper's numbers | not reproducible from this repository | raw data and scripts were not found (`data/extension/README.md`) |
| Seeded Stage 1 datasets (`data/original-research/regenerated`) | `python scripts/regenerate_stage1_data.py` (about 5 min); the grid file reproduced byte for byte on the assembling platform | new runs, not the original files |
| Critical Ising regeneration (`data/extension/regenerated`) | `bash experiments/extension/campaign/run_ising_campaign.sh …` then `analyse_ising_campaign.py` (hours); `L = 8` block file reproduced byte for byte | new runs, not the original files |
| Exact Curie–Weiss moments | `python experiments/extension/campaign/curie_weiss_exact.py OUT.csv` (seconds) | exact summation, no sampling |
| Notebooks | `make notebooks` | executed copies are stored with their outputs |

## Seeds

* Stored smoke results: root seed `20260818`, `R = 300` (manifests in `results/outputs/extension-cards-smoke`).
* Wolff runs: `./wolff L 200000 1` for `L = 16, 32`. `shuffle_sim.py`: `default_rng(12345)`. `field.py`: `default_rng(2026)`.
* `rsi_shuffle` tests: fixed integer seeds in each test.
* Stage 1 regeneration: root seed `20261115` (calibration), `+1` (trajectory), `+2` (grid). Ising campaign: `1000 + L` (random start), `2000 + L` (ordered start); bootstrap generator seed `20261115`. Notebooks: `default_rng(2026)`.

## Comparing stored and new runs

`scripts/compare_with_stored.py STORED NEW` compares Parquet outputs: numeric columns must be equal, `None` and `NaN` encodings are treated as the same missing value (a difference that appears only in non-numeric or null columns of the Parquet files), and manifests are compared on experiment id, seed and configuration only, because host, platform, library versions and timestamps legitimately differ.

## Environment

`pyproject.toml` (editable install), `requirements.txt`, `environment.yml`. The manifests of the stored runs record the library versions used (NumPy 2.5.0, SciPy 1.18.0, pandas 3.0.5, PyArrow 25.0.1, Python 3.14.5).
