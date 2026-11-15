# Original RSI study: experiments

The RSI paper's simulations were written as stand-alone scripts. They are preserved unchanged in `legacy-scripts/` and re-implemented in reusable, tested form in `src/rsi_shuffle`. This note says which script produced what, as far as can be determined from the scripts and the data; where it is an inference, it says so.

## Experiments

| Experiment | Model, parameters | Procedure | Original script | Reusable version | Data |
|---|---|---|---|---|---|
| **Equilibrium calibration** (RSI Fig. 2, Sec. 4.1) | Shuffle, `n = 2…4096` and validation sizes 10…3000, `T = t_mix = round(½ n ln n)`, `R = 1000` | Mean and sd of full-deck `AD` after `t_mix` against `AD_∞` | `transposition_shuffle.py` (distinct-pair chain; inferred, see `docs/known-issues.md` D1) | `rsi_shuffle.shuffle`, `rsi_shuffle.grid` (`m = n` rows) | `data/original-research/workbooks/figure2.xlsx`, `RSI Data.xlsx` |
| **Relaxation trajectories** (Figs. 3–4, Appx. Fig. 1) | `n = 512` (Fig. 3; `t` between 0 and `2 t_mix = 3194`, uneven spacing, `t/t_mix` = 0, 0.05, …, 1.1, 1.25, 1.5, 2), `n = 8…4096` (Fig. 4; `θ = t/t_mix` from 0.1 to 2.0 in steps of 0.1); `R = 1000` | Mean `AD(t)` and sd; normalised by `AD_∞` against `θ` for Fig. 4 | not identified; the `m = n` rows of `bayesian_simulation.py` compute the same quantity on a different `τ` grid | `rsi_shuffle.grid` | `figure3.xlsx`, `figure4.xlsx` (`appendix fig1.xlsx` is byte-identical to `figure4.xlsx` and is not duplicated) |
| **Tagged card** (Fig. 5) | `n = 512`, `T = t_mix`, `R = 1000` | Histogram of the final position of card 1 | `card1_position_v2.py` (and `card1_position.py` for `n = 1024`) | — | not stored; output was a figure |
| **Partial-estimator variance grid** (Figs. 6–7, S1–S5, Table 1) | `n = 2…4096`, `m` lists as in `grid.M_VALUES`, `τ = 0, 0.1, …, 2.0`, `R = 1000` | For each `(n, m, τ)`: shuffle `T = round(τ t_mix)` steps with independent pairs, draw a uniform `m`-subset without replacement, record the mean displacement of the subset; report mean, variance, sd, se, 95% interval, skewness, excess kurtosis | `bayesian_simulation.py` (no seed) | `python -m rsi_shuffle.grid` | `data/original-research/variance-grid` |
| **Single-configuration run** (Fig. 6) | `n = 512`, `m` set at the top of the file, `T` at `t_mix` | Same procedure for one `(n, m)` with histogram and fitted normal | `random_transposition_shuffle.py`, `code_runner.py` | `rsi_shuffle.grid.run_configuration` | `1.png`, `n=512, m=32, R=1000.png` (figures) |
| **First-`m` window** (discarded, Sec. 3.2, Appx. C.4) | As above but the "sample" is cards `0…m−1` | Demonstrates the bias | `legacy-scripts/first-m-window-invalidated/random_transposition_shuffle_first_m.py` | `rsi_shuffle.estimator.first_m_window_estimate` | `data/original-research/first-m-window-invalidated` |
| **Posterior and sample-size laws** (Figs. 8–9, S2, S4, S5) | `n = 512` (Fig. 8: `m = 128`, `σ₀ = 100`) | Closed-form conjugate update; `m*` from the empirical variances | no stand-alone script found | `rsi_shuffle.bayes`, `rsi_shuffle.theory` | derived from the variance grid |
| **Redrawing Fig. 2** | — | Plot of the sheet in `plot.xlsx` | `ad_plot.py` (reads a path outside the repository) | `scripts/make_original_figures.py` | `figure2.xlsx` |

Seeds: none of the legacy scripts fixes a seed (the paper says so). `rsi_shuffle.grid` and `rsi_shuffle.shuffle` accept a `numpy.random.Generator`.

## What is not here

* The script (if separate from those above) that generated `RSI Data.xlsx` sheet `1000 x 12 data` (1000 replicates of `AD` at `t_mix` for 12 sizes).
* The scripts for Figures 8, 9, S1–S5 and the AIC/BIC comparison of Figure S4. The numbers they used can be recomputed from `variance-grid`, but the original code was not found.

## Reproduction and checks

```bash
pytest tests                                       # exact enumeration and simulation tests of rsi_shuffle
python scripts/check_variance_grid.py              # new seeded runs against the stored grid and the exact law
python -m rsi_shuffle.grid --sizes 512 --taus 0 0.5 1 2 -R 1000 --seed 1 --out /tmp/n512.csv
python scripts/make_original_figures.py            # figures from stored data -> results/figures/original-research
```

`check_variance_grid.py` reruns nine `(n, m)` configurations at `τ = 1` in about two seconds.
