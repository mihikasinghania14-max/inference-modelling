# Seeded Stage 1 datasets

The Stage 1 data in `../workbooks` and `../variance-grid` were produced without a fixed seed, as the RSI paper states, and cannot be rebuilt bit-for-bit. The three files here are **seeded counterparts** of the same experiments. They are new runs, not the originals, and the figures of the paper were not drawn from them.

Produced by `scripts/regenerate_stage1_data.py` with root seed `20261115`, `R = 1000` replicates, `t_mix = round(½ n ln n)`. The complete script takes about five minutes.

| File | Counterpart of | Columns / contents |
|---|---|---|
| `calibration_pair_conventions.csv` | RSI Figure 2 (workbook `figure2.xlsx`, sheet `G1`) | For `n = 2 … 4096` and the validation sizes `10 … 3000`: `AD_∞`, and the mean and sd of `AD` after `t_mix` steps under **both** pair conventions (`independent`: the chain of RSI Definition 2.1; `distinct`: two different positions at every step). |
| `trajectory_n512.csv` | RSI Figure 3 (`figure3.xlsx`) | `n = 512`: mean and sd of `AD(t)` at `t = 0, 0.05 t_mix, …, 2 t_mix`, with `AD/AD_∞` and the closed-form mean law `φ(t) = 1 − (1 − 2/n)^t`. |
| `variance_grid_seeded.csv` | RSI Figures 6–7 (`../variance-grid/simulation_results.xlsx`) | The same 1617 `(N, M, τ)` rows and columns as the stored grid, generated with the original independent-replicate procedure and seed `20261117` (root seed + 2). |

## What the data show

* **Conventions.** At `n = 2` the distinct-pair chain gives mean `AD` exactly 1.0 with sd 0 (the value in the Figure 2 data); the independent-pair chain gives 0.502 with sd 0.500 (exact: 0.5 and 0.5). From `n = 64` upward the two conventions agree within Monte Carlo error. See `docs/known-issues.md`, D1 and R2.
* **Variance law.** At `τ = 1`, `n ≥ 64`, `m < n` (46 configurations) the ratio of simulated variance to the exact law (`rsi_shuffle.theory.estimator_variance`) has mean 1.002 and median 1.000, range 0.83–1.10, as expected from `R = 1000`. Against the leading-order `n²/(18m)` the ratio at `n ≥ 16` (57 configurations) has median 0.951 and range 0.80–1.12; the stored grid gave median 0.955 and range 0.725–1.098 (D3).
* **Mean law.** `AD(t)/AD_∞` follows `φ(t)` to within 0.002 along the whole trajectory.

## Reproduce

```bash
python scripts/regenerate_stage1_data.py                  # fully independent grid, ~5 min
python scripts/regenerate_stage1_data.py --share-decks    # decks reused across m: faster, same marginal distributions
```

Rerunning the script on the platform where the files were produced reproduced `variance_grid_seeded.csv` byte for byte. Floating-point results can differ in the last digits between NumPy versions and platforms.
