# Regenerated extension data

These files were produced for this repository by the code in `experiments/extension/campaign`. They are **new, seeded computations of the quantities that the manuscript reports**. They are not the original runs, which were not recovered (see `../README.md`), and nothing here has been tuned to agree with the manuscript.

| File | What it is | Produced by |
|---|---|---|
| `ising_critical_blocks/L{8,12,16,24,32,48,64,96,128}_start{0,1}.dat` | Raw output of the Wolff sampler at `T_c = 2/ln(1+√2)`: 100 block means of `(|M|, M², M⁴, E, E²)` per run, two runs per size (start 0: random configuration; start 1: ordered configuration). The `.err` files record `L`, the cluster flips between measurements `k`, the mean cluster size, the number of measurements, the seed and the start. | `run_ising_campaign.sh` |
| `ising_critical_sizes.csv` | Per size: relative floors `R_n[|M|]` and `R_n[E]` with jackknife standard errors, susceptibility `χ_M`, Binder cumulant, and the two-start agreement gate (|z| < 3 for both floors). | `analyse_ising_campaign.py` |
| `ising_critical_slopes.json` | Slopes of `log R_M`, `log R_E`, `log χ_M` against `log n` for `L ≥ 16`, with 95% intervals from 2000 block-bootstrap resamples (generator seed 20261115). | `analyse_ising_campaign.py` |
| `curie_weiss_exact.csv` | Exact (not sampled) equilibrium moments of the Curie–Weiss model, `N` from 100 to 50 000, `T` = 0.7, 0.9, 1.0, 1.5, 2.0 with `T_c = 1`. | `curie_weiss_exact.py` |

## Protocol

* **Sampler.** `ising_wolff.c`, single-cluster Wolff updates on the `L × L` torus. A fixed number `k` of cluster flips separates measurements, with `k = ⌈n / ⟨C⟩⌉`, where `⟨C⟩` is the mean cluster size estimated in a warm-up phase and `k` is then held fixed. A schedule that instead flips a fixed number of *spins* between measurements stops at size-biased times and biases the ensemble; the project notes record that mistake (`docs/known-issues.md`, X1), and this campaign does not use that schedule.
* **Sample sizes.** 4 × 10⁵ measurements for `L ≤ 32`; 2 × 10⁵ at `L = 48`; 1.5 × 10⁵ at `L = 64`; 8 × 10⁴ at `L = 96`; 6 × 10⁴ at `L = 128`; the same number for each of the two starts. These are the values in `run_ising_campaign.sh`.
* **Seeds.** `1000 + L` for the random start and `2000 + L` for the ordered start; the generator is xoroshiro128+.
* **Equilibration gate.** The two starts must agree on both floors to within 3 standard errors at every size. All nine sizes pass (`gate_pass` column).
* **Errors.** Jackknife over the 100 blocks of each run; slopes by block bootstrap.

## Result

| `L` | `n` | `R_n[|M|]` | `R_n[E]` | `χ_M` |
|---:|---:|---:|---:|---:|
| 8 | 64 | 0.0707 ± 0.0002 | 0.04153 ± 0.00009 | 2.73 |
| 16 | 256 | 0.0715 ± 0.0003 | 0.01430 ± 0.00004 | 9.31 |
| 32 | 1024 | 0.0723 ± 0.0003 | 0.00454 ± 0.00001 | 31.7 |
| 64 | 4096 | 0.0725 ± 0.0005 | 0.001370 ± 0.000006 | 106.9 |
| 128 | 16384 | 0.0724 ± 0.0007 | 0.000397 ± 0.000003 | 358.7 |

Slopes against `log n` for `L ≥ 16`: `R_M` +0.004 [0.000, 0.008]; `R_E` −0.863 [−0.866, −0.860]; `χ_M` +0.878 [0.875, 0.881]. The exact energy floor from Kaufman's partition function has an effective slope of −0.8637 over the same range (`results/outputs/extension-verification/ising_exact.txt`), and the universal susceptibility exponent is 7/8.

## Comparison with the integrated manuscript (Table 1)

The numbers the repository treats as the reference are those of the integrated manuscript *Bayesian Inference in Partially Observed Critical Systems* (5 November draft, `research/research-extension/integrated-manuscript`), not the earlier extension paper. The raw output behind the manuscript's table was not recovered, so the regenerated runs are an independent re-measurement, not a reproduction. They agree within Monte Carlo error:

| `L` | `R_n[E]` here | `R_n[E]` manuscript | `R_n[|M|]` here | `R_n[|M|]` manuscript | Binder ratio here | Binder ratio manuscript |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 0.01430 ± 0.00004 | 0.0142 | 0.0715 ± 0.0003 | 0.0712 | 0.8580 ± 0.0004 | 0.8586 |
| 32 | 0.004536 ± 0.000011 | 0.004519 | 0.0723 ± 0.0003 | 0.0717 | 0.8563 ± 0.0004 | 0.8571 |
| 64 | 0.001370 ± 0.000006 | 0.001363 | 0.0725 ± 0.0005 | 0.0721 | 0.8558 ± 0.0007 | 0.8562 |
| 128 | 0.000397 ± 0.000003 | 0.000394 | 0.0724 ± 0.0007 | 0.0720 | 0.8555 ± 0.0010 | 0.8565 |

The largest deviation is `R_n[E]` at `L = 16` (2.7 standard errors from the manuscript's 4-digit value; 0.4 from Kaufman's exact 0.014277). All `|M|` differences are positive and below 2.3 standard errors. The slopes agree: `R_M` +0.004 here and +0.002 in the manuscript; `χ_M` exponent 0.878 and 0.877; energy slope −0.863 and −0.864.

Differences in protocol: the manuscript reports `L = 4, 6` (covered here by the exact computations) and `4 × 10⁵` measurements at every size, whereas the regenerated runs use `4 × 10⁵` per start up to `L = 32` and fewer at larger sizes (see Protocol). The manuscript's Curie–Weiss sums run over `10³ ≤ N ≤ 10⁶` at `T = 1.5 T_c` and `0.8 T_c`; `curie_weiss_exact.csv` covers `100 ≤ N ≤ 5 × 10⁴` at five temperatures.

## How this relates to the earlier numbers

The extension paper (2 October) reports an `|M|` floor of about 0.066 and an energy slope of −0.82. The regenerated values, 0.072 and −0.863, agree with the exact small-lattice values and with the integrated manuscript's later runs, and not with those two numbers. See `docs/known-issues.md`, X1 and X2, and R1 below it. This repository does not explain why the original runs differed.

## What is not regenerated

Correlated fields, Algorithm 1 validation, the full `snr_floor` and `detection` campaigns and the P1–P8 campaign are not regenerated. `experiments/extension/verification/field.py` contains independent torus computations of the field and Algorithm 1 results, with stored output.

## Reproduce

```bash
bash experiments/extension/campaign/run_ising_campaign.sh data/extension/regenerated/ising_critical_blocks 2
python experiments/extension/campaign/analyse_ising_campaign.py data/extension/regenerated/ising_critical_blocks data/extension/regenerated/ising_critical
python experiments/extension/campaign/curie_weiss_exact.py data/extension/regenerated/curie_weiss_exact.csv
```

It needs a C compiler and runs for hours on a few cores; `L = 96` and `L = 128` dominate. The sampler's random number generator is self-contained, so rerunning `L = 8` with the same seed reproduced its block file byte for byte (checked when the campaign finished); other platforms and compilers should agree within the stated errors.
