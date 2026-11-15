# Independent verification scripts for the extension

These scripts were written during the preparation of the integrated manuscript (October 2026) to check the extension's algebra and numerics **independently of the original runs**. They are checks, not the experiments behind the extension paper's figures; the latter were not found (see `data/extension/README.md`).

| Script | What it checks | Method | Stored output |
|---|---|---|---|
| `shuffle_algebra.py` | `Σc_i`, `Σc_i²`, `Var F = (n+1)(2n²+7)/45`, `R_n`, `E[V_x]`, the `Δ_n` identity and the limits `ρ² → 2`, `χ → 2/45`, `n R_n → 2/5`; `E[F]` and `Var F` by enumeration of `S_n` for `n = 2…8`; finite-`n` `m*(ε)`; `κ`; predicted (analytic) coverage of the RSI posterior at `n = 100` and the effect of the prior width `σ₀` on the posterior sd | SymPy (exact), enumeration, closed forms | `shuffle_algebra.txt` |
| `cw.py` | Curie–Weiss limits `√2 − 1` and `π/2 − 1` at `T_c` and `1.5 T_c`, and the exponent `a` at `1.5 T_c`, `T_c`, `0.8 T_c` | exact sums up to `N = 10⁷` | `cw.txt` |
| `ising_exact.py` | 2D Ising at `T_c`: `R[e]`, `R[\|M\|]` and the Binder ratio by enumeration (`L = 4`) and transfer matrix (`L = 6, 8`); exact energy floor from Kaufman's torus partition function for `L = 4 … 128` and its log–log slope over `16 ≤ L ≤ 128` | enumeration, transfer matrix, mpmath | `ising_exact.txt` |
| `wolff.c`, `wolff_analyze.py` | `R[\|M\|]` and the Binder ratio on `16²` and `32²` tori from a Wolff single-cluster sampler (a fixed number `k` of cluster flips between measurements, `k = ⌈n/⟨C⟩⌉`) | Monte Carlo, 100 block means, jackknife errors | `wolff_L16.dat/.err`, `wolff_L32.dat/.err`, `wolff_summary.txt` |
| `shuffle_sim.py` | Mean AD, variance of AD and fixed points at `t_mix` and `2 t_mix` for `n = 512`; total variation distance between Poisson(2) and Poisson(1); coverage of the RSI posterior at `n = 100` | Monte Carlo, seed 12345 | `shuffle_sim.txt` |
| `field.py` | Gaussian random field on a torus: closed-form floor `σ²/(μ² n ΣG²)` against simulation; coverage collapse of the sample-mean interval as the correlation length grows; Algorithm 1 (batch-means floor estimator) on single fields | Monte Carlo, seed 2026 | `field.txt` |

## Running

```bash
bash experiments/extension/verification/run_all.sh            # runs everything, compares with the stored outputs
bash experiments/extension/verification/run_all.sh /tmp/out    # keep the outputs in /tmp/out
```

`run_all.sh` requires Python with NumPy, SciPy, SymPy and mpmath and a C compiler. It reports `OK` for each script whose output is identical to the stored output. When this repository was assembled, all seven comparisons were `OK` (the Wolff runs use seed `1` and 200 000 measurements: `./wolff L 200000 1`). Run time: about 25 seconds.

## Results at a glance (stored in `results/outputs/extension-verification`)

* `Var F = (n+1)(2n²+7)/45` and `R_n = (2n²+7)/(5(n−1)²(n+1))` agree with enumeration for `n ≤ 8`.
* `R[|M|]` on small tori: 0.06917 (`L = 4`), 0.06995 (`L = 6`), 0.07061 (`L = 8`); Wolff: 0.0703 ± 0.0004 (`L = 16`), 0.0717 ± 0.0005 (`L = 32`).
* Exact effective slope of `R_n[E]` over `16 ≤ L ≤ 128`: −0.8637.
* At `n = 512`, `t = t_mix = 1597`: mean AD 170.198 (equilibrium 170.666), mean fixed points 2.016 (equilibrium 1); at `2 t_mix` 170.694 and 0.988.

How these compare with the numbers printed in the extension paper is discussed in `docs/known-issues.md` (X1–X3, D6).
