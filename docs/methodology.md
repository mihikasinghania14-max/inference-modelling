# Methodology

## Stage 1 (RSI)

1. **System.** Random transposition shuffle on `n` cards, `n = 2…4096` (powers of two) plus validation sizes 10, 50, …, 3000 to rule out powers-of-two artefacts. Independent uniform pair draws (self-loops allowed).
2. **Calibration.** Mean `AD` after `t_mix = round(½ n ln n)` steps against `AD_∞ = (n²−1)/(3n)`, `R = 1000` replicates per size; error bars are ±1 sd. Trajectories extend to `2 t_mix`; normalising by `AD_∞` and `t_mix` tests collapse across sizes.
3. **Independent check.** Final position of a tagged card after mixing, tested for uniformity (a statistic independent of `AD`).
4. **Partial observation.** After mixing, reveal a uniform `m`-subset without replacement and form `ÂD_m`. A new subset is drawn in every replicate. The grid has 2 to 8 values of `m` per `n`, up to the census `m = n`, which isolates permutation randomness from sampling randomness.
5. **Variance law.** For every `(n, m)` the empirical variance of `ÂD_m` over replicates is compared with `n²/(18m)`; skewness and excess kurtosis assess normality (Hájek CLT).
6. **Inference.** Gaussian likelihood `ÂD_m | μ ~ N(μ, n²/18m)`, conjugate prior `N(AD_∞, σ₀²)` with `σ₀ = 100` at `n = 512` and a sensitivity check over `σ₀ ∈ {50, 75, 100, 150, 200}`; closed-form posterior and credible interval.
7. **Sample-size law.** `m*(ε)` read off the posterior criterion for `ε = 0.05, 0.10, 0.15, 0.20` at 95% across sizes; constancy in `n` compared with constant, power-law, logarithmic and linear fits by AIC and BIC.

## Stage 2 (extension)

1. **Decomposition.** Exact law of total variance separates a `1/m` sampling term from an `m`-independent floor `Var_ν(O)`; the relative floor `R_n` and susceptibility `χ_n = n Var_ν(O)` classify systems.
2. **Certification problem.** One-sided test that `φ = E[O]/O_∞` has reached 1 to relative tolerance `δ` at size `α` and power `1 − β`; cost `m*(δ) = κ²/[ρ²(δ² − κ²R_n)]`.
3. **Systems.** (S1) a product chain with exact marginals; (S2) the shuffle, with observables `AD` and the fixed-point fraction (a degenerate-SNR control); (S3) the 2D Ising model, equilibrium ensembles from a Wolff cluster sampler, relaxation with single-site (Glauber) dynamics; (S4) Curie–Weiss; correlated Gaussian fields on a torus.
4. **Measurement discipline (as implemented in `src/rsi_partial`).**
   * *Import boundary.* All theoretical predictions live in `analysis/predictions.py`; a test asserts no measuring module imports it.
   * *Veil.* Inference code receives only the `m` revealed values; an access audit asserts exactly `m` reads.
   * *Seed protocol.* One root seed; each `(configuration, replicate)` gets a child seed from `numpy.random.SeedSequence` keyed by a hash of the configuration, so any single configuration can be re-run independently.
   * *Equilibration gate.* Locally generated (MCMC) ensembles are accepted only if two independent initial conditions agree.
   * *Confidence intervals.* Every fitted slope carries a bootstrap 95% interval.
   * *Pre-registered bands.* Eight predictions P1–P8 carry support/falsify bands fixed before the measurements (`analysis/predictions.py`, `analysis/score.py`).
   * *Storage.* Parquet records plus a JSON manifest (seed, library versions, configuration, timestamp).
5. **Single-realisation floor estimation (proposed).** Algorithm 1 (batch means); validation protocol on the shuffle, Curie–Weiss and Ising with acceptance criteria H1–H3.

## What was *not* done

The full P1–P8 campaign at full replicate counts and the validation of Algorithm 1 are not available in this repository (`data/extension/README.md`).
