# Phase 4 — Refined Conjecture, Falsifiable Predictions, and Frozen Phase-5 Design

*A standalone scientific deliverable produced before any simulator is built. Its purpose is to make Phase 5 an execution phase, not another round of brainstorming. Epistemic rule enforced throughout: nothing that requires fresh simulation is reported as a result; such items are labelled **[PENDING COMPUTATIONAL VALIDATION]**. This document supersedes the loose "χ_n = o(n)" condition used in `Phase3_interpretation.md`, which Phase 4 corrects (see §1.4).*

---

## Notation (fixed for the remainder of the programme)

- Finite Markov chain with state `x_t`, unique stationary law `ν`, mixing time `t_mix`.
- `n` components; per-component contribution `D_i(x) ∈ [0, B]` bounded.
- Intensive additive observable `O(x) = (1/n) Σ_{i=1}^n D_i(x)`; equilibrium mean `O_∞ = E_ν[O]`, assumed `Θ(1)` unless stated.
- Uniform random `m`-subset `S ⊆ [n]` (without replacement); estimator `Ô_m = (1/m) Σ_{i∈S} D_i`.
- Mixing level (relative displacement) `ϕ = E[O(x_t)] / O_∞ ∈ (0,1]`; `ϕ=1` at equilibrium.
- Within-configuration population variance `V_x = (1/n) Σ_i (D_i(x) − O(x))²`.
- **Sampling signal-to-noise** `ρ²(ϕ) = (E[O])² / E_x[V_x]`, evaluated at mixing level `ϕ`; equilibrium value `ρ² := ρ²(1)`.
- **Observable susceptibility** `χ_n[O] = n · Var_ν(O) = Σ_{i,j} Cov_ν(D_i, D_j) / n` (per-component; `= Σ_j Cov(D_1,D_j)` when translation-invariant).
- **Relative variance floor** `R_n = Var_ν(O)/O_∞² = χ_n[O] / (n · O_∞²)`.
- Certification test: `H0: ϕ = 1` versus `H1: |1−ϕ| ≥ δ`; levels `α` (size), `β` (Type-II), `κ = z_{1−α} + z_{1−β}`.
- Certification sample size `m*(δ)`: smallest `m` meeting `(α, 1−β)` for deviation `δ`.

---

## 1. Refined central conjecture

### 1.1 The exact structural identity (mathematically established — see §2)

For any intensive additive observable estimated from a uniform `m`-subset without replacement, the law of total variance gives, exactly,

```
Var(Ô_m) = (1/m) · (n−m)/(n−1) · E_x[V_x]   +   Var_x(O).
           └──────── sampling term ────────┘     └── floor ──┘
```

Dividing by `O_∞²` and taking `m ≪ n`:

```
CV²(Ô_m)  =  1/(ρ²(ϕ) · m) · (1 + o(1))   +   R_n .          (★)
```

The first term is the RSI mechanism (vanishes as `1/m`); the second is an **`m`-independent floor** equal to the intensive observable's own equilibrium fluctuation. Identity (★) is the spine of the whole conjecture.

### 1.2 The susceptibility exponent

Define `a ∈ [0,1]` by `χ_n[O] ≍ n^{a}` (equivalently `R_n ≍ n^{a−1}` when `O_∞ = Θ(1)`). Three consequences follow directly from (★):

- **Intrinsic precision floor** `δ_floor(n) := √{R_n} ≍ n^{(a−1)/2}` — the best achievable relative precision at `m = ∞`.
- **Certification cost** for `δ > κ·δ_floor(n)`:
  `m*(δ) = κ² / [ ρ² ( δ² − κ² R_n ) ]`; **infeasible** (no `m` suffices) once `δ ≤ κ·δ_floor(n)`.
- **Effective information** `I_eff(m) = 1/CV² = ρ²m / (1 + ρ²m R_n)`, rising linearly (`≈ ρ²m`) then **saturating** at `I_max = 1/R_n ≍ n^{1−a}`, with crossover at `m_× = 1/(ρ² R_n) ≍ n^{1−a}`.

### 1.3 Refined conjecture (the two-condition principle, exponent form)

> **Conjecture C (refined).** Size behaviour of certification is governed by two exponents — the per-component SNR `ρ²` (condition i) and the observable-susceptibility exponent `a` (condition ii):
>
> - **(i) Non-degenerate SNR:** `ρ²(ϕ) → ρ²_∞(ϕ) ∈ (0,∞)`. If instead `ρ² ≍ n^{−2γ}` with `γ>0`, then `m*(δ) ≍ n^{2γ}/δ²` — cost grows with size regardless of correlations (the "inferability exponent" axis).
> - **(ii) Susceptibility regime** (given non-degenerate SNR and a finite-population CLT):
>   - `a = 0` (**summable correlations**): `δ_floor ≍ n^{−1/2}`, `I_max ≍ n`. **Fully scale-free:** `m*(δ) = κ²/(ρ²δ²)`, independent of `n` at every fixed `δ`.
>   - `0 < a < 1` (**sub-extensive but non-summable**): `δ_floor ≍ n^{−(1−a)/2}`, `I_max ≍ n^{1−a}`. **Asymptotically scale-free but pre-asymptotically size-limited:** `m*(δ)` converges to the size-free value only for `n ≫ (κ/δ)^{2/(1−a)}`; below that, a real precision floor blocks small-`δ` certification.
>   - `a = 1` (**extensive correlations**): `δ_floor = Θ(1)`, `I_max = Θ(1)`. **Size-limited:** a non-vanishing tolerance exists below which no sample of any size, at any `n`, can certify.

Condition (i) is RSI Proposition 6.1. Condition (ii) — and especially the sharp distinction between `a=0`, `0<a<1`, and `a=1` — is the genuinely new content. The operative quantity is the **observable's** susceptibility, not a property of the system alone: the same system at the same state can be scale-free in one observable and size-limited in another.

### 1.4 Why this supersedes the Phase-3 statement

`Phase3_interpretation.md` used "`χ_n = o(n)`" as condition (ii). That is **too weak**: 2-D Ising at `T_c` has `χ_n[M] ≍ n^{7/8} = o(n)`, so the crude condition would (wrongly) predict scale-free certification there, yet the pilots show clear degradation. The exponent form fixes this: `a = 7/8 ∈ (0,1)` puts critical magnetization in the *asymptotically-scale-free-but-size-limited* class, with `δ_floor ≍ n^{−1/16}` — a floor that decays so slowly it is effectively present at all reachable sizes. This correction is a Phase-4 result in its own right.

---

## 2. What Phase 3 actually justified

Strict separation. No empirical observation is promoted to a theorem.

**(a) Established mathematically (proven or exactly derivable):**
- The single-tagged-component kernel `P_t(i→j) = 1/n + (δ_{ij} − 1/n)(1−2/n)^t` and hence the exact mean law `E[AD(t)] = (1−(1−2/n)^t) AD_∞`, i.e. `ϕ(t) = 1−(1−2/n)^t`.
- The **product-chain results are exact** (independent `Bernoulli(ϕ/2)` coordinates): `ρ²(ϕ) = ϕ/(2−ϕ)`, `Var(ϕ̂_m) = ϕ(2−ϕ)/m`, `a = 0` exactly.
- The variance-decomposition identity (★) and `Var_ν(O) = χ_n[O]/n` — algebraic identities (law of total variance; definition of susceptibility).
- Equilibrium constants `ρ²=2` (cards, from RSI) and `ρ²=1` (chain, exact).

**(b) Established empirically (measured, robust across the tested range):**
- Cards `ρ²(ϕ) ≈ ϕ/(3/2−ϕ)` across `ϕ∈[0.47,1]` (matches to a few %).
- `n`-independence of the precision law over `n = 256 → 65 536` (chain) and `128 → 2048` (cards).
- Detection `m*` **identical across `n`** for `a=0` systems; `m*δ² → κ²/ρ²` as `δ→0` (cards `≈3.1`, chain `≈6.3`).
- Ising high-`T`: subsample variance `∝ 1/m` with ratio `→1` as `n` grows.
- Ising near-`T_c`: an `m`-independent floor dominates and extra sampling does not reduce it.

**(c) Strongly suggested (consistent, not yet pinned):**
- The cards floor scales as `R_n ≍ 1/n` (the `O(1/n)` correction; ratio drifting `0.8→1.0` with `n`).
- The detection constant is *exactly* `κ²/ρ²` rather than merely proportional.
- The Ising magnetization floor decays slowly with `n` (qualitatively `n^{a−1}`, `a` near `7/8`) — **but the near-`T_c` pilot used finite burn-in and the exponent is not measured.** [PENDING COMPUTATIONAL VALIDATION]

**(d) Conjectured (plausible, unproven):**
- The three-regime exponent principle (Conjecture C), both the `0<a<1` and `a=1` behaviours.
- Observable-dependence of the breakdown (energy vs magnetization at the same `T_c`).
- Information saturation `I_max = 1/R_n` and crossover `m_× = 1/(ρ²R_n)`.
- Self-calibration: `ρ²` and `R_n` estimable from the same subsample, giving `m*` without knowing the system.

**(e) Still unsupported (no evidence yet):**
- Mean-field / extensive-correlation `a=1` hard floor.
- Conserved-dynamics (SSEP) breakdown via slow modes.
- Head-to-head advantage over Gelman–Rubin / ESS.
- Heavy-tailed CLT failure.

---

## 3. Falsifiable predictions

Each turns Conjecture C into something Phase 5 can break. `κ = z_{1−α}+z_{1−β}` (for `α=0.05` one-sided, `β=0.20`: `κ ≈ 2.487`, `κ² ≈ 6.18`).

**P1 — The variance floor equals the whole-system fluctuation.**
- Statement: `lim_{m→∞} Var(Ô_m)/O_∞² = R_n = Var_ν(O)/O_∞²`, independent of `m`.
- Form: a plateau in `Var(Ô_m)` vs `m` at the level `Var_ν(O)`.
- Measure: `Var(Ô_m)` over `m` (many replicates at equilibrium) and, separately, `Var_ν(O)` from full-system samples.
- Support: plateau value within `±10%` of independently-measured `Var_ν(O)` across systems.
- Falsification: the large-`m` variance falls below `Var_ν(O)` (would break the identity ★) or plateaus at a different level not equal to `Var_ν(O)`.

**P2 — Susceptibility exponent controls the floor.**
- Statement: `R_n ≍ n^{a−1}` with `a` the observable-susceptibility exponent (`a=0` summable; `a=7/8` for 2-D Ising magnetization at `T_c`; `a≈0` up to logs for 2-D Ising energy at `T_c`).
- Form: straight line of slope `a−1` on `log R_n` vs `log n`.
- Measure: `Var_ν(O)` vs `n` at equilibrium (cluster-sampled ground truth).
- Support: fitted slope within `±0.1` of `a−1` predicted by the universality class.
- Falsification: magnetization floor slope `≈ −1` at `T_c` (would mean `a=0`, i.e. summable, contradicting criticality), or energy floor slope `≪ 0` (would mean energy is strongly critical).

**P3 — `m*` scaling in the summable regime is size-free.**
- Statement: for `a=0` systems, `m*(δ,n) = κ²/(ρ²δ²)` (`+ O(1/n)`), so `∂ log m*/∂ log n → 0`.
- Form: flat `m*` vs `n`; `m*·δ² → κ²/ρ²`.
- Measure: `m*(δ,n)` by power sweep.
- Support: fitted exponent `b` in `m* ≍ n^b` has `|b| < 0.1` with 95% CI excluding `0.25`; `m*δ²` within `±10%` of `κ²/ρ²`.
- Falsification: `b > 0.25` (CI excluding 0) for a system with `a=0`.

**P4 — `m*` scaling breaks predictably when `a>0`.**
- Statement: at fixed `δ`, `m*(δ,n) = κ²/[ρ²(δ² − κ²R_n)]` with `R_n ≍ n^{a−1}`; hence `m*` **rises** as `n` decreases toward the size where `κ²R_n → δ²`, and certification of `δ ≤ κ√R_n` is infeasible.
- Form: `m*` vs `n` bends upward / diverges as `n` drops; an infeasible region `δ < κ·δ_floor(n)`.
- Measure: `m*(δ,n)` at `T_c` for several `δ`; the largest `δ` still infeasible at each `n`.
- Support: observed infeasibility boundary tracks `δ = κ√R_n` within `±15%`; `m*` matches the floor-corrected formula.
- Falsification: `m*` at `T_c` is flat in `n` (indistinguishable from `a=0`), i.e. no floor effect.

**P5 — System-size independence is a property of `(ρ², a)`, not of the system.**
- Statement: two different systems with the same `(ρ², a=0)` have `m*δ²` agreeing to the same constant `κ²/ρ²`; a single system in two observables with different `a` shows different size-behaviour.
- Measure: `m*δ²` across {chain, cards, Ising high-T}; and Ising-at-`T_c` in magnetization vs energy.
- Support: constants agree across `a=0` systems; the two observables of one critical system fall in different regimes (energy scale-free, magnetization not).
- Falsification: `m*δ²` differs across `a=0` systems beyond the SNR difference, or energy and magnetization behave identically at `T_c`.

**P6 — Correlation range, not exchangeability, is the operative structural condition.**
- Statement: a non-exchangeable system with finite correlation length (`a=0`) is scale-free; an exchangeable-looking but long-range system (`a=1`) is not.
- Measure: Ising high-`T` (non-exchangeable, `a=0`) vs Curie–Weiss below `T_c` (mean-field, `a=1`).
- Support: high-`T` Ising scale-free; Curie–Weiss size-limited with `δ_floor=Θ(1)`.
- Falsification: Curie–Weiss below `T_c` scale-free (would refute that extensive susceptibility breaks it), or Ising high-`T` size-dependent.

**P7 — Ising temperature sweep traces a continuous crossover.**
- Statement: as `T ↓ T_c`, the correlation length `ξ(T)` grows and the floor exponent interpolates; `δ_floor(n,T) ≍ min(1, ξ(T)^{(2−η)}/n)^{1/2}`. Away from `T_c` (`ξ = O(1)`), `a=0`; at `T_c` (`ξ ~ L`), `a=7/8` (magnetization).
- Form: `m*(δ,n)` flat in `n` for `T ≳ 1.4 T_c`, bending up as `T → T_c`.
- Measure: `m*(δ,n)` on a `T`-grid straddling `T_c`.
- Support: monotone onset of size-dependence as `T→T_c`, quantitatively tracking `ξ(T)/L`.
- Falsification: abrupt (not continuous) onset, or size-dependence present even at high `T`.

**P8 — Effective information saturates; sample complexity and information are one curve.**
- Statement: `I_eff(m) = ρ²m/(1+ρ²m R_n)` — linear then saturating at `I_max = 1/R_n`, crossover `m_× = 1/(ρ²R_n)`.
- Form: `I_eff` vs `m` collapses onto the one-parameter curve; plateau `= 1/R_n`.
- Measure: `1/CV²(Ô_m)` vs `m` across `m` and `n`.
- Support: fitted `I_max` matches `1/R_n` (from P1) within `±15%`; `m_×` matches `1/(ρ²R_n)`.
- Falsification: `I_eff` grows without saturation at criticality (no plateau), or plateau `≠ 1/R_n`.

---

## 4. Phase-5 experimental design (frozen)

Designed to distinguish the four competing explanations (E1 "it's a mean", E2 "exchangeability", E3 "bounded SNR", E4 "summable correlations") and to test P1–P8 economically.

### 4.1 Systems and observables

| # | System | Generation of equilibrium ground truth | Off-eq. dynamics | Observable(s) | Predicted `(ρ², a)` |
|---|---|---|---|---|---|
| S1 | Product chain `{0,1}^n` | exact marginal | product Glauber | magnetization | `(1, 0)` |
| S2 | Random transposition (cards) | direct shuffle | random transposition | `AD` | `(2, 0)` |
| S2′| Cards, **fixed-point fraction** | direct shuffle | — | `F=(1/n)Σ1{π(i)=i}` | `(≍1/n, 0)` — SNR-degenerate control |
| S3 | 2-D Ising `L×L`, `T ≥ T_c` | **Wolff cluster** (bypasses critical slowing) | **single-site Glauber** | |M|, energy density `E`, block-mean | `(Θ(1), 0)` off-critical; `M: a=7/8`, `E: a≈0` at `T_c` |
| S4 | Curie–Weiss / mean-field Ising, `N` spins | exact / Metropolis | Glauber | magnetization | `T<T_c: a=1`; `T=T_c: a=1/2` |
| S5 | 1-D Ising, tunable `ξ(T)` | transfer-matrix exact / Glauber | Glauber | magnetization | `a` interpolates `0→` with `ξ/L` |
| S6 | Symmetric simple exclusion (SSEP), `N` sites | uniform (product Bernoulli) | conserved exchange | density in a window | conserved → slow modes; test `a` |

S1–S3 are the **primary** battery; S2′ and S4 are the **discriminating controls**; S5, S6 are **secondary** (crossover and conservation).

### 4.2 Grids

- **Sizes.** Exact systems (S1, S2): `n ∈ {2^7, 2^9, 2^11, 2^13, 2^15}` + two non-power-of-two validators `{3000, 12000}`. Ising (S3): `L ∈ {16, 24, 32, 48, 64, 96}` (`n = L²` up to `9216`; add `L=128` for the critical-observable comparison only). Mean-field (S4): `N ∈ {10^2, 10^3, 10^4, 10^5}`. S5: `L ∈ {256, 1024, 4096}` × `ξ ∈ {2,4,8,16,32,64}`. S6: `N ∈ {512, 4096, 32768}`.
- **Samples.** `m` geometric `{8,16,32,64,128,256,512,1024}` truncated to `m ≤ n/2`, plus the census `m=n` as a consistency check. For power sweeps, refine `m` locally around the bracketed `m*` (bisection to `±5%`).
- **Deviations.** `δ ∈ {0.20, 0.15, 0.10, 0.075, 0.05}`.
- **Temperatures (S3).** `T/T_c ∈ {2.0, 1.4, 1.15, 1.05, 1.00}` (approach from above) + `{0.85}` (broken-symmetry control). `T_c = 2/ln(1+√2) ≈ 2.269`.
- **Mixing levels (off-equilibrium).** `ϕ` targets via `θ = t/t_mix ∈ {0.1,…,2.0}` for S1–S3 to validate `ρ²(ϕ)` along the trajectory.

### 4.3 Replications, seeds, burn-in, autocorrelation, uncertainty

- **Replications.** `R = 2000` for variance/`ρ²` estimates; `R = 5000` for detection-threshold (`m*`) configs; `R = 400` for the expensive Ising `L≥96` (with widened error bars accepted). Rationale: a variance has relative SE `≈ √(2/R)` (≈3% at `R=2000`); a power near `0.8` has SE `≈ √(0.16/R)` (≈0.006 at `R=5000`).
- **Seed protocol.** One documented root seed per experiment; `numpy SeedSequence.spawn` a fresh child per `(system, regime, n, m, δ, replicate)`. Every output row carries its 128-bit seed and a config hash; a run manifest records root seed, code version, library versions. Independent replications must use non-overlapping streams (spawned children guarantee this). A `distribution_level` mode (unseeded) reproduces the RSI protocol for one cross-check.
- **Burn-in / ground truth.** Exact systems: no burn-in (direct generation). MCMC ground truth (S3, S4, S6): generate equilibrium ensembles with a **correct rapid sampler** — Wolff/Swendsen–Wang for Ising, exact product for SSEP density — *not* long local runs, so critical slowing never corrupts the ground truth (this fixes the Phase-3 near-`T_c` caveat). Off-equilibrium detection uses the **physical local dynamics** (single-site Glauber; conserved exchange for SSEP), which is the process whose convergence is being certified.
- **Autocorrelation.** For any locally-generated sample, measure the integrated autocorrelation time `τ_int` of the observable (batch-means or windowed ACF); require replicate spacing `≥ 2τ_int` or independent chains; report effective sample size `R_eff = R/(1+2τ_int/spacing)`. Cluster-generated ensembles are checked for near-independence (`τ_int ≈ 0.5`).
- **Equilibration check.** Every MCMC ground-truth ensemble validated by (a) agreement of the observable from two initial conditions (ordered vs random) to within Monte Carlo error, and (b) stationarity of a running observable trace. Configs failing either are flagged, not silently used.
- **Statistical-uncertainty requirement.** Report all key quantities (`R_n`, `ρ²(ϕ)`, `m*`, floor slope `a`) with bootstrap 95% CIs (≥1000 bootstrap resamples). A prediction "passes" only if its CI lies inside the support band; "fails" only if its CI lies inside the falsification band; otherwise "inconclusive → increase `R`."

### 4.4 Primary vs secondary comparisons

**Primary (must run; decide the project):**
1. `m*(δ,n)` flat vs `n` for S1, S2, S3-high-`T` → tests P3, P5 (E4 vs E1/E2).
2. `m*(δ,n)` vs `n` for S3 at `T_c` (magnetization) → tests P4, P7.
3. **Same system, two observables** at `T_c`: magnetization vs energy → tests P2, P5 (that observable-susceptibility, not the system, is operative).
4. **2×2 factorial** `{ρ² bounded (S2) vs degenerate (S2′)} × {χ summable (S2) vs extensive (S4 below T_c)}` → discriminates all four explanations at once.
5. Floor identity: large-`m` plateau `= Var_ν(O)` across S1–S4 → tests P1.

**Secondary (strengthen; not decisive):**
6. `ρ²(ϕ)` trajectory curves for S1–S3 (confirm the off-equilibrium SNR functions).
7. Information-saturation curve `I_eff(m)` and crossover `m_×` (P8).
8. `ξ`-crossover (S5) and conserved-dynamics (S6).
9. Self-calibration: estimate `ρ̂², R̂_n` from the *same* subsample and predict `m*`; compare to truth.

### 4.5 What counts as confirmation / falsification (global)

- **Confirmation of the refined principle:** Primary-1 flat (`|b|<0.1`), Primary-2 bending per the floor formula, Primary-3 showing the two observables in different regimes, Primary-4 scale-free only in the `(bounded ρ², summable χ)` cell, Primary-5 plateaus at `Var_ν(O)`. All five with CIs in the support bands.
- **Meaningful falsification:** any of — an `a=0` system showing `m*∼n^{b}`, `b>0.25`; a `T_c` magnetization case that is flat in `n` (no floor); energy and magnetization behaving identically at `T_c`; the large-`m` plateau not equal to `Var_ν(O)`; the 2×2 scale-free in a cell the conjecture forbids. Any one of these forces revision of Conjecture C, not merely of a constant.

---

## 5. Control experiments (ruling out spurious "scale-free" effects)

Each named confound gets a dedicated control.

1. **Finite-size artefacts.** Include non-power-of-two sizes (`3000, 12000`) and require the fitted `m*(n)` exponent to be stable when the smallest and largest `n` are dropped (leave-one-size-out). A scale-free verdict that depends on a single size range is rejected.
2. **Inadequate burn-in.** Ground truth via rapid (cluster/exact) samplers, plus the two-initial-condition agreement check (§4.3). A separate *negative* control: deliberately under-burn a chain and confirm the certificate flags it as non-equilibrium (the method should *detect* under-burn-in, not be fooled by it).
3. **Autocorrelation.** Report `R_eff`; recompute every key statistic with replicate spacing doubled and require agreement within CI. If halving the effective sample changes `m*` beyond CI, the uncertainty budget was wrong.
4. **Estimator bias.** Verify `E[Ô_m] = O(x)` empirically at each `(n,m)` (RSI Lemma 2.10) — the sampling distribution must be centred on the full-system value; a nonzero mean-shift signals an implementation or access error. Also confirm `m=n` gives `Ô_m = O` identically.
5. **Hidden full-system information.** Audit that the estimator touches only `m` components: assert the count of `D_i` evaluations equals `m` per replicate; shuffle component labels between generation and observation so no positional leakage can help; confirm results are invariant to a random relabelling of components.
6. **Insufficient repetitions.** Every headline number carries a bootstrap CI; "inconclusive" is a permitted verdict that triggers more `R` rather than a forced call. Power estimates near threshold use `R≥5000`.
7. **Gaussian-approximation failure.** At each `(n,m,ϕ)` record skewness and excess kurtosis of the sampling distribution (RSI reports `<0.19, <0.18` at `m=32`); flag `|skew|>0.3` or `|kurt|>0.5`. In the small-`m` regime (`m ≤ 8`) compare the Gaussian test to an exact/bootstrap calibration of the null; if the size (`α`) is mis-calibrated, report the corrected `m*`. This directly bounds where the CLT-based constant `κ²/ρ²` is trustworthy.

---

## 6. Counterexamples (the ones that could actually falsify Conjecture C)

Ordered by how sharply each targets a specific clause.

**CX1 — Curie–Weiss / mean-field Ising (targets the `a=1` clause).**
Below `T_c` the complete-graph model has extensive magnetization susceptibility, `χ_N[M] = Θ(N)` (`a=1`), across a whole phase, not a point. Prediction: `δ_floor = Θ(1)` — a non-vanishing tolerance below which certification is impossible at any `N`. *Could falsify:* if `m*` is flat in `N` here, condition (ii) is wrong and extensive susceptibility does **not** break scale-freeness — the central refinement collapses. This is the highest-value counterexample because it decouples "extensive correlation" from "criticality."

**CX2 — Sub-extensive-but-problematic: critical magnetization itself (targets the crude `o(n)` condition).**
2-D Ising `M` at `T_c` has `χ_n = o(n)` yet certification degrades. This is already the reason §1.4 refines the condition; Phase 5 must confirm the `0<a<1` behaviour (slow floor, pre-asymptotic size-limitation) rather than either `a=0` or `a=1`. *Could falsify:* if the floor is absent (behaves as `a=0`), the exponent refinement is unnecessary; if it behaves as a hard `Θ(1)` floor (`a=1`), the sub-extensive class is illusory.

**CX3 — Long-range / tunable-`ξ` (targets the crossover clause P7).**
1-D Ising (or a power-law-coupled chain) with `ξ` dialled from `O(1)` to `O(L)`. Prediction: continuous interpolation of size-dependence with `ξ/L`. *Could falsify:* a discontinuous jump, or size-dependence at `ξ=O(1)`, would mean the floor is not governed by correlation range.

**CX4 — Unusual-correlation observables (targets "observable, not system").**
On the *same* critical Ising: energy density (`a≈0`, log-divergent specific heat) vs magnetization (`a=7/8`) vs a conserved/block observable. Prediction: different regimes for different observables of one system. *Could falsify:* if all observables of a system share its fate, the operative quantity is a system property, not `χ_n[O]` — a major simplification/relocation of the theory.

**CX5 — Conserved dynamics with sub-extensive static `χ` but slow modes (SSEP; targets the static-vs-dynamic gap).**
SSEP has a product-Bernoulli stationary law (static `χ` summable, `a=0`) but a conservation law that makes the *dynamics* mix slowly (diffusive hydrodynamic modes). *Could falsify a hidden assumption:* if certification is problematic here despite `a=0`, then a **third** condition is needed — the static susceptibility criterion is not sufficient once the off-equilibrium relaxation carries long-lived modes, and the principle must be split into a *static* (achievable-floor) part and a *dynamic* (time-to-equilibrate) part. This is the counterexample most likely to reveal something the current conjecture misses.

---

## 7. Theorem programme (frozen statements)

**T1 — Anchor theorem (realistically provable; publishable alone).**
*For the random transposition shuffle, `E[AD(t)] = (1−(1−2/n)^t)·AD_∞` exactly, and for a uniform `m`-subset (`1≤m≤n`),*
`CV²(ϕ̂_m) = (3/2−ϕ)/ϕ · (1/m) · (1 + R_n(ϕ,m))`, with `|R_n(ϕ,m)| ≤ C(m/n + 1/n)` *uniformly on `ϕ∈[ϕ_0,1]` for fixed `ϕ_0>0`.*
Hence the certification cost `m*(δ) = κ²/(2δ²)` holds with an explicit, size-uniform error.
- Hypotheses: the shuffle only; the exact kernel; a closed-form bound on `Cov(D_i,D_j)` at equilibrium (computable).
- Status: the mean law and kernel are done; the covariance bound is the remaining work.

**T2 — General sufficiency theorem (ambitious but realistic).**
*Let `{D_i^{(n)}}` be a triangular array of uniformly bounded contributions of an equilibrated chain with `ρ_n²(ϕ)→ρ_∞²(ϕ)∈(0,∞)`, summable susceptibility `χ_n[O]=O(1)`, and a Lindeberg/negligibility condition. Then `CV²(Ô_m) = 1/(ρ_∞²m) + o(1/m) + O(1/n)`, and `m*(δ) → κ²/(ρ_∞²δ²)`, independent of `n`.*
- Hypotheses: uniform boundedness; convergent SNR; **`a=0`** (summable `χ`); Hájek-type finite-population CLT (Lindeberg). Delivers the `a=0` clause of Conjecture C in generality.

**T3 — Necessity / floor theorem (partial; the achievable half is provable).**
*Under identity (★), for any intensive additive observable, `inf_m CV²(Ô_m) ≥ R_n = χ_n[O]/(n O_∞²)`. Consequently, if `χ_n[O] = Θ(n)` (`a=1`), then `inf_m CV(Ô_m) ≥ c > 0` for all `n`, so no size-free certification below tolerance `c` exists.*
- Hypotheses: only the variance identity (already established) plus `χ_n=Θ(n)`. This is the *rigorous* half of "condition (ii) is necessary" — an immediate corollary of (★), hence a **secure** theorem. The full `0<a<1` rate and the two-sided "iff" remain conjectural.

**Meaningful contribution even if T1/T2 stall.**
Identity (★) plus T3 already establish, rigorously, that (a) partial-observation certification of an intensive observable has an information floor equal to the whole-system fluctuation, and (b) extensive susceptibility makes size-free certification impossible. That alone — a clean necessary condition and an exact floor — is a genuine, self-contained theoretical result, independent of whether the sharp sufficiency law is proven. The empirical programme then supplies the `a=0` sufficiency and the `0<a<1` rates as validated phenomenology.

---

## 8. Simulator requirements derived from Phase 4 (handoff basis)

Translating the frozen design into computational requirements (implementation deferred):

1. **`System` interface** exposing: `sample_equilibrium()` (via a *correct rapid* sampler, distinct from the physical dynamics), `step()` (the physical dynamics whose convergence is certified), `true_phi(t)` where analytic, and metadata `O_infinity`, `rho2(phi)`, `susceptibility_exponent_a` (where known). Must support both exact-marginal and MCMC ground truth, tagged per config.
2. **Dual ground-truth path.** Cluster/exact equilibrium generation separated from local-dynamics relaxation, so critical slowing never corrupts ground truth (fixes the Phase-3 caveat). Wolff for Ising; exact product for chain/SSEP-density.
3. **`Observable` interface** with `per_component(x) → D_i[]`, `population_value(x)`, and analytic `E_infinity`, `rho2(phi)`, `a` where available. Multiple observables per system (magnetization, energy, block, fixed-point fraction) selectable independently.
4. **Memory-aware subsampler**, `O(Rm)` not `O(Rn)` (the naïve `argsort(random(R,n))` OOM-ed at `n≈10⁵`); exploit exchangeability (first-`m`) where valid.
5. **Detection/power module** as first-class: sweep `m`, estimate power and Type-I from replicate pools, locate `m*` by bisection to `±5%`.
6. **Diagnostics per config**, emitted automatically: `Var_ν(O)` (floor), `χ_n[O]`, `ρ²(ϕ)`, skewness/kurtosis of the sampling law, `τ_int`, `R_eff`, and the two-initial-condition equilibration check.
7. **Self-calibration hooks**: estimate `ρ̂²` and `R̂_n` from the observed subsample and predict `m̂*`, stored alongside the ground-truth `m*`.
8. **Seed/reproducibility layer**: `SeedSequence.spawn` per `(config, replicate)`; per-row seed + config hash; run manifest; optional unseeded `distribution_level` mode.
9. **Access audit**: assert exactly `m` component reads per replicate; component-relabelling invariance test (control #5).
10. **Storage + bootstrap-CI scorer + figure factory**: tidy per-`(config,replicate)` records; bootstrap CIs on every headline quantity; automatic support/falsification/inconclusive verdicts against §3–§4.

## 9. Open questions (smallest set that would materially change the project)

1. **Is the operative quantity the observable-susceptibility exponent `a`, confirmed by the energy-vs-magnetization contrast at one `T_c`?** If yes, the theory is about observables; if no, it relocates to a system property. (Primary-3.)
2. **Does an extensive-`χ` system (`a=1`, mean-field below `T_c`) exhibit the predicted hard, non-vanishing precision floor?** This is the single cleanest test of condition (ii)'s necessity. (CX1.)
3. **Do conserved dynamics (SSEP) break certification despite summable static `χ`?** If so, a distinct *dynamic* condition must join the static one, and the principle splits in two. (CX5.)
4. **Is the method self-calibrating** — can `ρ²` and `R_n` be recovered from the same subsample to predict `m*` without prior system knowledge? This decides whether it is a deployable diagnostic or only an analyzed phenomenon.

## 10. Final research claim (end of Phase 4)

> Partial observation of a large equilibrated stochastic system, through a uniform random sample of `m` of its `n` components, certifies the equilibrium **mean of a chosen intensive additive observable** — a necessary, not sufficient, condition for full equilibrium. The relative error of this certificate decomposes exactly (law of total variance) into a sampling term `1/(ρ²m)` and an `m`-independent floor `R_n = χ_n[O]/(nO_∞²)` equal to the whole-system fluctuation of the observable. Consequently the sample size to certify agreement to relative tolerance `δ` is `m*(δ)=κ²/[ρ²(δ²−κ²R_n)]`, and the number of observations is **independent of system size precisely when the observable's susceptibility is bounded (`a=0`)**: the per-component signal-to-noise ratio `ρ²` sets the size-free constant, while an extensive susceptibility (`a=1`, e.g. mean-field) imposes an irreducible precision floor that no sample of any size overcomes, and a sub-extensive-but-non-summable susceptibility (`0<a<1`, e.g. 2-D critical magnetization) makes certification only asymptotically size-free with a slowly vanishing floor. *Established now:* the exact variance identity and floor lower bound (T3), the exact shuffle and product-chain laws, and the size-independence of `m*` in `a=0` systems across up to four orders of magnitude (empirical). *Pending computational validation:* the `a`-dependent floor exponents, the extensive-`χ` hard floor, the static-vs-dynamic distinction under conservation, and self-calibration.
