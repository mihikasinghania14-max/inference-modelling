# Mathematical background

Notation follows the two papers. Statements marked *(verified)* are checked by a test or script in this repository; the rest are quoted from the papers.

## 1. The shuffle and the observable

A deck of `n` cards starts ordered; at each step two positions `a_t, b_t` are drawn independently and uniformly from `[n]` and the cards there are swapped (identity if `a_t = b_t`). The chain is irreducible, aperiodic and symmetric, so the uniform law `ν` on `S_n` is stationary. Mixing is sharp near `t_mix = ½ n ln n` (Diaconis & Shahshahani 1981); the papers use `⌊½ n ln n⌉`.

`D_i(π) = |π(i) − i|`, `AD(π) = (1/n) Σ_i D_i`, and `F = n·AD` is Spearman's footrule.

At equilibrium `π(i)` is uniform, so `E|π(i) − i| = [(i−1)i + (n−i)(n−i+1)] / (2n)` (smallest at the centre, largest at the ends) and

```
AD_∞ = (n² − 1)/(3n)                                   (verified by enumeration, n ≤ 7)
E[D²] = (n² − 1)/6,   Var(D) = n²/18 · (1 + O(n⁻²))    (single card)
Var(F) = (n + 1)(2n² + 7)/45                           (verified by enumeration, n ≤ 8; Diaconis & Graham 1977)
```

## 2. Partial observation

`ÂD_m = (1/m) Σ_{i∈S} D_i` with `S` a uniform `m`-subset. Because `P(i ∈ S) = m/n`, `E_S[ÂD_m | π] = AD(π)`, hence `E[ÂD_m] = AD_∞` for every `m` (RSI Lemma 2.10). A fixed window of positions is biased because of the position dependence of `E D_i`.

**Exact variance.** Conditioning on `π` and using the finite-population formula `Var_S(ÂD_m | π) = V_π/m · (n−m)/(n−1)` with `V_π = (1/n) Σ (D_i − AD)²` gives (RSI Eq. E.3; extension Lemma 1)

```
Var(ÂD_m) = (n − m)/(m(n − 1)) · E_ν[V_π]  +  Var_ν(AD)
E_ν[V_π]  = (n − 2)(n + 1)(5n² + n + 12) / (90 n²)       (verified by enumeration, n ≤ 6)
Var_ν(AD) = Var(F)/n²
```

Both closed forms are in `rsi_shuffle.theory`. For `m ≪ n` the first term is `n²/(18m)(1 + O(m/n))` (RSI Theorem 5.1); the floor `Var_ν(AD) = O(n)` is lower order in that regime. For large `n` the ratio of the exact variance to `n²/(18m)` tends to `1 − m/(5n)` (tested), which accounts for the departures of the observed ratios from 1 in RSI Figure S3.

## 3. Scale invariance (Stage 1)

Per-card signal-to-noise ratio `ρ² = (E D)²/Var(D) → 2`. The Fisher information about `μ = E[AD]` in the Gaussian model `ÂD_m | μ ~ N(μ, n²/18m)` is `18m/n²`; with `φ = μ/AD_∞`, `I(φ) = AD_∞² · 18m/n² = 2m(1 + O(n⁻²))`, so `CV(ÂD_m) = 1/√(2m)`, independent of `n`. For a two-sided credible interval with relative half-width `ε` at level `1 − α`: `m*(ε) = z²_{1−α/2}/(2ε²)`, giving 48, 192 and 768 at `ε = 0.20, 0.10, 0.05` (verified, `tests/test_rsi_shuffle_simulation.py`).

**Posterior.** Prior `μ ~ N(μ₀, σ₀²)`, likelihood as above: precision `1/σ_post² = 1/σ₀² + 1/σ_m²`, mean `μ_post = σ_post²(μ₀/σ₀² + ÂD_m/σ_m²)` (`rsi_shuffle.bayes`). RSI Figure 8: `n = 512`, `m = 128`, `ÂD_m = 168.2`, `σ₀ = 100` gives `N(168.23, 10.61²)` (verified).

**General principle (RSI Proposition 6.1).** If `n` exchangeable components carry bounded contributions obeying a finite-population CLT and the per-component SNR `ρ` stays bounded away from 0 and ∞, then `CV = ρ⁻¹ m^{−1/2}(1 + o(1))` and `m* = z²/(ρ²ε²)`, independent of `n`.

## 4. The floor (Stage 2)

For any intensive additive observable `O = (1/n) Σ D_i` and uniform `m`-subset: `Var(Ô_m) = (n−m)/(m(n−1)) E_ν[V_x] + Var_ν(O)` (extension Lemma 1, the same law-of-total-variance argument as in section 2). Dividing by `O_∞²`:

```
CV²(Ô_m) ≥ R_n = Var_ν(O)/O_∞²      for every m                           (Corollary 1)
χ_n = n Var_ν(O)
χ_n = O(1)  ⇒  R_n = Θ(1/n) → 0      (size-free certification possible)
χ_n = Θ(n)  ⇒  R_n = Θ(1)            (no m certifies a tolerance δ < κ√R_n,  κ = z_{1−α} + z_{1−β})
```

**Shuffle.** `R_n = (2n² + 7)/(5(n−1)²(n+1))`, `n R_n → 2/5` (verified, enumeration `n ≤ 7` and symbolic derivation). Certification cost including the floor: `m*(δ) = κ²/[ρ²(δ² − κ²R_n)]`, diverging as `δ ↓ κ√R_n` (in a finite system `m*` is capped at `n`); for the shuffle `ρ² = 2` and `m*δ² → κ²/2 ≈ 3.09` with `κ = 2.486` at `α = 0.05`, `β = 0.20` (verified). The mean law of `AD` from the identity under the independent-pair chain is `E[AD(π_t)] = AD_∞ (1 − (1 − 2/n)^t)` (verified by exact transition-matrix calculation at `n = 4`), so the mixing level is `φ(t) = 1 − (1 − 2/n)^t`, and `ρ²(φ) = φ/(3/2 − φ)` to leading order (extension working notes; not independently verified in this repository).

**Size-free sufficiency (extension Proposition 1).** Under a triangular-array CLT (Hájek 1960), bounded contributions, `ρ² → ρ²_∞ ∈ (0,∞)` and summable susceptibility `χ_n = O(1)`: `m*(δ) → κ²/(ρ²_∞ δ²)`.

## 5. Observable dependence at criticality

`χ_n ≍ n^a`. For the 2D Ising model at `T_c = 2/ln(1+√2)` on an `L × L` torus (`n = L²`): `⟨|M|⟩ ~ L^{−β/ν}` and `Var(|M|) ~ L^{−2β/ν}`, so `R[|M|]` tends to a constant and the susceptibility exponent is `a = γ/(2ν) = 7/8` (with `γ = 7/4`, `ν = 1`, `β = 1/8`). The energy density has a vanishing floor: its exact torus floor `R_n[E]` (Kaufman 1949; computed in `experiments/extension/verification/ising_exact.py`) has effective log–log slope −0.864 over `16 ≤ L ≤ 128` and carries logarithmic corrections. Curie–Weiss: `a = 0` above `T_c` (`Var(M) ∝ 1/N`) and `a = 1` below (`Var(M) → m₀²`); the limits `√2 − 1` and `π/2 − 1` at `T_c` and `1.5 T_c` are checked in `cw.py`.

## 6. Correlated fields and Algorithm 1

A stationary Gaussian field on an `L × L` torus with correlation length `ξ` has relative floor `σ²/(μ² n Σ_x G(x)²)` (checked numerically in `field.py`). Algorithm 1 estimates `Var_ν(O)` from a single realisation by batch means: partition into blocks of size `b`, take the plateau of `b·Var(block means)` as `b → ∞`, inflate to a one-sided upper bound `R̂⁺_n`, refuse to certify when `δ ≤ κ√R̂⁺_n`. It is a *proposed* method; hypotheses H1–H3 in the extension state what would validate it.

## References

Diaconis & Graham (1977) *J. R. Stat. Soc. B* 39(2); Diaconis & Shahshahani (1981) *Z. Wahrsch. verw. Gebiete* 57(2); Hájek (1960); Kaufman (1949) *Phys. Rev.* 76; Onsager (1944) *Phys. Rev.* 65; Wolff (1989) *Phys. Rev. Lett.* 62; Cochran (1977) *Sampling Techniques*; Levin, Peres & Wilmer (2009) *Markov Chains and Mixing Times*. Full references are in the two papers.
