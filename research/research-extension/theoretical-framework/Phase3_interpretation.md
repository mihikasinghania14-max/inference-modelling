# Phase 3 → Rigorous Claim: A Mathematical Interpretation of the Three-System Pilots

*No code. This document converts the Phase-3 discoveries into the strongest claim that is actually defensible, scrutinises the working slogan, and separates what is proven, what is a theoretical consequence, and what remains conjecture. It is written to be read alongside the RSI paper and `Phase3_findings.md`.*

---

## 0. The object under the slogan: what is actually certified

The working phrase is **"scale-free partial-observation equilibrium certificate."** Three of its four words are load-bearing and one is an overclaim. Take them apart.

The procedure observes one realisation `π_t` of a Markov chain, reveals a uniform random `m`-subset `S` of its `n` components, forms the additive estimator
`Ô_m = (1/m) Σ_{i∈S} D_i(π_t)`
of the intensive observable `O(π_t) = (1/n) Σ_i D_i(π_t)` (for the shuffle, `D_i = |π_t(i)−i|` and `O = AD`), and compares `Ô_m` to the **known** equilibrium mean `O_∞`. From the sampling distribution it either reports an interval for the current value or decides `H0: E[O]=O_∞` versus `H1: |E[O]−O_∞| ≥ δ·O_∞`.

Now test the slogan against each candidate meaning of "certify":

| Candidate | Verdict | Why |
|---|---|---|
| (1) equilibrium itself (`law(π_t)=ν`) | **No** — overclaim | Matching one scalar functional is a projection. Configurations with `O=O_∞` but `‖law(π_t)−ν‖_TV = Ω(1)` exist; the RSI paper says as much in §6. |
| (2) closeness of the observable's *distribution* to equilibrium | **No** | Only the **first moment** is tested. The distribution of `O` (variance, tails) is never matched. |
| (3) a confidence **bound on the observable's mean deviation** | **Yes** | The interval `ϕ̂ ± z·CV` is exactly a `(1−α)` probabilistic bound on `ϕ = E[O]/O_∞`. |
| (4) **detectability of a prescribed deviation** `δ` in the observable's mean | **Yes** | The power statement: any relative shift `≥ δ` of `E[O]` is caught with controlled Type-I/Type-II error at `m ≥ m*(δ)`. |
| (5) something else | **This is the precise object** | A **fixed-sample hypothesis test / interval estimate for a single intensive equilibrium functional (the mean of a chosen additive observable), with sample complexity independent of system size in a specified correlation regime.** Options (3) and (4) are its estimation and decision faces. |

So the honest reading is: **not an "equilibrium certificate" but an equilibrium-*mean* monitor for one observable** — a *necessary-condition screen*. Passing it is necessary for mixing, not sufficient. The alternative hypothesis is not "any non-equilibrium law" but the **one-parameter relaxation family** `ϕ = 1−δ` traced by the dynamics; against adversarial alternatives that preserve `E[O]` the test has power equal to its size (i.e., none), by construction.

**Strongest defensible claim (the sentence the whole project should be able to stand behind):**

> For an intensive additive observable `O = (1/n)Σ_i D_i` of a finite Markov chain with **known** stationary mean `O_∞`, a uniform random sample of `m` components gives an estimator whose relative error — and hence the sample size `m*(δ)` needed to test the observable's mean against a relative deviation `δ` at fixed error levels — is **independent of the system size `n`**, provided (i) the per-component signal-to-noise ratio `ρ² = E[D]²/Var(D)` is bounded away from `0` and `∞`, and (ii) the summed per-component covariance (susceptibility) is **sub-extensive**, `χ_n = Σ_j Cov(D_1,D_j) = o(n)`. The test certifies agreement of *the observable's mean* with its equilibrium value — a necessary, not sufficient, condition for equilibrium. When (ii) fails (e.g. at a critical point) the estimator variance acquires an `m`-independent floor set by the susceptibility, and size-independence is lost.

Everything below argues that this claim, and only this claim, is what the pilots earn.

---

## 1. The central new research question

*Under partial observation of a large equilibrated stochastic system, is the sample complexity of certifying — to a prescribed relative tolerance `δ` and confidence `1−α` — that an intensive additive observable has reached its equilibrium mean **independent of the system size**, and **which structural property of the components' correlation structure** decides whether it is?*

This is genuinely distinct from the RSI paper, which asked only for the *value* of `AD_∞` at equilibrium and proved an `n`-independent estimation `CV`. The new question is about **testing** (has it arrived?) and about the **structural condition** governing when size-independence holds.

## 2. The central hypothesis

A **two-condition principle**. Size-independent certification holds *iff* both

- **(i) bounded per-component SNR:** `ρ² = E[D]²/Var(D) = Θ(1)` as `n→∞` (RSI Prop. 6.1's condition), and
- **(ii) sub-extensive susceptibility / summable correlations:** `χ_n = Σ_j Cov(D_1,D_j) = o(n)`,

in which regime the cost is `m*(δ) = (z_{1−α}+z_{power})² / (ρ² δ²)`, with `ρ²` the single system-specific constant. When (ii) fails, the estimator variance picks up an irreducible floor `≈ χ_n/n` that subsampling cannot reduce, and `m*` becomes `n`-dependent.

The off-equilibrium refinement (from the exact single-component kernel): along the relaxation the constant `ρ²` becomes a **function** `ρ²(ϕ)`, equal to `ϕ/(3/2−ϕ)` for the shuffle and `ϕ/(2−ϕ)` for the product chain, recovering the equilibrium constant at `ϕ=1`.

## 3. The contribution beyond the RSI paper — stated exactly

RSI establishes: at equilibrium, estimation of `AD_∞` has `CV=1/√(2m)`, `n`-independent; Prop. 6.1 gives a **sufficient** condition (**exchangeable** components + bounded `ρ²`) for size-independent **estimation**.

The extension adds four things RSI does not contain:

1. **Estimation → certification.** Reframes the problem as a hypothesis test on the equilibrium functional — which is what the paper's own opening question ("has the system equilibrated?") actually demands, and which it never solved. Yields `m*(δ)=(z_α+z_β)²/(ρ²δ²)`.
2. **A trajectory, not a point.** From the exact single-card kernel `P_t(i→j)=1/n+(δ_{ij}−1/n)(1−2/n)^t`, the mean law `E[AD(t)]=(1−(1−2/n)^t)AD_∞` (which *derives* the paper's empirical Fig. 4 residual) and the SNR **function** `ρ²(ϕ)`. RSI has only the `ϕ=1` value.
3. **The hidden second condition.** By leaving the exchangeable world (Ising), the pilots expose that bounded SNR is **not sufficient**: one also needs summable correlations. RSI's exchangeability *tacitly guaranteed* (ii); the extension makes it explicit, separable, and shows it is the binding one in the interesting failure.
4. **A physically meaningful, universal failure boundary.** The breakdown sits at criticality and is governed by susceptibility (hence by critical exponents), not by system idiosyncrasy — turning "a limitation" into a predictive statement about *which* systems the method serves.

## 4. Status of each Phase-3 result — empirical / theoretical / conjectural

**Exact theoretical consequences (provable now, some already derived):**
- Single-component kernel and the mean law `ϕ(t)=1−(1−2/n)^t`; corollary `ϕ(θ)≈1−n^{−θ}`. *(Rigorous — complete-graph spectral fact + linearity.)*
- Product-chain results are **exact**: independent `Bernoulli(ϕ/2)` coordinates give `ρ²(ϕ)=ϕ/(2−ϕ)` and `Var(ϕ̂)=ϕ(2−ϕ)/m` with no correction. *(Exact, not asymptotic.)*
- The variance **decomposition** `Var(Ô_m)=E_σ[V_σ]/m·FPC + Var_σ(O)` is an identity (law of total variance); the second term `= χ_n/n` in the intensive scale. *(Rigorous structural identity — this is the mathematical heart of §6.)*

**Theoretical consequences at leading order (derivable, not yet with controlled error):**
- Cards `ρ²(ϕ)=ϕ/(3/2−ϕ)` and `Var(ϕ̂)=ϕ(3−2ϕ)/(2m)`. *(Leading order in `n`; uses equilibrium image law for the off-equilibrium variance — exact for the mean, heuristic for the covariance term.)*
- Detection cost `m*=(z_α+z_β)²/(ρ²δ²)`. *(Follows from Gaussian/finite-population-CLT approximation; finite-`m` corrections unquantified.)*

**Empirical observations (measured, not yet derived):**
- Numerical match of `ρ²(ϕ)` and the precision law across `ϕ` and across `n` (2.4 orders of magnitude).
- `m*` identical across `n`; `m*δ² → (z_α+z_β)²/ρ²` as `δ→0`.
- The Ising floor: `m`-independent, dominates near `T_c`, and decays slowly with `n`.

**Still conjectural:**
- The `O(1/n)` **order** of the cards covariance correction (measured ≈0.8→1.0; not bounded).
- The **general two-condition principle**, in either direction; the **necessity** ("iff") half especially.
- The **critical-exponent scaling** of the Ising floor (`≈ χ/n ~ L^{2−η}/L² = L^{−1/4}` in 2D) — consistent with 3 lattice sizes and *approximate* (finite burn-in near `T_c`), not established.
- That the Gaussian constant is *exactly* `(z_α+z_β)²/ρ²` rather than `(1+o(1))` times it.

## 5. What the cards → product chain → Ising progression demonstrates

It is a **controlled sweep of correlation structure with everything else held fixed** (intensive additive observable, bounded per-component SNR). The dependence between components runs:

`product chain` (independent, `χ_n=0`) → `cards` (exchangeable, weak, `χ_n=O(1)` with a small negative sum-constraint coupling) → `Ising high-T` (short-range positive, finite `ξ`, `χ_n=Θ(1)`) → `Ising critical` (long-range, `ξ~L`, `χ_n→∞`).

The demonstration: **bounded SNR is not the operative condition** — all four systems have `ρ²=Θ(1)`, yet only the first three are scale-free. The variable that flips the outcome is **the summability of correlations**. The chain and cards alone could never reveal this (both satisfy (ii) trivially); Ising is *required*, and criticality is *required*, to make (ii) bind. This is why the progression is a designed experiment, not an assortment of systems.

## 6. Why the Ising critical failure is scientifically important, not a mere limitation

Four reasons, in increasing depth:

1. **It is a positive identification, not an absence.** It converts "the method works" into "the method works *exactly when* `χ_n=o(n)`, and here is the mechanism when it doesn't: the estimator variance is floored by the whole-system fluctuation `Var_σ(O)=χ_n/n`, which subsampling cannot touch." A boundary with a mechanism is a result.
2. **The failure is universal, not idiosyncratic.** Its size-dependence is set by critical exponents (`γ, η, ν`), so the *rate* at which `m*` grows should be predictable from the universality class alone — a far stronger and more transferable statement than any single system's success.
3. **It coincides with the hardest real problem.** Critical slowing down is precisely where practitioners' MCMC convergence diagnostics struggle. Showing that partial-observation certification fails *there*, for the *same structural reason* (diverging susceptibility / non-summable correlations), explains why no cheap subsample can rescue convergence checking at a phase transition — an explanatory, not merely cautionary, contribution.
4. **It disciplines the claim.** A method that never fails in any tested regime is usually one whose failure regime hasn't been found. Finding it *by design* is evidence of rigor and is what separates a defensible claim from an overclaim.

## 7–8. Can the results support a general principle? — and its precise (conjectural) form

Yes. The law of total variance gives the exact decomposition for any intensive additive observable estimated from a uniform `m`-subset without replacement:

```
Var(Ô_m) = (1/m)·(n−m)/(n−1)·E_σ[V_σ]   +   Var_σ(O),      V_σ = (1/n)Σ_i (D_i − O)².
                └── sampling term, ∝ 1/m ──┘        └─ floor, m-independent ─┘
```

Dividing by `O_∞²` and taking `m ≪ n`:

```
CV²(Ô_m)  =  1/(ρ² m)·(1+o(1))   +   χ_n/(n·O_∞²),        χ_n := Σ_j Cov(D_1, D_j).
```

The first term is the RSI mechanism; the second is the **whole-system fluctuation**, equal to the susceptibility per component over `n`. This yields the principle:

> **Conjecture (Scale-free certification principle).**
> Let `{D_i^{(n)}}_{i=1}^n` be a triangular array of bounded per-component contributions of an equilibrated chain, with intensive mean `O_∞^{(n)} = Θ(1)`, per-component SNR `ρ_n² = E[D]²/Var(D)`, and susceptibility `χ_n = Σ_j Cov(D_1,D_j)`. Estimate `O_∞` from a uniform `m`-subset. Then the sample size `m*(δ,α,β)` to test `H0: E[O]=O_∞` against `H1:|E[O]−O_∞|≥δO_∞` at levels `α,β` satisfies
>
> `m*(δ) → (z_{1−α}+z_{1−β})² / (ρ_∞² δ²)`,  independent of `n`,
>
> **if and only if** (i) `ρ_n² → ρ_∞² ∈ (0,∞)` **and** (ii) `χ_n = o(n)` (equivalently `Var_σ(O)→0`, a law of large numbers across components), together with a Lindeberg/uniform-integrability condition ensuring the finite-population CLT. If (ii) fails with `χ_n = Θ(n^a)`, `0<a≤1`, then `CV²` has an irreducible floor `Θ(n^{a−1})` and `m*` is `n`-dependent (diverging as `a→1`).

Sufficiency of (i)+(ii) is plausibly provable (below). Necessity — the "iff" — is the ambitious part and should be advertised as conjecture; (ii)'s necessity is what the Ising data support but do not prove.

## 9. The strongest theorem we could realistically prove

Two tiers; aim for both, promise only the first.

**Anchor theorem (realistically provable, publishable on its own).** *For the random transposition shuffle, `E[AD(t)] = (1−(1−2/n)^t)·AD_∞` exactly, and for a uniform `m`-subset with `m≤n`,*
`CV²(ϕ̂_m) = (3/2−ϕ)/(ϕ) · (1/m) · (1 + R_n(ϕ,m))`, with `|R_n| ≤ C(m/n + 1/n)` *uniformly on `ϕ∈[ϕ_0,1]`,*
*so the certification cost `m*(δ)=(z_α+z_β)²/(2δ²)` holds with an explicit, size-uniform error.* The kernel is exact; the covariance term is computable in closed form for the shuffle; this is a bounded, self-contained result.

**Principle theorem (ambitious but realistic — sufficiency direction).** *Under a `φ`-mixing or finite-susceptibility hypothesis on `{D_i}` with `ρ_n²→ρ_∞²∈(0,∞)`, `χ_n=o(n)`, and a Lindeberg condition, `Var(Ô_m)/O_∞² = 1/(ρ_∞²m)+o(1/m)+o(1)` and hence `m*(δ)→(z_α+z_β)²/(ρ_∞²δ²)`.* Provable via a Hájek finite-population CLT plus a covariance bound; delivers the *sufficiency* half of the principle in full generality.

The **necessity** half and the **critical-exponent floor law** (`floor ~ n^{a−1}` with `a` from the universality class) are the right *conjectures to state and test*, not to promise as proofs.

## 10. The strongest counterexample / failure regime to hunt for

The pilots found failure *at* criticality; the sharper test is to **decouple criticality from the operative condition (ii)** and to **isolate condition (i)**:

- **Extensive correlations away from any critical point.** Mean-field **Curie–Weiss Ising on the complete graph**, or any long-range model, has `χ_n=Θ(n)` across an entire low-temperature phase, not just at a point. If `m*` grows with `n` *throughout a phase* (not only at `T_c`), it shows the binding condition is **summable correlations**, not "criticality" — a cleaner counterexample than the lattice.
- **Conserved dynamics.** A conservation law (e.g. **symmetric simple exclusion**) induces slowly-decaying, non-summable correlations even off-critical; a failure here separates "phase transition" from "slow modes / conservation."
- **SNR-degenerate observable, correlations fine.** The **fixed-point count** on the *same* shuffle has `ρ_n²→0` (Poisson mean 1) while (ii) holds; predicted `m*~n` (inferability exponent `α=1`). This isolates condition (i) and is the clean control that (i) binds independently of (ii).
- **Heavy tails.** A per-component contribution with infinite variance breaks the CLT itself: the Gaussian constant is meaningless even if a "mean" exists — the failure mode of the *distributional* assumption rather than of (i) or (ii).

The decisive object is a **`2×2` factorial**: `{ρ² bounded, ρ²→0} × {χ summable, χ extensive}`, predicting scale-free certification in exactly one cell.

## 11. Experiments that discriminate between competing explanations

Four rival explanations for the observed scale-freeness, and the design that kills each:

| Explanation | Discriminating manipulation | Predicted outcome if true |
|---|---|---|
| (E1) "It's just that `O` is a sample **mean** (LLN/CLT)." | Keep `O` a mean but send `ρ²→0` (fixed points). | If E1 were the whole story it would still be scale-free; it is **not** (`m*~n`) ⇒ E1 insufficient. |
| (E2) "It requires **exchangeability**." | Use a non-exchangeable but short-range system (Ising high-T). | Already scale-free ⇒ **E2 falsified** (exchangeability not necessary). |
| (E3) "It requires **bounded per-component SNR**." | Hold `ρ²` bounded, make `χ` extensive (mean-field / critical). | Fails ⇒ E3 alone insufficient; (ii) binds. |
| (E4) "It requires **summable correlations** `χ=o(n)`." | The `2×2` factorial above. | Scale-free **only** in the (bounded-`ρ²`, summable-`χ`) cell ⇒ supports the two-condition principle. |

The `2×2` is the experiment that isolates the mechanism; the individual rows are the falsifiers already partly in hand.

## 12. The additional systems that are actually informative (each tied to a hypothesis)

Chosen to *test a stated hypothesis*, not for variety:

- **Curie–Weiss / complete-graph Ising** — tests whether the operative failure is *criticality* or *extensive correlation per se*. (Hypothesis: it's `χ=Θ(n)`, so failure spans a whole phase.)
- **Symmetric simple exclusion process (SSEP)** — tests whether a **conservation law** breaks (ii) without a phase transition, and extends certification to genuinely off-equilibrium, conserved dynamics.
- **Off-critical Ising with tunable `ξ`** (or the 1-D Ising line) — the correlation length is a **continuous knob**; sweeping `ξ` from `O(1)` to `O(L)` should trace the crossover `floor ~ ξ^{2−η}/n` and test the *quantitative* prediction, not just pass/fail.
- **Fixed-point / rare-event observable on the shuffle** — same system, degrade only `ρ²`; isolates axis (i) and pins the inferability exponent `α=1`.
- **Expander / random-regular-graph walk** — rapid mixing, bounded degree, summable correlations, no geometry; tests that it is dependence *range*, not spatial dimension, that matters.

Deliberately **excluded**: riffle, top-to-random, adjacent-transposition shuffles. All share the uniform stationary law and exchangeability, so they satisfy (ii) automatically and test nothing about the operative condition; they would only re-confirm mechanism-universality (a minor corollary).

## 13. What an expert mathematician would object to (and the honest reply)

- **"Certificate" is too strong.** Matching `E[O]` is a one-dimensional projection; equilibrium requires matching the full law. *Reply: agreed — reframe as a necessary-condition monitor for a chosen functional; state the projection explicitly.*
- **The alternative hypothesis is narrow.** Power is only against the relaxation family `ϕ=1−δ`; against `E[O]`-preserving deviations the power is nil. *Reply: define the alternative space precisely; claim detection only within it.*
- **The variance law is leading-order and semi-heuristic off equilibrium.** It uses the equilibrium image law inside an off-equilibrium computation and drops the inter-component covariance. *Reply: the mean is exact; the covariance term must be bounded (the anchor theorem does this for the shuffle) before any uniform claim.*
- **The Gaussian/finite-population CLT is an approximation**, worst exactly where it's cheapest (small `m`) and where `m→n` (FPC + floor). *Reply: quantify finite-`m` corrections; state the constant as `(1+o(1))`.*
- **The susceptibility condition is imprecise and the "iff" is unproven.** Which covariances, what normalisation, and is uniform integrability needed? *Reply: only sufficiency is claimed as a theorem target; necessity stays a conjecture; a Lindeberg condition is included.*
- **The Ising exponent rests on three sizes with possibly-unequilibrated `T_c` samples.** *Reply: acknowledged as approximate; the crossover experiment with measured autocorrelation time is required to make it quantitative.*

## 14. What an ISEF judge would ask that we cannot yet answer

- **"If your test passes, is the system mixed?"** — No; only this observable's mean is at equilibrium. *We do not yet have a principled recipe for choosing observables that make the screen closer to sufficient.*
- **"Can you predict, from a measurable quantity, whether your method will work on a given system?"** — We have the criterion `χ_n=o(n)`, but have not demonstrated predicting success/failure from an *a priori* measurement of susceptibility. *(A proposed experiment.)*
- **"Is it better than Gelman–Rubin / effective sample size, head-to-head?"** — Not yet benchmarked on a real sampler.
- **"Which claims are proven versus simulated?"** — Currently mostly simulated; the anchor theorem is provable but unproven.
- **"Is the method self-calibrating — can you estimate `ρ²` from the same small sample and thereby predict `m*` without knowing the system?"** — Plausibly yes (`ρ²` is estimable from the subsample), and this would be a strong selling point, *but we have not shown it works or quantified its own sampling error.*
- **"Name a real deployment where it beats current practice."** — We have candidates (large-lattice MCMC, distributed monitoring) but no demonstration.

These are the gaps that convert the project from "promising" to "grand-award"; every one is a concrete Phase-4/5 task, not a dead end.

---

## Proposed final project (built only from currently defensible claims)

**Title.**
*How Much Must You See to Know a Large System Has Settled? Size-Independent Certification of an Equilibrium Observable, and the Correlation Condition That Governs It.*

**Research question.**
For a large stochastic system observed only through a uniform random subset of its components, is the number of observations needed to certify — to a prescribed tolerance and confidence — that an intensive additive observable has reached its equilibrium mean independent of the system size, and what property of the components' correlation structure determines when this holds?

**Hypothesis.**
Size-independent certification holds exactly when two conditions are met: the per-component signal-to-noise ratio `ρ²=E[D]²/Var(D)` is bounded away from `0` and `∞`, and the summed component covariance (susceptibility) is sub-extensive, `χ_n=o(n)`. In that regime the required sample is `m*(δ)=(z_{1−α}+z_{1−β})²/(ρ²δ²)`, set by the single constant `ρ²`; when susceptibility becomes extensive (as at a critical point) the estimator variance is floored by the whole-system fluctuation `χ_n/n`, which subsampling cannot reduce, and size-independence is lost.

**Abstract (≈220 words).**
Many large stochastic systems are observed only in part, yet the questions asked of them are global: has an aggregate quantity settled to its equilibrium value? Building on a random-transposition model with an exactly known equilibrium, we study estimation of an intensive additive observable — the average displacement — from a uniform random sample of `m` of `n` components, and recast the problem as a test of whether the observable's mean has reached equilibrium. We show, analytically for the shuffle and a product chain and empirically across three systems of increasing coupling, that the sample needed to certify agreement to a relative tolerance `δ` is `m*≈(z_α+z_β)²/(ρ²δ²)`, independent of system size, where `ρ²` is the per-component signal-to-noise ratio; the equilibrium constant extends to a signal-to-noise *function* along the entire relaxation. A law-of-total-variance decomposition identifies the operative condition: size-independence requires not only bounded `ρ²` but **summable inter-component correlations**. When susceptibility is sub-extensive (independent, exchangeable, or short-range systems) the certificate is size-free; at a critical point, where susceptibility diverges, the estimator variance acquires an irreducible floor that no additional sampling removes, and the cost grows with size. The certificate is a necessary-condition monitor for one observable, not a guarantee of full equilibrium. The contribution is a precise sample-complexity law for equilibrium-observable certification and the correlation criterion that bounds its validity.

*Note on honesty of framing: the abstract claims a size-independent sample-complexity law (defensible: derived exactly for two systems, measured for a third), an SNR function along the trajectory (defensible: exact mean law + leading-order variance, pilot-confirmed), a variance-decomposition criterion (defensible: an identity), and a demonstrated critical-point failure (defensible: measured, mechanism identified). It does **not** claim full-equilibrium certification, a proven general theorem, or a benchmarked deployment — those are labelled as the work ahead.*
