# When can a glimpse certify equilibrium? An exact certification law and the correlation boundary that bounds it

**Mihika Singhania**
*Extension of "Partial Observation of Random Systems" (Research Science Initiative India). Draft manuscript; status tags in brackets mark what is proven, measured, or pending.*

---

## Abstract

A large stochastic system is often observed only in part, yet the question asked of it is global: has it reached equilibrium? We study when this can be certified from a uniform random sample of `m` of a system's `n` components, using a random-transposition shuffle as an exactly-solvable testbed with a known equilibrium. Writing the estimator's relative variance with the law of total variance splits it exactly into a sampling term `1/(ρ²m)` and an `m`-independent floor `R_n = Var_ν(O)/O_∞²` equal to the system's own equilibrium fluctuation. For the shuffle we obtain the certification law in closed form, including `Var(F) = (n+1)(2n²+7)/45` and floor constant `2/5`, and show that the sample needed to certify to relative tolerance `δ` is set by `δ`, not by `n`, precisely when the floor vanishes. The floor is a strict lower bound, so when it does not vanish, certification below a fixed tolerance is impossible at any size. We map the boundary with controlled experiments (bootstrap confidence intervals): certification is size-free exactly when the per-component signal-to-noise ratio is bounded and correlations are summable; a degenerate observable (fixed points; slope −1.00) or extensive correlation (mean-field below `T_c`; a constant floor) each break it. The sharpest case is a two-dimensional Ising critical point, where the order parameter is uncertifiable (a relative floor that no sample beats; exponent `a ≈ 0.86 ≈ 7/8`) while the energy density stays certifiable, so the observable, not the system, decides. An out-of-domain illustration on spatially correlated fields reproduces the same boundary and shows the practical payoff: where the floor does not vanish, an ordinary independence interval is badly overconfident (true coverage far below nominal) while the certificate correctly refuses. The certificate is a necessary-condition monitor for one observable, and this correlation criterion bounds its reach.

**Significance.** Statistical practice reports how confident an estimate is, but not whether the available data can justify the claim at all. This work isolates an exact, measurable quantity, the relative floor, that settles the second question for a partial observation of a large system: when it vanishes a fixed-size sample certifies a system of any size, and when it does not, no sample of any size can. The result is proved in closed form for an exactly-solvable model, its sharpest consequence is demonstrated at a physical critical point, and its practical value is shown where partial, correlated measurements are the norm.

---

## 1. Introduction

Decisions are increasingly made about whole systems from partial views of them, and the cost of misplaced confidence tends to fall on whoever did not produce the estimate. Many large stochastic systems evolve through repeated local interactions while an observer sees only a fraction of the state, yet the questions of interest are global: has the system equilibrated, and what is an aggregate statistic of its configuration? The foundational work of this programme (the RSI paper) studied the second question for the random-transposition shuffle and its average displacement, and showed that at equilibrium the relative uncertainty of the partial estimator is `1/√(2m)`, independent of the system size `n`.

That size-independence, taken alone, is a classical fact: finite-population sampling has long known that to estimate a mean to fixed relative precision one needs approximately a fixed number of observations regardless of population size (the finite-population correction). The contribution here is not size-independence. It is three things the classical picture does not supply. First, an exactly-solvable instance, a Markov chain whose equilibrium and every moment are theorems rather than estimates, which lets us write the certification law in closed form. Second, a reframing from estimation to certification: rather than reporting a value, we ask "has the observable reached equilibrium?" as a hypothesis test with a sample-complexity guarantee. Third, and most importantly, the correlation boundary. Classical finite-population theory assumes independent or exchangeable units; the entire question here is what happens when the units are correlated. We show that certification feasibility is governed by an exact, measurable quantity, the observable's relative floor, and that its failure is sharpest and most surprising at a critical point, where the same system is uncertifiable in one observable and certifiable in another.

**Relation to prior work.** The size-independence of relative precision is the finite-population correction. The decomposition of estimator variance into a sampling term and a between-configuration floor is the law of total variance and cluster sampling, and the inflation of variance by correlation is the effective-sample-size theme of spatial statistics and correlation-robust inference. Distortion of estimates by subsampling near criticality has been documented in neural dynamics (subsampling scaling). Our work synthesizes these into an exact, ground-truth-checkable certification law for equilibrating Markov systems and identifies the observable-dependent boundary at criticality. We differ from Markov-chain-Monte-Carlo convergence diagnostics, which use full-state summaries across multiple chains, by using a single chain and a random subset of components with a size-independent budget.

**Scope, stated up front.** The method certifies the mean of one chosen observable. Agreement there is necessary, but not sufficient, for full equilibrium, and the alternative it has power against is the one-parameter relaxation family, not an arbitrary non-equilibrium law. We regard this precision as a feature, and keep it visible throughout.

## 2. Setup and notation

Let a finite Markov chain have stationary law `ν` and `n` components with bounded per-component contributions `D_i(x) ∈ [0,B]`. The intensive observable is `O(x) = (1/n) Σ_i D_i(x)`, with equilibrium mean `O_∞ = E_ν[O]`. An observer draws a uniform random subset `S ⊆ [n]` of size `m` without replacement and forms `Ô_m = (1/m) Σ_{i∈S} D_i`. Write `V_x = (1/n) Σ_i (D_i(x) − O(x))²` for the within-configuration population variance, `χ_n[O] = n·Var_ν(O)` for the observable's susceptibility, and `R_n = Var_ν(O)/O_∞²` for the relative floor. The mixing level is `ϕ = E[O(x_t)]/O_∞ ∈ (0,1]`. The certification test is `H0: ϕ = 1` against `H1: |1−ϕ| ≥ δ`, at size `α` and power `1−β`, with `κ = z_{1−α} + z_{1−β}`. Our running systems are the random-transposition shuffle (`D_i = |π(i)−i|`, `O = AD`), a product chain of independent coins, the two-dimensional Ising model, and, for the illustration of Section 7, a Gaussian random field.

## 3. The floor identity and its consequences

**Lemma 1 (exact variance decomposition). [proven]** For a uniform `m`-subset without replacement,
> `Var(Ô_m) = [(n−m)/(m(n−1))]·E_ν[V_x] + Var_ν(O)`.

The first term is the sampling term, which shrinks as `1/m`; the second is an `m`-independent floor equal to the whole-system fluctuation, `Var_ν(O) = χ_n[O]/n`. (Proof in Appendix A.)

**Theorem 1 (floor lower bound). [proven]** For every `m`, `Var(Ô_m) ≥ Var_ν(O)`, hence `CV²(Ô_m) ≥ R_n`. Consequently, if `χ_n[O] = O(1)` the floor is `Θ(1/n)` and vanishes; if `χ_n[O] = Θ(n)` the relative floor is bounded below by a positive constant for all `n`, so there is a tolerance below which no sample of any size, at any `n`, can certify. The proof is one line from Lemma 1 (drop the non-negative sampling term). This converts "correlation limits certification" from folklore into a rigorous impossibility in the extensive-correlation limit.

## 4. The exactly-solvable anchor

Here `O = AD`, `O_∞ = AD_∞ = (n²−1)/(3n)`.

**Proposition 2 (single-card kernel and mean law). [proven]** The position of a tagged card is a lazy walk on the complete graph with `P_t(i→j) = 1/n + (δ_{ij}−1/n)(1−2/n)^t`. By linearity, `E[AD(t)] = (1−(1−2/n)^t)·AD_∞`, i.e. `ϕ(t) = 1−(1−2/n)^t`, which tends to `1−n^{−θ}` at `t = θ·t_mix` (Figure 2). This derives, analytically, the empirically-observed residual of the RSI relaxation curves: the collapse across sizes is asymptotic, not exact. (Appendix B.)

**Proposition 3 (footrule variance). [proven]** For uniform `π`, the Spearman footrule `F = Σ_i|π(i)−i|` has `Var(F) = (n+1)(2n²+7)/45` (inclusion-exclusion; Appendix C), verified by direct enumeration at `n=3`. Hence `Var_ν(AD) = Θ(n)` and the relative floor `R_n = (2n²+7)/(5(n−1)²(n+1))` has exact constant `c = lim n·R_n = 2/5`.

**Theorem 2 (exact certification law). [proven, MC-verified]** Combining Lemma 1 with Propositions 2 and 3 gives `Var(Ô_m)`, and hence `CV²(ϕ̂_m)`, in closed form. To leading order `CV(ϕ̂_m) = 1/√(2m)·(1+o(1))`, size-independent, with the `Θ(1/n)` floor of constant `2/5`. The exact expression matches Monte Carlo to 0.05-0.4% across `(n,m)` (Figure 1). A corollary gives the exact certification cost `m*(ε) = 5κ²n/[2(5ε²n − 2κ²)]`, which reduces to `κ²/(2ε²)` for large systems and diverges to infeasibility as `ε ↓ κ√R_n`, the exact instance of the Theorem-1 obstruction. That the independently-derived constant `2/5` matches the simulator's measured `n·R_n = 0.37-0.41` is a cross-check of theory and instrument.

## 5. Off-equilibrium and the general law

**Proposition 4 (trajectory signal-to-noise). [proven, leading order; measured]** At mixing level `ϕ`, `ρ²(ϕ) = ϕ/(3/2 − ϕ)` for the shuffle and `ϕ/(2 − ϕ)` for the coins, recovering `2` and `1` at equilibrium and matching simulation across `ϕ` (Figure 3). The equilibrium constant is thus the endpoint of a smooth trajectory law, not an isolated fact.

**Corollary 5 (certification cost). [measured]** Under a finite-population central limit theorem, the sample to certify a deviation `δ` is `m*(δ) → κ²/(ρ²δ²)`, size-independent; the simulator measures `m*δ² → 3.09` as `δ→0` (from 4.23 at `δ=0.20`), the finite-`δ` overshoot decreasing toward the constant, with bootstrap confidence intervals. The deployable (self-calibrated) test is the honest one; the idealized-null cost is an upper bound.

**Theorem 3 (general sufficiency). [proven under stated hypotheses]** For a triangular array of bounded contributions with per-component SNR converging to `ρ²_∞ ∈ (0,∞)`, summable susceptibility `χ_n = O(1)`, and a Lindeberg condition, `m*(δ) → κ²/(ρ²_∞δ²)`, independent of `n`. (Appendix D.) This is the size-free clause under checkable hypotheses; the shuffle and the coins are instances.

## 6. The correlation boundary

The proven results say size-free certification requires a bounded per-component SNR and a vanishing floor. Two experiments isolate these as the operative conditions, and a third pushes to a genuine critical point.

**The 2×2 discriminator (Figure 4). [measured, bootstrap CIs]** Holding the observable intensive and additive, we vary two structural axes independently: the per-component SNR (bounded vs degenerate) and the correlation structure (summable vs extensive). Cards' average displacement has `ρ²` flat in `n` (slope −0.00); the fixed-point fraction has degenerate `ρ² ∝ 1/n` (slope −1.00 [−1.02, −0.99]), so its certification cost grows with size. Curie-Weiss magnetization has `Var(M) ∝ 1/N` above `T_c` (slope −0.99 [−1.01, −0.97], certifiable) but a constant floor below `T_c` (slope +0.00 [0.000, 0.001]; floor ≈ 0.82), where extensive correlation makes certification impossible below a fixed tolerance, Theorem 1 realized. Only the (bounded-SNR, summable-correlation) cell is size-free; the other three break it, each for a distinct, measured reason.

**Criticality (Figure 5). [measured]** At the two-dimensional Ising critical point (Wolff ground truth, `L` to 80), the order parameter's relative floor `R_n[|M|]` is flat (slope +0.03; an approximately 26% certification floor that no sample beats), with susceptibility exponent `a[M] ≈ 0.86` consistent with `7/8`. The energy density at the same point has a floor that vanishes as `n^{−0.84}` and is therefore certifiable. The operative quantity is the observable's relative floor `R_n`, which combines its susceptibility and its mean; the raw susceptibility exponent alone can mislead, because at criticality the order parameter's own mean collapses. The choice of observable, not the system, decides whether a small sample can certify.

**The headline (Figure 6).** Across seven orders of magnitude in `n`, the certification cost `m*(δ)` is flat for a summable-correlation observable, against the rising curve intuition expects, until, at a critical order parameter, it is unbounded. That flatness with a sharp break is the result in one picture.

## 7. An illustrative application: correlated fields and the overconfidence the floor catches [measured; illustrative]

*This section is a demonstration, not a load-bearing result. Its purpose is to show that the floor identity of Section 3 governs a system built by none of the machinery above, and that reading the floor has immediate practical value. The theory of Sections 3 to 6 stands regardless of what follows.*

Spatially correlated fields are the natural home of partial observation: a district's air quality, a region's soil moisture, a population's exposure, are each summarised by an average that someone certifies from a handful of sensors. We take the canonical model of such a field, a stationary Gaussian random field on an `L×L` grid with tunable correlation length `ξ`, positive mean `O_∞`, and observable `O` the spatial mean, and read the floor `R_n = Var_ν(O)/O_∞²` off it directly, with no shuffle and no Ising structure.

**The floor obeys the theory out of domain (Figure 7a).** When the correlation length is fixed as the grid grows (summable correlation), the floor falls as `R_n ∝ n^{−0.98}`, matching the `Θ(1/n)` law of Theorem 3 in a system that shares none of the anchor's structure. This decay is the non-trivial confirmation, since nothing about the field was tuned to produce it. When the correlation length instead scales with the system (`ξ = L/8`, extensive correlation), the floor is flat in `n` at `√R_n ≈ 13%` (slope −0.01), the obstruction of Theorem 1 made visible in a spatial field rather than a magnet. The floor here is set by the field's effective sample size, which the correlation length controls, not by Gaussianity; the Gaussian field is simply the cleanest case in which it can be computed and checked exactly.

**The payoff: the floor names an overconfidence the standard interval hides (Figure 7b).** Fix a strongly-correlated field and let an analyst estimate the region's mean from `m = 200` sites with the ordinary 95% interval that assumes independent samples. As the correlation length grows, that interval keeps reporting a confident half-width of a few percent, yet its true coverage of the equilibrium mean falls from 91% to about 12%: an interval that should contain the truth 19 times in 20 in fact does so about 1 time in 8, with no outward sign. An agency clearing a district as safe at nominal 95% confidence would, in this regime, be wrong roughly seven times in eight. The floor is exactly the quantity that flags this, since `√R_n` has grown past the demanded tolerance, so certification to that tolerance is impossible and the honest verdict is to refuse.

**What this does and does not show.** The refusal above is computed from the floor estimated across many independent realisations of the field, an ensemble the single-field analyst does not possess. So the demonstration establishes that the floor exists in a realistic correlation structure and that ignoring it produces a concrete, quantified error; it does not yet show the floor can be read from one field in the wild. Estimating `R_n` from a single realisation with a guarantee is precisely the open problem of Section 10, and validation on measured fields, air-quality monitor networks with their documented spatial correlation, is the first target of the deployable instrument. A Gaussian, stationary field is a model; real concentration fields are non-stationary and heavier-tailed, which the same effective-sample-size logic should cover but which remains to be checked.

## 8. Discussion

The certificate is a necessary-condition monitor for one observable's mean, made cheap (`Θ(1/δ²)` reads, independent of `n`) where a full audit costs `Θ(n)`. It is naturally self-calibrating: `ρ²` and `R_n` can be estimated from the observed subsample, so the method can declare its own infeasibility from the data when the floor is high, as it does for Curie-Weiss below `T_c` and for the correlated field of Section 7. Its relation to MCMC convergence diagnostics is complementary: it is cheaper and size-free in the summable regime, and it fails exactly where those diagnostics struggle (criticality), for the same structural reason, an explanatory rather than merely cautionary statement. The limits are honest and three-fold: one observable's mean (necessary, not sufficient); one alternative family; and a Gaussian test that must be replaced by an exact or bootstrap null at very small `m`.

## 9. Methods

Systems, a leak-audited observation layer, a batched generation engine, frequentist and Bayesian inference, per-configuration diagnostics, and bootstrap-CI scoring are implemented in a Python instrument (`rsi_partial`) validated against analytical ground truth (exact shuffle moments; Onsager energies) and against the Phase-3 pilot findings, with a seed protocol for reproducibility. Ising equilibrium ensembles use a Wolff cluster sampler (so critical slowing does not bias ground truth); off-equilibrium relaxation uses single-site Glauber. The Gaussian random field of Section 7 is generated by filtering white noise to a chosen correlation length, normalised to unit marginal variance by a fixed (not per-realisation) constant, with the floor and the interval coverage estimated over independent realisations. The 2×2, criticality, and application figures use matched samplers and bootstrap confidence intervals on all fitted slopes.

## 10. Conclusion

From a few random pieces one can certify a huge system's equilibrium at a cost set by the demanded precision, not the system's size, but only for observables whose relative floor vanishes. That floor equals the system's own equilibrium fluctuation; it is an exact lower bound; and at a critical point it makes the order parameter uncertifiable while other observables remain fine. The size-independence is classical; the exact law, the certification framing, and the observable-boundary at criticality are the contribution, and a correlated-field illustration shows the same boundary governs realistic data and prevents a concrete inferential error.

**Still open [pending].** A publication-grade critical-exponent fit (larger `L`, corrections-to-scaling, CI on the slope, in the reproducible instrument); the full eight-prediction verdict table against the pre-registered bands; a data-only estimator of the floor with a guarantee, validated on a live correlated stream; and conserved-dynamics (exclusion processes), where slow modes may demand a dynamic condition beyond the static floor.

---

## Related work (expanded)

Finite-population correction and the sample-size-independent-of-`N` fact are standard survey sampling. Effective sample size under correlation and correlation-robust inference are central to spatial statistics (e.g. spatial-correlation-robust inference; optimal spatial subsample size for variance estimation). Subsampling bias near criticality has been studied in neural systems (subsampling scaling). MCMC convergence diagnostics (Gelman-Rubin; nested `R̂`; the convergence-diagnostics literature) assess convergence from full-state summaries. The random-transposition shuffle, its mixing (Diaconis-Shahshahani), and the footrule (Diaconis-Graham) are classical; the exact footrule variance is re-derived here in-house.

## Appendices (proof sketches)

**A. Lemma 1.** Condition on `x`; the SRSWOR sample mean has `E_S[Ô_m|x]=O(x)` and `Var_S(Ô_m|x)=V_x/m·(n−m)/(n−1)`; apply the law of total variance.

**B. Proposition 2.** The one-step matrix has off-diagonal `2/n²`; on the mean-zero subspace its eigenvalue is `1−2/n`, giving the kernel. The mean law follows because `Σ_j|j−i|P_t(i→j) = (1−r)·(1/n)Σ_j|j−i|` (the `j=i` term vanishes).

**C. Proposition 3.** `E[F²] = Σ_i E[D_i²] + Σ_{i≠j} E[D_iD_j]`; the off-diagonal sum over distinct pairs equals `P² − 2Q + R₂` with `P = n(n²−1)/3`, `R₂ = n²(n²−1)/6`, `Q = Σ_i g(i)²`, `g(i) = i² − (n+1)i + n(n+1)/2`; evaluating `Q = n(n−1)(n+1)(7n²−8)/60` and subtracting `(E[F])²` gives `Var(F) = (n+1)(2n²+7)/45`.

**D. Theorem 3.** Lemma 1 gives `CV² = (1/ρ_n²)/m·(1+O(m/n)) + χ_n/(nO_∞²)`; the hypotheses send the first coefficient to `1/ρ²_∞`, the floor to `O(1/n)`, and (Lindeberg/Hájek) make the detection statistic asymptotically normal, inverting the design equation to `m*→κ²/(ρ²_∞δ²)`.

---

### Figures
- **Figure 1** `fig1_exact_cv2.png`, exact certification law vs Monte Carlo.
- **Figure 2** `fig2_trajectory.png`, exact mean trajectory and its asymptotic collapse.
- **Figure 3** `fig3_rho2.png`, per-piece SNR along the trajectory.
- **Figure 4** `fig4_discriminator.png`, the 2×2 discriminator (centerpiece).
- **Figure 5** `fig5_criticality.png`, Ising at criticality: order parameter vs energy.
- **Figure 6** `fig6_headline.png`, flat certification cost with a critical break.
- **Figure 7** `fig7_application.png`, illustrative application on a correlated field: (a) the floor vanishes (short-range) or stays flat (long-range); (b) as correlation grows the naive interval's true coverage collapses while the certificate floor rises past the demanded tolerance and refuses.
