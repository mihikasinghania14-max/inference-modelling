# Phase 7 — Paper Skeleton & ISEF Submission Plan

*The integration phase. This is the scaffold for the sequel paper: section-by-section structure with every result slotted by status, so Phase-5 numbers drop into marked placeholders as they finish. It is a **skeleton**, not the final draft — the full paper needs the complete P1–P8 campaign (only the `a=1` floor is done so far). Status tags: **[READY]** proven or measured now; **[PARTIAL]** some evidence, needs the full run; **[PENDING]** awaiting Phase-5. Sources: `phase1…phase6` project docs, `CRITIQUE_simulator.md`, `RESULTS.md`, the RSI paper.*

---

## 0. Front matter

- **Working title:** *How Little Must You See? Size-Independent Certification of Equilibrium in Partially Observed Random Systems.*
  (Alt, more technical for the ISEF board: *A Sample-Complexity Law for Certifying Equilibrium Observables, and the Correlation Condition That Bounds It.*)
- **Author:** Mihika Singhania. Advisor acknowledgement: Dr Arvind Ayyer (IISc); RSI India.
- **Framing sentence (opening):** the RSI paper proved that *at equilibrium* the relative uncertainty of the average-displacement estimator is `1/√(2m)`, independent of system size; this paper drops "at equilibrium" and asks what is required for that size-independence to hold, and where it fails.
- **One-line contribution:** partial-observation certification of an intensive equilibrium observable has an exact information floor equal to the whole-system fluctuation; it is size-free **iff** the observable's susceptibility is bounded, with the per-component signal-to-noise ratio setting the constant — proven for the shuffle, and confirmed to break at extensive susceptibility.

### Abstract (draft — defensible claims only) [READY to draft, numbers to finalize]

> Many large stochastic systems are seen only in part, yet the questions asked of them are global: has an aggregate settled to equilibrium? We recast estimation of an intensive additive observable from a uniform random sample of `m` of `n` components as a test of whether the observable's mean has reached equilibrium, and ask when the required sample is independent of `n`. A law-of-total-variance identity splits the estimator's relative variance exactly into a sampling term `1/(ρ²m)` and an `m`-independent floor `R_n = χ_n/(nO_∞²)` equal to the system's own equilibrium fluctuation. This yields three results. (i) For the random-transposition shuffle we derive the exact certification law in closed form, with `Var(F)=(n+1)(2n²+7)/45` giving floor constant `2/5`, so `CV = 1/√(2m)` to leading order and the certification cost `m* ≈ κ²/(ρ²δ²)` is size-free. (ii) A general sufficiency theorem shows this holds for any bounded, summable-susceptibility observable obeying a finite-population CLT. (iii) The floor is a strict lower bound, so at extensive susceptibility (`a=1`) certification below a fixed tolerance is impossible at any `n`; we confirm this on the mean-field Ising model, where the same observable is size-free above `T_c` and floor-limited below. The certificate is a necessary-condition monitor for one observable, and the correlation criterion `χ_n=o(n)` bounds its validity. [FINAL numbers: critical-exponent behaviour PENDING.]

---

## 1. Introduction [READY]

- The partial-observation problem; the two mature theories it joins (finite Markov chain equilibrium + finite-population sampling). *(Reuse RSI §1 framing.)*
- What the RSI paper established (`1/√(2m)`, `ρ²=2`, `m*=z²/2ε²`, Prop 6.1). One paragraph, cite it as the foundation.
- The gap this paper closes: RSI answered "what is the value" at equilibrium but assumed the answer to "has it settled"; and Prop 6.1's exchangeability tacitly hid a second condition. State the three contributions (identity+floor, exact shuffle law, the correlation boundary).
- Scrupulous scope sentence: we certify the *mean of one observable* — necessary, not sufficient, for full equilibrium.

## 2. Setup and notation [READY]

Pull the fixed notation block from `phase6_theory.md`: chain, `ν`, `D_i`, `O`, `O_∞`, `Ô_m`, `V_x`, `ρ²(ϕ)`, `χ_n[O]`, `R_n`, `ϕ`, the certification test `H0: ϕ=1` vs `H1: |1−ϕ|≥δ`. Define the three systems (shuffle, product chain, Ising) as the running examples.

## 3. The variance identity and the floor [READY — proven]

- **Lemma 1** (exact decomposition, `phase6` Part I): `Var(Ô_m) = FPC·E_ν[V_x]/m + Var_ν(O)`.
- **Theorem 1 = T3** (floor lower bound + `a=1` impossibility). This is the conceptual spine; lead with it.
- Figure 3.1: schematic of "sampling term shrinks, floor does not" — the two-term picture.

## 4. The random-transposition anchor [READY — exact]

- **Prop 4.1** single-card kernel `P_t(i→j)=1/n+(δ−1/n)(1−2/n)^t` (proof).
- **Prop 4.2** exact mean law `E[AD(t)]=(1−(1−2/n)^t)AD_∞`; Corollary: derivation of the RSI Fig-4 residual (`ϕ(θ)≈1−n^{−θ}`).
- **Prop 4.3** exact moments incl. `Var(F)=(n+1)(2n²+7)/45` (inclusion–exclusion proof, `Q` computed).
- **Theorem 2 = T1** exact `CV²` closed form; floor `R_n=(2n²+7)/(5(n−1)²(n+1))`, constant `c=2/5`; **Corollary T1′** exact certification cost with infeasibility threshold.
- Figure 4.1: `CV²` exact vs Monte Carlo (ratios ≈1.00); Figure 4.2: mean-trajectory collapse `1−n^{−θ}`.
- Table 4.1: exact vs measured (`AD_∞`, `ρ²`, floor slope, `c=2/5` vs measured `n·R_n=0.37–0.41`).

## 5. Off-equilibrium and the general law [READY leading-order; T2 conditional]

- **Prop 5.1** off-equilibrium SNR `ρ²(ϕ)=ϕ/(3/2−ϕ)`, precision `Var(ϕ̂)=ϕ(3−2ϕ)/2m` (leading order).
- **Corollary 5.2** detection cost `m*_det=κ²/(2δ²)`, size-free (under the Hájek CLT; note the oracle-vs-deployable null distinction from `RESULTS.md`).
- **Theorem 3 = T2** general sufficiency under bounded contributions + `ρ²→ρ²_∞` + summable `χ` + Lindeberg.
- Figure 5.1: `ρ²(ϕ)` curve, cards + chain, measured vs theory. Figure 5.2: `m*δ²` → `κ²/ρ²` as `δ→0` with bootstrap CIs (from `RESULTS.md`).

## 6. The correlation boundary [PARTIAL → PENDING]

The heart of the *new* science; the axis RSI never explored.
- **6.1 The exponent framework** (`phase4` §1): `χ_n≍n^a`; three regimes `a=0` / `0<a<1` / `a=1`; `δ_floor≍n^{(a−1)/2}`; `I_max≍n^{1−a}`.
- **6.2 `a=1` extensive susceptibility — the hard floor.** [READY] Curie–Weiss: same system `a=0` above `T_c` (slope −0.996), `a=1` below (slope +0.005), floor unremovable (`Var(M̂_m)/floor≈1.0` for all m). Confirms Theorem 1's impossibility. Figure 6.1: the two-regime floor slopes with bootstrap CIs. This is the falsification-grade result.
- **6.3 `0<a<1` — criticality.** [PENDING] 2-D Ising at `T_c`: predicted `a=7/8` (magnetization) vs `a≈0` (energy). The one genuinely conjectural claim; needs the Wolff-ground-truth exponent fit at scale. Figure 6.2 (placeholder): floor slope vs observable at `T_c`.
- **6.4 The 2×2 discriminator.** [PARTIAL] `{ρ² bounded/degenerate} × {χ summable/extensive}` — scale-free only in one cell. Cards + fixed-points + Curie–Weiss give three cells; the fourth is quick.
- **6.5 Beyond static: conservation (SSEP).** [PENDING] whether a conservation law breaks certification despite summable static `χ` — the CX5 question, may reveal a needed dynamic condition.

## 7. Discussion [READY to draft]

- Necessary-not-sufficient; observable choice; self-calibration (the method estimates `ρ²`, `R_n` from its own subsample and declares infeasibility — from `RESULTS.md`).
- Relation to MCMC convergence diagnostics: cheaper (`Θ(1/δ²)` reads vs `Θ(n)`), size-free in the summable regime, and it fails exactly where Gelman–Rubin does (criticality) for the same structural reason — an explanatory, not merely cautionary, contribution.
- Limits and honest scope (the three RSI-style limitations, updated).

## 8. Methods [READY]

The `rsi_partial` instrument: systems / veil (audited access) / engine / inference / diagnostics / analysis; Wolff ground truth vs Glauber relaxation; seed protocol; bootstrap CIs; the equilibration gate. Reproducibility statement. *(Condense `Phase5_architecture.md` + `RESULTS.md` "corrections" section.)*

## 9. Conclusion [READY]

Restate the final claim (`phase4` §10 / `phase6`): exact floor identity, size-free iff bounded susceptibility, the shuffle law in closed form, and the proven extensive-`χ` impossibility — with the critical-exponent rate as the open frontier.

### Appendices
- A. Single-card kernel spectral proof. B. Footrule variance (inclusion–exclusion, `Q`). C. Exact `CV²` derivation + MC verification. D. Hájek CLT hypotheses for T2. E. Curie–Weiss floor measurement. F. The `rsi_partial` API + reproducibility.

---

## Figure & table inventory (build order)

| # | Content | Status | Source |
|---|---|---|---|
| 3.1 | Two-term floor schematic | [READY] | Lemma 1 |
| 4.1 | Exact `CV²` vs MC | [READY] | `exact_cv.py` verification |
| 4.2 | Mean-trajectory `1−n^{−θ}` collapse | [READY] | Phase-3/6 |
| 4.3 | Table: exact vs measured constants (`c=2/5`) | [READY] | `phase6` + `RESULTS.md` |
| 5.1 | `ρ²(ϕ)` cards + chain | [READY] | Phase-3 |
| 5.2 | `m*δ²→κ²/ρ²`, bootstrap CIs | [READY] | `RESULTS.md` |
| 6.1 | Curie–Weiss two-regime floor slopes | [READY] | `RESULTS.md` (a=1) |
| 6.2 | Ising floor exponent vs observable at `T_c` | [PENDING] | Phase-5 run |
| 6.3 | 2×2 discriminator grid | [PARTIAL] | Phase-5 run |
| 6.4 | `m*` vs `n` flat (a=0) / bending (critical) heatmap | [PENDING] | Phase-5 run |

## Pre-submission checklist (what must finish before the paper is done)

1. **[PENDING]** Ising critical floor exponent at scale (Wolff GT): magnetization `a≈7/8`, energy `a≈0`, each with a bootstrap-CI slope fit. *This is the single most important remaining run — it is the only claim that can still surprise us.*
2. **[PENDING]** Complete the 2×2 discriminator (add the degenerate-SNR × summable cell = fixed-points, already built) and the P1–P8 verdict table.
3. **[PARTIAL]** Tighten `m*δ²→κ²/ρ²` with more `δ` points to show the `δ→0` limit cleanly.
4. **[OPTIONAL]** SSEP (S5/S6) — strengthens the discussion; not required for the core claims.
5. Finalize `Var(F)` footrule constant citation vs in-house derivation (done — derived).

## ISEF submission materials (parallel to the paper)

- **Abstract (250 words):** the §0 draft above, tightened once Fig 6.2/6.4 land.
- **Board layout:** hero = the flat `m*`-vs-`n` line with the critical break (Fig 6.4); the identity as the "why"; the Curie–Weiss two-regime slopes as the falsification-grade result; **the PROBE game as the live interactive** (judges play it — the flat trace *is* the thesis).
- **The two "genuinely new" claims to foreground for judges:** (a) the floor is an *exact information limit*, proven; (b) size-free inference is bounded by *summable correlations*, not exchangeability — and we found the boundary, at criticality.
- **Anticipated judge questions + answers:** already drafted in `phase3_interpretation.md` §14 and the critique's ISEF-hat — port them.
- **Intellectual-ownership narrative:** RSI paper = foundation; this work = the trajectory extension, the exact floor theorem, the correlation boundary, the instrument, and the game. The theory and the simulator cross-validate (`c=2/5` predicted and measured; `a=1` proven and measured).

---

### Status summary

Ready to write in full **now**: §§1–5, §6.1–6.2, §7–9, and all [READY] figures — that is the majority of the paper, resting on proofs and already-measured results. The paper's *conclusions do not depend on the pending runs*; those sharpen §6.3–6.5 (the critical exponent and conservation), which are the frontier, not the foundation. Recommended: draft the [READY] paper now (it stands on its own as a complete result), and treat §6.3–6.5 as the "future work / in progress" that the full campaign will fill.
