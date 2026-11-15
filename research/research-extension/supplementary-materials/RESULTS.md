# RESULTS — RSI Phase-5 Simulator

*Running results log. Every number is tagged: **[EST]** established mathematically ·
**[EMP]** established empirically (measured here or in Phase 3) · **[CONJ]** conjectured ·
**[PCV]** pending computational validation. Theory lives only in
`analysis/predictions.py`; measurements are produced by the estimators/engine and compared
to theory only in the scoring layer. This file deliberately **under-claims**: where a
result is a re-confirmation of Phase 3, it says so; where it is new, it is stated with a CI.*

---

## Honest status

**Stage-0 machinery is built and validated; one genuinely new result (the `a=1` hard
floor) has been produced and scored. The full P1–P8 campaign has not been run.**

- The instrument (systems, veil, engine, inference, diagnostics, analysis, storage) is
  complete and passes an analytical + regression + structural gate (**74 tests green**).
- The `a=0` measurements (`ρ²`, size-free `m*`, cards floor `∝1/n`) **re-confirm Phase 3**
  — they are the validation base, not the contribution.
- **New:** the Curie–Weiss `a=1` extensive-susceptibility **hard floor** (CX1) is measured
  with a bootstrap CI and scored SUPPORT — the same system flips from `a=0` (above `T_c`)
  to `a=1` (below). This is the first falsification-grade result beyond the paper.
- **Not yet done:** the full P1–P8 verdict campaign; Ising critical-exponent (`a=7/8`) at
  scale; systems S5 (tunable-ξ) and S6 (SSEP). The instrument is ready to run them.

---

## Validation gate (re-confirmations of known truth)

Measured with the correct rapid/direct samplers (not the Phase-3 finite-burn-in proxy),
fixed seeds, R = 4000–15000.

| Quantity | Measured | Theory | Tag |
|---|---|---|---|
| Cards `AD_∞` (n=512) | 170.60 ± 0.08 | `(n²−1)/3n` = 170.67 | [EST] |
| Cards mean law `E[AD(t)]` | matches across t | `(1−(1−2/n)^t)·AD_∞` | [EST] |
| Cards `ρ²(1)` (n=256…8192) | 2.000–2.008 (flat in n) | 2 | [EST] |
| Exact estimator variance | matches < 6% | RSI App. E identity ★ | [EST] |
| Cards **floor slope** `d log R_n/d log n` | **−1.037 [−1.117, −0.957]** | `a−1 = −1` (a=0) | [EMP] |
| Chain `ρ²(ϕ)`, `Var(ϕ̂_m)` | match < 3–6% | `ϕ/(2−ϕ)`, `ϕ(2−ϕ)/m` | [EST] |
| Ising energy/\|M\| (Wolff) | < 2% at T=2.0 | Onsager | [EST] |
| Glauber (checkerboard) stationarity | `<E>` preserved (diff/se 0.2) | detailed balance | [EST] |

---

## Detection cost `m*` — honestly bracketed (item 3)

Cards, n=4096, calibrated (oracle) null, **bootstrap 95% CIs**, δ-sweep. Target
`κ²/ρ² = 3.091`. The finite-δ overshoot shrinks toward the target as δ→0 (the `σ₁>σ₀`
correction), so the constant is **measured**, not merely bracketed. [EMP]

| δ | oracle `m*` [95% CI] | `m*δ²` [CI] | held-out Type-I | deployable `m*δ²` (Type-I) |
|---|---|---|---|---|
| 0.20 | 106 [97, 110] | 4.23 [3.87, 4.40] | 0.045 | 3.85 (0.051) |
| 0.15 | 180 [161, 190] | 4.05 [3.62, 4.27] | 0.045 | 3.79 (0.051) |
| 0.10 | 367 [349, 395] | 3.67 [3.49, 3.95] | 0.056 | 3.49 (0.061) |
| 0.075| 605 [557, 694] | 3.40 [3.13, 3.90] | 0.053 | 3.17 (0.060) |

- **Two nulls, reported separately.** The *oracle* (calibrated) null assumes you know the
  equilibrium sampling law; its Type-I is measured on a **held-out** H0 pool (so ≈0.05 is a
  real check, not tautological). The *deployable* Gaussian null is what you'd have in the
  field; its Type-I **inflates** to ~0.06 as `m` grows (the finite-`m` self-normalized
  effect). The headline size-free story holds under both; the deployable path is the honest
  one for a real deployment.
- Size-free check [EMP]: at fixed δ, `m*` is flat in `n` (fitted exponent `b≈0`); the CIs
  overlap across n. (Earlier "168/168/168" was a coarse-grid artefact; with CIs the honest
  statement is "consistent with `b=0`," not "identical.")

---

## NEW: the `a=1` hard floor (Curie–Weiss, CX1 / P2) — item 4

Mean-field Ising, exact equilibrium sampler, floor `Var_ν(M)` vs `N ∈ {100…3000+}`,
slope = `a−1` with a bootstrap CI. **[EMP, scored SUPPORT for the CONJ]**

| Regime | `a` estimate | slope (a−1) [95% CI] | P2 verdict |
|---|---|---|---|
| **T = 0.7 (below `T_c`)** | **1.005** | +0.005 [+0.004, +0.005] | **SUPPORT** (a=1) |
| T = 2.0 (above `T_c`) | 0.004 | −0.996 [−1.002, −0.990] | SUPPORT (a=0) |

- The **same system** is `a=0` above `T_c` (self-averaging, `Var(M)∝1/N`) and `a=1` below
  (`Var(M)→m0²=Θ(1)`). Below `T_c` the floor is a genuine barrier: `Var(M̂_m)/floor ≈ 1.0`
  for every m from 16 to 1024 — **subsampling cannot reduce it**. This is the predicted
  hard, non-vanishing precision floor; a flat `m*` here would have *falsified* the refined
  conjecture, and it did not.
- (Caveat: below-`T_c` the slope CI is very tight and can exclude *exactly* 0 by mild
  finite-`N` drift as `Var(M)→m0²`; the verdict is by CI-within-band, and the asymptotic
  slope is 0.)

---

## Deployability (item 5)

- **Self-calibration recovers the floor.** From repeated subsamples alone (no system
  knowledge), fitting `CV²(m)=1/(ρ²m)+R_n` recovers `ρ̂²=1.978` (cards, true 2) and, for
  Curie–Weiss below `T_c`, `R̂_n≈1.0` (the hard floor) → predicted `m*=∞`: the method
  **declares certification infeasible from the data**. [EMP]
- **Equilibration gate detects under-burn-in.** A two-initial-condition check flags a chain
  that has not relaxed and refuses to record it (the Phase-4 negative control). Below `T_c`,
  local Glauber is correctly flagged as slow even at 120 sweeps — which is exactly why Wolff
  is the ground-truth sampler. [EMP]

---

## Corrections to earlier over-claims (from the audit)

- **Access audit is now real.** The engine and detection driver select subsets through
  `veil.reveal_batch`; the read count is derived from the returned array, not hard-coded.
- **Memory claim corrected.** The honest statement is **peak memory `O(batch·n)`, never
  `O(R·n)`** — not "`O(Rm)`" (the configuration must be realized to be sampled; only the
  *selection* is `O(batch·m)`).
- **Ising off-equilibrium now exists.** Checkerboard-vectorized Glauber (~1400 sweeps/s at
  L=96) + a relaxation generator; the `phi<1`-silently-returns-equilibrium bug is closed.
  Ising detection certifies the **energy density** (nonzero `O_∞`), for which `ϕ` is
  well-posed; `|M|` is used only for the equilibrium floor exponent.

---

## Not yet established (Phase-5 targets — deliberately unrun)

All **[PCV]**: the Ising critical floor exponent (`M: a=7/8`, `E: a≈0`) at scale with the
Wolff ground truth; the full P1–P8 verdict campaign; the tunable-ξ crossover (S5); the
static-vs-dynamic gap under conservation (SSEP, S6). S5/S6 register onto the existing
interfaces without pipeline changes; the scorer already emits support/falsify/inconclusive
verdicts.
