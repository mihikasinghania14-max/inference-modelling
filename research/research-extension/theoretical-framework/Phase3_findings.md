# Phase 3 — Pilot Experiments: Findings and Simulator Requirements

**Project direction (locked):** *A scale-free partial-observation equilibrium certificate* — turn the RSI scale-invariance principle into a **practical, system-size-independent convergence diagnostic**: certify whether a large stochastic system has reached equilibrium from a random subset of its components, in a number of observations that does not grow with the system. *Mostly applied* (a working diagnostic validated across real systems and benchmarked), *partly theoretical* (the off-equilibrium SNR law and detection threshold).

**What Phase 3 is.** Deliberately lightweight, throwaway simulations run to *discover the patterns* before building the instrument — not the Phase-2 simulator itself. Every result below is from a short vectorised script; the point is to learn exactly what the real simulator must do. Runtimes were seconds to a couple of minutes each.

**Three testbeds, chosen to span coupling strength:**

| Testbed | Coupling | Ground truth | Role |
|---|---|---|---|
| Random transposition (cards) | weak, non-local | exact | theory anchor (RSI system) |
| Product / hypercube chain | none (independent coords) | exact | cleanest applied case |
| 2D Ising–Glauber | strong, local | long-run proxy | flagship applied + breakdown test |

---

## 1. What the pilots found

### 1.1 The off-equilibrium SNR law holds along the *entire* trajectory (not just at equilibrium)

The RSI paper's per-component signal-to-noise ratio `ρ²=2` is an equilibrium object. The pilots confirm it is the endpoint of a smooth **off-equilibrium law** `ρ²(ϕ)`, where `ϕ = AD(t)/AD_∞ ∈ (0,1]` is the mixing level:

- **Cards:** `ρ²(ϕ) = ϕ/(3/2 − ϕ)` — empirical vs predicted across `ϕ = 0.47…1.0`: `0.52/0.45, 1.04/0.91, 1.34/1.30, 1.76/1.76, → 2.0/2.0`. (Panel A, blue.)
- **Product chain:** `ρ²(ϕ) = ϕ/(2 − ϕ)` — empirical vs predicted match to **three decimals** across `ϕ`: `0.176, 0.333, 0.538, 0.818, 1.000`. (Panel A, orange.)

Both recover their equilibrium constant (`2` and `1`) at `ϕ=1`. This is the trajectory-wide generalisation of RSI Theorem 5.3 — new knowledge the paper does not contain.

### 1.2 The precision law is genuinely system-size-independent

`Var(ϕ̂_m)` predicted by the single-component law, empirical/predicted ratio:

- **Cards** `= ϕ(3−2ϕ)/(2m)`: ratio `0.80–1.04`, drifting up to `1.0` as `n` grows from 128→2048 — the predicted **O(1/n) covariance correction** from weak inter-card dependence.
- **Product chain** `= ϕ(2−ϕ)/m`: ratio `0.97–1.03`, essentially **1.00 and flat** across `n = 256 → 65 536` (no covariance correction, as independent coordinates predict).

(Panel B.) Scale invariance is confirmed over **2.4 orders of magnitude** in `n`.

### 1.3 The applied headline: certifying equilibrium costs `m* ≈ c/δ²`, independent of `n`

Detection problem `H0: ϕ=1` (mixed) vs `H1: ϕ=1−δ`, minimum `m` for 80% power at 95% one-sided. Results (Panel C):

| δ | cards m\* (n=1024 / 4096) | chain m\* (n=4096 / 65 536) |
|---|---|---|
| 0.20 | 90 / 90 | 160 / 160 |
| 0.15 | 150 / 150 | 280 / 280 |
| 0.10 | 315 / 330 | 640 / 620 |
| 0.075 | 540 / 555 | — |

**m\* is identical across system sizes**, and `m*·δ² → (z_α+z_β)²·ρ²(1)⁻¹`: cards cluster at `≈3.1 = (z_α+z_β)²/2`, chain at `≈6.3 = (z_α+z_β)²`. The number of components you must observe to certify equilibrium to tolerance `δ` is a **constant of the question**, and the constant is exactly the per-component SNR. This is the RSI cost law `m*=z²/(ρ²ε²)`, now demonstrated for the *detection/certification* problem the paper motivated but never solved, and shown to transfer across systems.

### 1.4 Where it breaks — the applicability boundary (the most interesting applied result)

2D Ising–Glauber at equilibrium, variance of magnetisation estimated from `m` random sites (Panel D):

- **High-T (T=3.5, short correlation length):** for `m ≪ n`, `Var(M̂_m) ≈ Var(site)/m`, with ratios `→ 1.0` as `n` grows (`1.20, 1.06, 0.97` for `m=16` at `n=256,1024,2304`). The magnetisation floor `Var(M_full)` decays `∝ 1/n`. **Scale-free certification holds.**
- **Near-Tc (T=2.27, diverging susceptibility):** `Var(M̂_m)` is pinned to an **irreducible floor** — `0.58, 0.56, 0.56` at `m = 16, 64, 256`. Observing more sites **does not help** (ratios vs `Var(site)/m` blow up to 100+), and the floor decays only `∝ n^{-1/8}`, far too slowly.

The certificate works precisely when correlations are short-ranged (finite susceptibility ⇒ the finite-population CLT and bounded-SNR conditions of RSI Prop. 6.1 hold) and **fails at criticality** — the physical realisation of the paper's own stated failure mode (non-CLT / strong dependence). This is also exactly where classical MCMC convergence diagnosis is hardest, which is what makes it a compelling applied finding rather than a mere caveat.

---

## 2. Synthesis — one principle, a measurable constant, a sharp boundary

Across all three systems the picture is coherent:

> From a uniform random `m`-subset you can estimate the mixing level and certify/reject equilibrium in an observation budget that is **independent of system size**, provided the observable has **bounded per-component SNR and short-range (finite-susceptibility) dependence**. The budget is `m*(δ) ≈ (z_α+z_β)² / (ρ² δ²)`, and `ρ²` — the single-component signal-to-noise ratio — is the one system-specific number that sets the cost (`2` for cards, `1` for the product chain). The principle **breaks at criticality**, where susceptibility diverges and subsampling hits an irreducible, size-dependent floor.

This is genuinely "mostly applied, partly theoretical": a working, benchmarkable diagnostic whose applicability is delimited by a clean, physically meaningful condition.

**Benchmark framing (vs Gelman–Rubin / R̂).** Classical diagnostics evaluate a scalar over the *full* state across multiple chains — `Θ(n)` reads per evaluation. The partial-observation certificate needs one chain and `Θ(m*) = Θ(1/δ²)` reads, independent of `n`. Both share the critical-slowing breakdown (neither escapes diverging susceptibility), so the honest claim is *cheaper, size-independent certification in the short-correlation regime* — not a universal replacement. A head-to-head on a real sampler is a Phase-5 deliverable.

---

## 3. What the pilots imply for the Phase-2 simulator (the point of sequencing this first)

The pilots surfaced concrete, non-obvious requirements. The simulator design from Phase 1 stands, with these refinements:

1. **`System` must expose an exact `true_phi(t)` and `rho2(phi)` where analytic, and a `susceptibility(state)` hook.** The card and chain systems have closed-form `ϕ(t)` and `ρ²(ϕ)`; Ising does not, and needs susceptibility/correlation-length measurement to *predict its own breakdown*. The interface must accommodate both "exact ground truth" and "measured ground truth (long-run proxy)" systems.
2. **Sampling must be memory-aware.** The naïve `argsort(random(R,n))` subset draw OOM-ed at `n≈10⁵` (a 4 GB allocation). The instrument needs an O(Rm) reservoir/partial-selection sampler, and should exploit coordinate exchangeability (first-`m`-columns) where the system permits — this alone unlocked `n` up to 65 536.
3. **Two ground-truth modes, cleanly separated:** *exact-marginal generation* (fast, faithful for product-form systems — let us reach 65 536 in seconds) vs *full-dynamics generation* (required whenever coordinates are coupled — Ising). The scorer must know which was used.
4. **An equilibration/autocorrelation guard.** Near Tc, fixed burn-in is untrustworthy (my Ising near-Tc numbers are approximate for this reason). The simulator must measure integrated autocorrelation time and either adapt burn-in or flag "not equilibrated," so ground truth is never silently wrong.
5. **A detection/power module is first-class, not an afterthought.** The headline metric is `m*(δ)` at a target power — the instrument must sweep `m`, estimate power/Type-I from replicate pools, and locate `m*` with a bracketing search (the coarse grid cost me precision in the first cards run).
6. **A dependence/susceptibility diagnostic per configuration** so every run reports whether it is inside or outside the applicability regime — turning "where it breaks" into an automatic label rather than a manual post-hoc.
7. **Seed discipline** (`SeedSequence.spawn` per `(config, replicate)`) — the pilots used ad-hoc seeds; the instrument needs reproducible, independent streams for defensible ISEF methodology.
8. **The observable and system layers must be independent** (confirmed: the same estimator/detector code ran unchanged across cards, chain, and Ising by swapping only the per-component readout).

---

## 4. Honest status

**Confirmed numerically (pilot-level, not yet proven):** the off-equilibrium SNR laws `ρ²(ϕ)` for cards and chain; n-independence of the precision and detection laws; the `m*≈c/δ²` certification cost with `c` set by `ρ²`; the criticality breakdown.

**Caveats:** Ising near-Tc figures use finite burn-in and are qualitatively (not quantitatively) reliable — the *floor-dominates-and-doesn't-fall-with-m* pattern is robust; the exact floor values are approximate. All pilots are exploratory (`R = 3000–6000`), adequate for pattern discovery, not for final publication error bars.

**To prove in Phase 6:** `ρ²(ϕ)` for cards with a controlled error term; the leading `O(1/n)` covariance correction; and a precise statement of the applicability boundary in terms of susceptibility / correlation length (a general condition: *scale-free certification ⇔ susceptibility bounded as `n→∞`*).

**Updated hypothesis to carry into the full experiment:** the certificate's cost is `m*(δ) = (z_α+z_β)²/(ρ² δ²)` with `ρ²` the equilibrium per-component SNR, valid iff the observable's susceptibility stays bounded; at a continuous phase transition the cost acquires an irreducible, size-dependent floor set by the critical susceptibility exponent.

---

*Figure `phase3_summary.png` accompanies this report: (A) off-equilibrium SNR laws, (B) n-independent precision, (C) size-independent detection cost `m*≈c/δ²`, (D) the Ising high-T success / near-Tc breakdown.*
