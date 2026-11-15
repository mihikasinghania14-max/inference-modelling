# Phase 5 — Ising Criticality: the last conjecture, measured

*Run in the cloud container with a Wolff cluster sampler at `T_c = 2/ln(1+√2)`, sizes `L ∈ {16,24,32,48,64,80}` (`n` up to 6400), `R = 1600–2500` decorrelated equilibrium samples per size (decorrelation flips scaled with `L`; `⟨E⟩` converges toward the Onsager value, variances monotone — equilibration is clean). This resolves the one claim that was still **[CONJECTURE / PENDING]** in Phases 4–7.*

## Result

| Quantity | Measured slope (log–log vs `n`) | Predicted | Verdict |
|---|---|---|---|
| `χ[M] = n·Var(|M|)` | **+0.858** | +7/8 = +0.875 | **CONFIRMED** `a[M] ≈ 7/8` |
| `Var(|M|)` | −0.142 | −1/8 = −0.125 | confirmed |
| **`R_n[M] = Var(|M|)/⟨|M|⟩²`** | **−0.023 (flat, ≈0.066)** | — | **relative floor is `Θ(1)`** |
| `Var(E)` (energy density) | −0.867 | ≈ −1 (+log) | confirmed |
| **`R_n[E] = Var(E)/⟨E⟩²`** | **−0.855** (0.0135 → 0.0009) | — | **floor vanishes** |

## What it means (and how it sharpens the conjecture)

1. **The magnetization susceptibility exponent `a[M] = 7/8` is confirmed.** The critical order parameter is strongly, sub-extensively correlated, exactly as the universality class predicts.

2. **But the operative quantity is the *relative* floor `R_n`, not the raw exponent `a` — and for the order parameter `R_n[M] = Θ(1)`.** The Phase-4 framework wrote `R_n ≍ n^{a−1}` assuming `O_∞ = Θ(1)`. At criticality the magnetization's own mean collapses (`⟨|M|⟩ ~ n^{−1/16}`), so both numerator and denominator vanish at the same rate and the relative floor is **constant** (measured ≈ 0.066, flat across a 25× range in `n`). **Practical consequence: certifying the magnetization at criticality is impossible below a fixed tolerance (~√0.066 ≈ 26%), for any sample size at any `n`** — the same fate as the extensive `a=1` case, reached by a different route (a vanishing signal rather than an extensive variance).

3. **The same system certifies fine in a different observable.** Energy density has `O_∞ = Θ(1)` and a floor that **vanishes** as `n^{−0.86}`, so it is asymptotically size-free at the very same critical point. This is the cleanest possible confirmation of the central refined thesis: **size-free certification is a property of the observable's relative floor, not of the system** — one system, one critical point, opposite outcomes for magnetization vs energy.

## Status changes

- Phase-7 §6.3 (Ising critical exponent): **PENDING → READY.** Figure 6.2 (floor slope vs observable at `T_c`) can be drawn from these numbers.
- The refined statement for the paper: at a continuous transition the *order-parameter* certificate acquires a **non-vanishing relative floor** (impossible, like `a=1`), while *non-critical extensive observables* (energy density) stay certifiable. The operative criterion is `R_n → 0`, i.e. `χ_n[O] = o(n·O_∞²)`, which correctly reduces to "summable susceptibility" only when `O_∞ = Θ(1)`.
- This also cross-checks **Theorem T3**: a non-vanishing relative floor forbids size-free certification, and here it is measured directly for the critical order parameter.

## Caveat

Finite-size study to `L=80`; slopes carry the usual few-percent fit uncertainty and mild corrections-to-scaling (the `χ[M]` slope drifts toward 7/8 as small-`L` points are dropped). A production run for the paper should extend to `L≈128–256` with bootstrap-CI slope fits and an explicit corrections-to-scaling term, but the qualitative result — `R_n[M]=Θ(1)`, `R_n[E]→0` — is robust and already decisive.
