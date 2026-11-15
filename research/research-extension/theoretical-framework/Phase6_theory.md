# Phase 6 — Theory: Proofs and Derivations

*The mathematical track of the programme, developed in parallel with the Phase-2/5 compute campaign (which runs separately). Every result is tagged: **[PROVEN]** rigorous here; **[PROVEN — leading order]** rigorous to leading order in `n` with the correction's order controlled; **[cited]** rests on a classical result named in place; **[CONJECTURE / PCV]** not proven, awaiting Phase-5 evidence. Numerical spot-checks referenced below were run during Phase 6 and agree with every closed form claimed. Notation follows `Phase4_conjecture_and_design.md`.*

Throughout: finite Markov chain, stationary law `ν`; `n` components; per-component contribution `D_i(x)∈[0,B]`; intensive observable `O(x)=(1/n)Σ_i D_i(x)`, equilibrium mean `O_∞=E_ν[O]`; uniform `m`-subset `S` without replacement; estimator `Ô_m=(1/m)Σ_{i∈S}D_i`; within-configuration population variance `V_x=(1/n)Σ_i(D_i(x)−O(x))²`; observable susceptibility `χ_n[O]=n·Var_ν(O)`; relative floor `R_n=Var_ν(O)/O_∞²=χ_n[O]/(n O_∞²)`.

---

## Part I. The variance decomposition and the floor theorem (T3)

This part needs **no** structural assumption beyond finite second moments — no exchangeability, no CLT. It is the secure core of the programme.

### Lemma 1 (Exact variance decomposition). [PROVEN]

For a uniform random `m`-subset `S⊆[n]` drawn without replacement, independent of the chain, the estimator satisfies

```
Var(Ô_m) = (n−m)/(m(n−1)) · E_ν[V_x]  +  Var_ν(O).            (★)
```

**Proof.** Condition on the configuration `x`. Given `x`, `Ô_m` is the mean of a simple random sample without replacement from the finite population `{D_1(x),…,D_n(x)}`. Classical SRSWOR theory (Cochran 1977, Thm 2.2) gives the two conditional moments
`E_S[Ô_m∣x]=O(x)` and `Var_S(Ô_m∣x)=V_x/m·(n−m)/(n−1)`.
By the law of total variance,
`Var(Ô_m)=E_ν[Var_S(Ô_m∣x)]+Var_ν(E_S[Ô_m∣x]) = (n−m)/(m(n−1))·E_ν[V_x] + Var_ν(O).` ∎

The two terms are the **sampling term** (which one can drive down by observing more components) and the **floor** `Var_ν(O)` (which one cannot). Note `Var_ν(O)=χ_n[O]/n` by definition of `χ_n`.

### Theorem T3 (Floor lower bound and the extensive-susceptibility obstruction). [PROVEN]

For every `1≤m≤n`,
`Var(Ô_m) ≥ Var_ν(O) = χ_n[O]/n`, hence `CV²(Ô_m) ≥ R_n`,
with equality approached as `m→∞` (formally, as the FPC-weighted sampling term vanishes; exactly at `m=n`, where `Ô_n=O` and `Var(Ô_n)=Var_ν(O)`). Consequently:

- If `χ_n[O]=O(1)` (summable correlations, `a=0`), the floor is `R_n=Θ(1/n)→0`: certification is not floor-limited.
- If `χ_n[O]=Θ(n)` (extensive correlations, `a=1`), then `R_n=Θ(1)` is **bounded below by a positive constant for all `n`**, so there exists a tolerance `δ_0>0` such that **no sample of any size `m≤n`, at any `n`, can certify the observable's mean to relative precision `δ_0`.**

**Proof.** Both terms in (★) are non-negative and the second equals `Var_ν(O)`; drop the first for the inequality. Divide by `O_∞²` for the `CV²` form. The dichotomy substitutes `χ_n[O]=n·Var_ν(O)` into `R_n=Var_ν(O)/O_∞²=χ_n[O]/(nO_∞²)` with `O_∞=Θ(1)`. For `a=1`, `R_n→c>0`; any `δ_0<κ√c` is uncertifiable since even `m=n` gives `CV=√R_n≥√c`. ∎

**Remark (why this matters).** T3 converts the informal claim "condition (ii) is necessary" into a *proven* statement in the extensive limit: the partial-observation certificate has an information floor equal to the whole-system fluctuation, and extensive susceptibility makes size-free certification impossible. This holds for the Curie–Weiss counterexample (CX1) by construction, independent of any Phase-5 run. The sub-extensive rate `0<a<1` is **not** covered by T3 and remains a scaling conjecture (Part VI).

---

## Part II. The random-transposition anchor (T1)

Here `D_i=|π(i)−i|`, `O=AD`, `O_∞=AD_∞=(n²−1)/(3n)`.

### Proposition 2.1 (Single-tagged-component kernel). [PROVEN]

Under the random transposition shuffle, the position of one tagged card is a Markov chain on `[n]` with

```
P_t(i→j) = 1/n + (δ_{ij} − 1/n)·r_t,     r_t = (1 − 2/n)^t.
```

**Proof.** One step draws `(a,b)` uniform and independent and swaps positions `a,b`; a card at position `p` moves to `b` if `p=a`, to `a` if `p=b`, else stays. Hence the one-step matrix `M` has, for `q≠p`, `M_{pq}=P(\{a=p,b=q\}∪\{a=q,b=p\})=2/n²`, and `M_{pp}=1−2(n−1)/n²`. Write `M=(1/n)J+λ(I−(1/n)J)` with `J` the all-ones matrix: the uniform vector `1` is fixed with eigenvalue `1` (row sums are `1`), and any `v⊥1` gives `Mv=(M_{pp}−2/n²)v=(1−2n/n²)v=(1−2/n)v`, so `λ=1−2/n` with multiplicity `n−1`. Thus `M^t=(1/n)J+λ^t(I−(1/n)J)`, which is the stated kernel. ∎ *(Numerically confirmed: e.g. `n=16,t=10`, `P_home` 0.3089 vs 0.3091.)*

Equivalently: **with probability `r_t` the card sits exactly at its origin; with probability `1−r_t` its position is uniform on `[n]`.**

### Proposition 2.2 (Exact mean trajectory). [PROVEN]

`E[AD(t)] = (1 − r_t)·AD_∞`, i.e. `ϕ(t)=1−(1−2/n)^t`. In normalised time `t=θ·t_mix`, `t_mix=½ n ln n`, `ϕ(θ)=1−n^{−θ}(1+o(1))`.

**Proof.** By linearity (which needs only single-card marginals, so inter-card correlations are irrelevant to the mean),
`E|π_t(i)−i|=Σ_j|j−i|P_t(i→j)=Σ_j|j−i|[1/n+r_t(δ_{ij}−1/n)] = (1−r_t)·(1/n)Σ_j|j−i|`,
using `|i−i|=0`. Averaging over `i` gives `E[AD(t)]=(1−r_t)·(1/n)Σ_i(1/n)Σ_j|j−i|=(1−r_t)AD_∞`. The `θ` form uses `(1−2/n)^{θ·½n\ln n}=\exp(−θ\ln n(1+o(1)))=n^{−θ}(1+o(1))`. ∎

**Corollary 2.3 (Derivation of the RSI Figure-4 residual). [PROVEN]** The normalised curve `AD(t)/AD_∞=1−n^{−θ}` is *not* exactly size-independent; at fixed `θ>0` it increases toward 1 with `n`, and the family collapses only as `n→∞`. This reproduces, analytically, the empirically-observed ordering in the RSI paper (larger decks relax slightly faster in normalised time, the spread narrowing with `n`), which the paper reported without derivation.

### Proposition 2.4 (Equilibrium second moments — fully exact). [PROVEN]

At equilibrium (`π` uniform), writing `F=Σ_i|π(i)−i|` for Spearman's footrule (`AD=F/n`):
(i) `E_π[(1/n)Σ_i D_i²]=(n²−1)/6`.
(ii) **`Var(F)=(n+1)(2n²+7)/45`**, hence `Var_π(AD)=Var(F)/n²=(n+1)(2n²+7)/(45n²)=Θ(n)` with leading coefficient `2/45`.
(iii) `E_π[V_π]=(n²−1)(n²+2)/(18n²)−Var_π(AD)=(n−2)(n+1)(5n²+n+12)/(90n²)`.

**Proof.** (i) For uniform `π`, `π(i)` is uniform on `[n]`, so `E[D_i²]=(1/n)Σ_j(j−i)²`; averaging, `E[(1/n)ΣD_i²]=(1/n²)Σ_{i,j}(j−i)²=E[(X−Y)²]=(n²−1)/6` for `X,Y` i.i.d. uniform. (ii) Write `F=Σ_iD_i`. Then `E[F²]=Σ_iE[D_i²]+Σ_{i≠j}E[D_iD_j]`. The diagonal is `Σ_iE[D_i²]=n(n²−1)/6`. For the off-diagonal, `(π(i),π(j))` is uniform over ordered pairs of distinct values, so `E[D_iD_j]=\frac{1}{n(n−1)}Σ_{k≠l}|k−i||l−j|`. Summing over `i≠j` and removing the two distinctness constraints by inclusion–exclusion gives `Σ_{i≠j}Σ_{k≠l}|k−i||l−j|=P²−2Q+R_2`, where `P=Σ_{i,k}|i−k|=n(n²−1)/3`, `R_2=Σ_{i,k}(i−k)²=n²(n²−1)/6`, and `Q=Σ_i g(i)²` with `g(i)=Σ_k|i−k|=i²−(n+1)i+n(n+1)/2`. Evaluating `Q=n(n−1)(n+1)(7n²−8)/60` and assembling `E[F²]−(E[F])²` with `E[F]=(n²−1)/3` yields `Var(F)=(n+1)(2n²+7)/45`. (iii) `E_π[V_π]=(n²−1)/6−AD_∞²−Var_π(AD)` and `(n²−1)/6−AD_∞²=(n²−1)(n²+2)/(18n²)`; subtracting (ii) gives the stated closed form. ∎ *(Enumeration check `n=3`: `Var(F)=100/45`, matching the 6 permutations of `S_3` exactly. Symbolic `Q`, `Var(F)`, `E_π[V_π]` computed and cross-checked against Monte Carlo.)*

### Theorem T1 (Anchor: exact size-uniform relative variance). [PROVEN — exact closed form]

At equilibrium, for a uniform `m`-subset, combining Lemma 1 with Prop. 2.4 gives the **exact** law

```
Var(Ô_m) = (n−m)/(m(n−1))·E_π[V_π] + Var_π(AD),
CV²(ϕ̂_m) = Var(Ô_m)/AD_∞²
         = −(m n³ −5m n² −4m n −10m −5n⁴ +9n³ −10n² +24n) / (10 m (n−1)³ (n+1)).
```

It decomposes exactly into sampling term plus floor,

```
CV²(ϕ̂_m) = (n−m)(n−2)(5n²+n+12) / (10 m (n−1)³ (n+1))  +  R_n,
R_n = (2n²+7) / (5(n−1)²(n+1)),        n·R_n → 2/5.
```

Hence for `m≪n`, `CV²(ϕ̂_m)=1/(2m)·(1+O(1/n)+O(m/n)) + R_n` with **`R_n=(2/5)·n^{−1}(1+O(1/n))`**, so `CV(ϕ̂_m)=1/√(2m)(1+o(1))` is **independent of `n` to leading order**, recovering RSI Theorem 5.3 as the `n→∞` limit and exhibiting the T3 floor with its **exact constant `c=2/5`**.

**Proof.** Substitute the exact `E_π[V_π]=(n−2)(n+1)(5n²+n+12)/(90n²)` and `Var_π(AD)=(n+1)(2n²+7)/(45n²)` (Prop. 2.4) and `AD_∞²=(n²−1)²/(9n²)` into Lemma 1 and simplify; the floor `R_n=Var_π(AD)/AD_∞²=9\,Var(F)/(n²−1)²=(2n²+7)/(5(n−1)²(n+1))`, whose leading term is `2/(5n)`. ∎ *(The exact `CV²` was checked against Monte Carlo at `(n,m)=(256,32),(512,64),(2048,128)`: ratios `1.004, 1.000, 1.001`. The constant `c=2/5=0.40` matches the Phase-5 measurement `n·R_n=0.37–0.41`.)*

### Corollary T1′ (Exact certification cost with floor). [PROVEN — modulo the Gaussian test]

For estimation to relative tolerance `ε` at confidence `1−α` (`κ_est=z_{1−α/2}`) or detection of a deviation `δ` (`κ=z_{1−α}+z_{1−β}`), the leading floor-corrected sample is
`m*(ε) = \dfrac{5κ²n}{2(5ε²n − 2κ²)} = \dfrac{κ²}{2ε²}\big(1+O(1/(ε²n))\big)`,
which **exceeds** the size-free value `κ²/(2ε²)` by a floor penalty and **diverges to infeasibility as `ε↓κ√(2/(5n))=κ√R_n(1+o(1))`** — the exact instance of the T3 obstruction. For `ε≫√R_n` (the large-system regime) the cost is `κ²/(2ε²)`, independent of `n`.

---

## Part III. Off-equilibrium signal-to-noise (leading order)

### Proposition 3.1 (Mixture form of the off-equilibrium single-card law). [PROVEN]

At time `t`, writing `r=r_t`, `ϕ=1−r`: a tagged card's contribution is `D_i(t)=0` with probability `r`, and `D_i(t)=|U−i|` with `U` uniform on `[n]` with probability `1−r`. Hence `E[D_i(t)]=ϕ\,μ_i` and `E[D_i(t)²]=ϕ\,ν_i`, where `μ_i=(1/n)Σ_j|j−i|`, `ν_i=(1/n)Σ_j(j−i)²`.

**Proof.** Immediate from Prop. 2.1 in the form `P_t(i→j)=r\,δ_{ij}+(1−r)(1/n)`. ∎

### Proposition 3.2 (Trajectory SNR and precision). [PROVEN — leading order]

Define the sampling SNR `ρ²(ϕ)=(E[O])²/E_x[V_x]` at mixing level `ϕ`. Then, to leading order in `n`,

```
ρ²(ϕ) = ϕ / (3/2 − ϕ),        Var(ϕ̂_m) = ϕ(3−2ϕ)/(2m)·(1+O(1/n)),
```

recovering `ρ²(1)=2` and `Var(ϕ̂)=1/(2m)` at equilibrium.

**Proof.** With `W` the contribution of a uniformly-chosen card, Prop. 3.1 gives `E[W]=ϕ\,AD_∞` and `E[W²]=ϕ·(n²−1)/6`. To leading order `E_x[V_x]=E[W²]−(E[O])²+o(n²)=ϕ n²/6−ϕ²n²/9`, so
`ρ²(ϕ)=(ϕ²n²/9)/(ϕn²/6−ϕ²n²/9)=ϕ/(3/2−ϕ)`. The precision follows from `Var(ϕ̂_m)=E_x[V_x]/(m\,AD_∞²)·(1+o(1))=[ϕ/6−ϕ²/9]·9/m=ϕ(3−2ϕ)/(2m)`. ∎ *(Numerically confirmed across `ϕ∈{0.26,0.55,0.87}`: e.g. `ρ²` 1.363 vs 1.363.)*

The dropped term is the permutation-level fluctuation `Var_x(AD(t))`, of order `n` (Part II), giving the `O(1/n)` correction — the same floor as at equilibrium.

---

## Part IV. Detection corollary

### Corollary 4.1 (Size-free certification cost). [PROVEN — under the finite-population CLT]

Assume the finite-population CLT (Hájek 1960) applies to `Ô_m` at equilibrium, so the test statistic is asymptotically standard normal (verified in RSI/Phase-3 via skewness/kurtosis at `m≥32`). Testing `H0: ϕ=1` against `H1: |1−ϕ|≥δ` at size `α`, power `1−β`, the minimum sample is
`m*_det(δ)=κ²/(2δ²)·(1+o(1))`, `κ=z_{1−α}+z_{1−β}`, **independent of `n`**.

**Proof.** Under `H1`, the standardised gap is `δ\,AD_∞/\mathrm{SD}(Ô_m)`; by T1, `\mathrm{SD}(Ô_m)=AD_∞·1/√(2m)(1+o(1))`, so the gap is `δ√(2m)(1+o(1))`. The design equation `δ√(2m)=κ` yields `m*=κ²/(2δ²)`; `n` cancels because both signal and noise carry the factor `AD_∞`. ∎

*Caveat (flagged in the critique).* This is the idealised-null cost. At small `m` the self-normalised statistic is skewed and the Gaussian size is mis-calibrated; the deployable procedure needs an exact/bootstrap null there, and the constant should be read as `(1+o(1))`. This is a Phase-5 measurement, not a theorem.

---

## Part V. General sufficiency (T2)

### Theorem T2 (Scale-free certification, sufficiency). [PROVEN — under stated hypotheses; Hájek cited]

Let `{D_i^{(n)}}_{i=1}^n` be a triangular array of uniformly bounded contributions of an equilibrated chain with intensive mean `O_∞^{(n)}=Θ(1)`, satisfying

1. **(Non-degenerate SNR)** `ρ_n²:=(O_∞)²/E_ν[V_x]→ρ_∞²∈(0,∞)`;
2. **(Summable susceptibility, `a=0`)** `χ_n[O]=n\,\mathrm{Var}_ν(O)=O(1)`;
3. **(Lindeberg/negligibility)** the finite-population CLT of Hájek (1960) holds for `Ô_m` as `m,n−m→∞`.

Then `CV²(Ô_m)=1/(ρ_∞²m)+o(1/m)+O(1/n)`, and the certification cost is
`m*(δ)→κ²/(ρ_∞²δ²)`, independent of `n`.

**Proof.** By Lemma 1, `CV²(Ô_m)=(n−m)/(m(n−1))·E_ν[V_x]/O_∞² + \mathrm{Var}_ν(O)/O_∞²`. Hypothesis 1 gives `E_ν[V_x]/O_∞²=1/ρ_n²→1/ρ_∞²`; hypothesis 2 gives the floor `\mathrm{Var}_ν(O)/O_∞²=χ_n/(nO_∞²)=O(1/n)`. For `m≪n` the FPC is `1+O(m/n)`, so `CV²=1/(ρ_∞²m)(1+o(1))+O(1/n)`. Hypothesis 3 makes the standardised detection statistic asymptotically normal, so the design equation `δ/CV=κ` inverts to `m*→κ²/(ρ_∞²δ²)`. ∎

T2 is exactly the `a=0` clause of the refined conjecture, now a theorem under checkable hypotheses. The random-transposition shuffle (`ρ_∞²=2`) and the product chain (`ρ_∞²=1`, where 1–3 hold *exactly*) are instances.

---

## Part VI. What is proven, and what remains open

**Proven now (independent of any Phase-5 run):**
- The exact variance decomposition (Lemma 1) and the **floor lower bound** (T3), including the **`a=1` impossibility** — the rigorous necessity half of condition (ii).
- The **single-card kernel** and the **exact mean law** (Props 2.1–2.2), and with them the analytic explanation of the RSI Figure-4 residual (Cor 2.3).
- The **exact equilibrium moments** (Prop 2.4: `Var(F)=(n+1)(2n²+7)/45`, enumeration-verified) and the **fully exact anchor law** T1, including the closed-form `CV²`, the exact floor `R_n=(2n²+7)/(5(n−1)²(n+1))` with **constant `c=2/5`**, and the exact certification cost T1′.
- The **off-equilibrium SNR** `ρ²(ϕ)=ϕ/(3/2−ϕ)` and precision law (Prop 3.2), leading order.
- The **general sufficiency theorem** T2 under summable-susceptibility + CLT hypotheses.

**Open / conjectural (Phase-5 targets, PCV):**
- The **sub-extensive rate** `0<a<1` (critical magnetization): T3 gives a floor but not its `n`-exponent; the claim `R_n≍n^{a−1}` with `a` from the universality class is a scaling **conjecture**, to be measured (P2).
- The **necessity of the summable condition in the borderline `a<1` regime** (T3 proves it only at `a=1`).
- The **exact small-`m` (non-CLT) law** and the true finite-`m` detection constant (Cor 4.1 caveat).
- The **exact off-equilibrium correction constant** in `ρ²(ϕ)` (the `O(1/n)` sub-leading term, needing the off-equilibrium footrule-type variance; the leading `ρ²(ϕ)` and equilibrium constants are exact).

**How Phase 5 checks the theory.** When the compute campaign returns: the `a=0` systems must reproduce `m*δ²→κ²/ρ_∞²` and a flat floor slope `−1` (confirming T1/T2); the mean-field system must show the T3 hard floor (`δ_floor→const`); the critical-Ising floor slope tests the *conjectured* `a=7/8` rate that no theorem yet fixes. A clean disagreement on the last would revise the conjecture, not the theorems — the theorems above stand on the proofs, not the simulator.

---

### Honesty note

Nothing here is claimed beyond its tag. T3, Props 2.1–2.2, Prop 2.4 and **T1 (now an exact closed form, MC- and enumeration-verified)** are unconditional. Prop 3.2 is leading-order with the correction's order controlled and numerically confirmed. Cor 4.1, T1′ and T2 are conditional on the Hájek finite-population CLT, standard for bounded populations with non-degenerate variance and explicitly checked empirically at `m≥32`. The critical-exponent rate is the one genuinely new physical claim that remains a conjecture pending measurement.
