# Phase 5 — Architecture Proposal for the Partial-Observation Research Simulator

*Deliverable produced before any code is written, per the build specification §"Your explicit first action". Maps the frozen Phase-4 scientific specification onto a concrete, modular, reproducible, memory-efficient Python instrument. Nothing here reopens the science: the conjecture, predictions (P1–P8), experimental matrix (§4), controls (§5), and counterexamples (§6) of `Phase4_conjecture_and_design.md` are treated as frozen. This document is the code design + test plan + experiment map, plus a short list of modelling subtleties I found while cross-reading the four authoritative documents.*

Status legend used throughout, mirroring the epistemic guardrails: **[EST]** established mathematically · **[EMP]** established empirically (Phase-3) · **[CONJ]** conjectured · **[PCV]** pending computational validation. The simulator **measures**; comparison to any [EST]/[EMP]/[CONJ] value happens **only** in the analysis/scoring layer.

---

## 0. Design invariants (the non-negotiables, encoded structurally)

These are enforced by architecture, not convention:

1. **The veil is physical.** The inference engine receives *only* the `m` revealed contributions — never the full state. Enforced by an `Observation` dataclass that carries exactly `m` values plus `(n, m, seed, config_hash)`; the full state is not a field. An **access-audit** counter asserts exactly `m` component reads per replicate. (Handoff §3, §Controls; Phase-4 §5.5.)
2. **No prediction is ever hard-coded into an estimator.** Estimators compute raw moments (`Ô_m`, `V_x`, `Var_ν(O)`, skew/kurtosis…). Predicted values (`ρ²(ϕ)=ϕ/(3/2−ϕ)`, `κ²/ρ²`, `a`, floors) live only in `analysis/predictions.py` and are used only for scoring. Grep-able separation: nothing in `systems/`, `observables/`, `inference/`, `measure/` imports from `analysis/predictions.py`.
3. **Ground truth ≠ dynamics.** Every `System` exposes `sample_equilibrium()` (a *correct rapid* sampler) as a method distinct from `step()` (the *physical dynamics* whose convergence is certified). For Ising these are literally different algorithms (Wolff vs single-site Glauber). The scorer records which path produced each ensemble.
4. **Relative/dimensionless quantities are canonical.** All headline quantities are dimensionless (`ϕ`, `ρ²`, `R_n`, `CV²`, `I_eff`). The code never assumes `D_i ≤ B` for a constant `B`, nor `O_∞ = Θ(1)` — see modelling note M1. This lets cards (`O_∞=Θ(n)`) and the fixed-point fraction (`O_∞→0`) flow through the same pipeline.
5. **Reproducibility is structural.** Every stored row carries its 128-bit seed + config hash; every run writes a manifest (root seed, code version, library versions, host, timestamp). A row can be regenerated from `(config, seed, code version)` alone.
6. **Measure, then score, then plot.** Figures are generated from stored tables, never from live runs.

---

## 1. Package layout

```
rsi_partial/                          # pip-installable (pyproject.toml, NumPy/SciPy/pandas/pyarrow/matplotlib)
├── core/
│   ├── interfaces.py      # System, Observable, GroundTruth enum, RegimeTag ABCs
│   ├── rng.py             # SeedSequence.spawn protocol, per-config child streams, config hashing
│   └── types.py           # Observation, ConfigKey, Record dataclasses (frozen); the "veil" boundary
├── systems/
│   ├── product_chain.py   # S1  (exact marginal; product Glauber)
│   ├── cards.py           # S2  (direct shuffle; random-transposition dynamics)
│   ├── ising2d.py         # S3  (Wolff/SW cluster GT; single-site Glauber dynamics)
│   ├── curie_weiss.py     # S4  (Stage 3)
│   ├── ising1d.py         # S5  (transfer-matrix exact GT; Glauber) (Stage 4)
│   └── ssep.py            # S6  (product-Bernoulli GT; conserved exchange) (Stage 4)
├── observables/
│   ├── cards_ad.py            # AD = (1/n)Σ|π(i)-i|
│   ├── cards_fixedpoint.py    # S2' F = (1/n)Σ 1{π(i)=i}  (Stage 0 — degenerate-SNR interface proof)
│   ├── chain_mag.py           # magnetization
│   └── ising_obs.py           # |M|, energy density, block-mean
├── observe/
│   └── veil.py            # ObservationModel: memory-aware O(Rm) subsampler (no full-state leak)
├── inference/
│   ├── estimator.py      # Ô_m, sampling variance (with/without FPC), variance split
│   ├── bayes.py          # conjugate Gaussian posterior (RSI §5 / App F), credible intervals
│   ├── detect.py         # H0:ϕ=1 vs H1:|1-ϕ|≥δ; power & Type-I from replicate pools
│   ├── mstar.py          # m*(δ) finder: estimation-m* AND detection-m*, bisection to ±5%
│   └── selfcalib.py      # ρ̂², R̂_n from the observed subsample → m̂*   (Stage 4)
├── measure/
│   └── diagnostics.py    # per-config metric bundle (see §4); τ_int, R_eff, equilibration checks
├── storage/
│   ├── schema.py         # canonical Record schema + validation
│   └── io.py             # Parquet (CSV fallback), deterministic filenames, JSON manifests
├── analysis/
│   ├── predictions.py    # ALL theory values (P1–P8 support/falsify bands). Isolated.
│   ├── bootstrap.py      # bootstrap 95% CIs (≥1000 resamples) on every headline quantity
│   ├── score.py          # per-prediction support/falsify/inconclusive verdict tables
│   └── figures.py        # figure factory (from stored tables only)
├── experiments/
│   ├── registry.py       # declarative config grids (Phase-4 §4.2), staged
│   ├── stage0_validation.py  stage1_afree.py  stage2_critical.py  stage3_factorial.py  stage4_secondary.py
│   └── run.py            # CLI: build configs → run → store; resumable; smoke (reduced R) vs full
└── tests/
    ├── test_analytical.py    # asserts against closed forms only (no self-referential fixtures)
    ├── test_veil.py          # access audit, relabelling invariance, O(Rm) memory
    ├── test_identity.py      # variance identity = sampling term + floor
    ├── test_ising_gt.py      # Wolff vs Onsager energy/|M|; Glauber detailed balance
    └── test_regression_phase3.py  # Phase-3 numbers within MC error (the gate)
```

---

## 2. Core interfaces

### 2.1 `System` (ABC)

```python
class System(ABC):
    name: str
    ground_truth: GroundTruth            # EXACT_MARGINAL | MCMC_CLUSTER | TRANSFER_MATRIX
    def sample_equilibrium(self, n, rng, *, params) -> State: ...   # correct RAPID sampler
    def step(self, state, rng, *, params) -> State: ...            # PHYSICAL dynamics (certified)
    def initial_state(self, n, rng, *, kind) -> State: ...         # 'ordered' | 'random' (equilibration check)
    def true_phi(self, t, n, *, params) -> float | None: ...       # analytic where known, else None
    # metadata (None where unknown — never faked):
    def O_infinity(self, n, *, params) -> float | None: ...
    def rho2(self, phi, *, params) -> float | None: ...            # analytic SNR function where known
    def susceptibility_exponent_a(self, *, params) -> float | None: ...
    def t_mix(self, n, *, params) -> int | None: ...
```

`State` is a plain NumPy array (permutation for cards; `{0,1}^n` for chain; `±1` `L×L` for Ising) — the system owns its encoding. `params` carries e.g. temperature `T`, target mixing level.

- **S1 product chain:** state `∈{0,1}^n`. `sample_equilibrium`: each coord iid `Bernoulli(ϕ/2)` (exact). `step`: product Glauber (each coord independently relaxes toward 1/2). `true_phi(t)=` analytic (independent-coordinate relaxation). `O_infinity=1/2`, `rho2(ϕ)=ϕ/(2−ϕ)` **[EST]**, `a=0` **[EST]**.
- **S2 cards:** state = permutation array `π`. `sample_equilibrium`: `rng.permutation(n)` (exact uniform). `step`: draw `a,b ~ U{0..n-1}` independently (incl. self-loops, per Def 2.1 / App C.2), swap. `true_phi(t)=1−(1−2/n)^t` **[EST]**. `O_infinity=(n²−1)/(3n)` **[EST]** (AD observable). `rho2(1)=2` **[EST]**; `rho2(ϕ)=ϕ/(3/2−ϕ)` is **[EMP]/leading-order** → lives in predictions, **not** returned as truth (return the leading-order form only tagged, or `None`; see M1). `a=0` (via `R_n≍1/n`) **[EMP→CONJ]**.
- **S3 Ising 2D:** state `±1` `L×L`, periodic BC. `ground_truth=MCMC_CLUSTER`. `sample_equilibrium`: **Wolff** cluster (bypasses critical slowing); Swendsen–Wang available as cross-check. `step`: **single-site Glauber** (the certified physical dynamics). `T_c=2/ln(1+√2)`. Metadata mostly `None` (measured, not analytic) except energy/|M| comparison to Onsager in tests.

### 2.2 `Observable` (ABC)

```python
class Observable(ABC):
    name: str
    def per_component(self, state) -> np.ndarray:   # D_i[]  — length n, the ONLY thing the veil samples
    def population_value(self, state) -> float:     # O(x) = mean(D_i)  (used for GT + m=n check)
    def E_infinity(self, n, *, params) -> float | None: ...
    def rho2(self, phi, *, params) -> float | None: ...
    def a(self, *, params) -> float | None: ...
```

Observables are decoupled from systems: `cards→{AD, fixed_point_fraction}`, `chain→{magnetization}`, `ising→{abs_magnetization, energy_density, block_mean}`. The veil calls only `per_component`; `population_value` is used solely for ground-truth diagnostics and the `Ô_n=O` identity check, never fed to inference.

**Ising `per_component` definitions** (each returns a length-`n=L²` array; the veil samples `m` of them):
- `abs_magnetization`: `D_i = s_i` (signed spin); population value `|mean(s)|`. Note |M| is a function of the *global* sum, so the estimator target is the mean spin and |·| is applied at the population level — the per-component readout is `s_i`. (Sign-symmetry handled in analysis.)
- `energy_density`: `D_i = −s_i·(Σ_{j∼i} s_j)/2` (per-site bond energy, double-counting corrected); population value = mean.
- `block_mean`: `D_i` = mean spin of the block containing site `i` (coarse-grained observable for the conserved/block comparison).

### 2.3 The veil — `ObservationModel` (`observe/veil.py`)

```python
@dataclass(frozen=True)
class Observation:
    values: np.ndarray      # exactly m revealed D_i  (the ONLY payload)
    n: int; m: int; seed: int; config_hash: str
    # NO full-state field. Inference sees this object only.
```

Memory strategy for the `O(Rm)` hot path (M-strategy, resolves the Phase-3 OOM):

- **Naïve `argsort(random(R,n))` is banned** (`O(Rn)`, OOM at `n≈10⁵`).
- For a single replicate: draw an `m`-subset without replacement via **partial Fisher–Yates / `rng.choice(n, m, replace=False)`** → `O(m)` indices, then gather `D[idx]` → `O(m)`. Across `R` replicates the loop is `O(Rm)` in the sampling, never materializing an `R×n` array.
- **Exchangeability fast-path (valid only where mathematically justified):** when the observable's components are exchangeable under `ν` (product chain; cards `AD` at equilibrium — every card has common inclusion probability `m/n`, Lemma 2.10), the *first-`m`* components of a freshly generated equilibrium state are a valid uniform `m`-subset, so we can generate only what we read. This is the trick that unlocked `n=65 536` in Phase-3. **Guarded by a per-system/observable `exchangeable_under_nu` flag** and cross-checked once against explicit `choice`-based sampling in tests. Not used for Ising (non-exchangeable) or off-equilibrium sampling where position matters.
- **Component relabelling** (control §5.5): an optional random permutation applied between generation and observation, asserting invariance and precluding positional leakage.

### 2.4 Inference engine (`inference/`)

- **`estimator.py`** — `Ô_m = mean(values)`; sampling variance estimate both raw (`s²/m`) and finite-population-corrected (`s²/m · (n−m)/(n−1)`); and the **variance split**: given a replicate pool, decompose `Var(Ô_m)` into sampling term `(1/m)(n−m)/(n−1)·E_x[V_x]` and floor `Var_x(O)` (identity ★). Both terms measured independently → P1/P8.
- **`bayes.py`** — conjugate Gaussian posterior exactly per RSI §5 / App F: weak prior `µ0=O_∞`, configurable `σ0` (default sweepable `{50,75,100,150,200}` in raw scale, or scale-relative); closed-form `1/σ²_post=1/σ²_0+1/σ²_m`, `µ_post=σ²_post(µ0/σ²_0 + Ô_m/σ²_m)`; credible interval `µ_post ± z_{1−α/2}σ_post`. `σ²_m` is the **measured** sampling variance (or the analytic `n²/18m` for cards as a labelled cross-check, never as the estimator's own input).
- **`detect.py`** — `H0:ϕ=1` vs `H1:|1−ϕ|≥δ`, one-sided `α`. Power and Type-I estimated from replicate pools at each `(n,m,δ,ϕ)`. Both the Gaussian/CLT test and, for small `m≤8`, a bootstrap/exact null calibration (control §5.7).
- **`mstar.py`** — **two distinct m\*** (see modelling note M2), each by bisection to ±5% on a geometric `m` grid refined locally:
  - **detection-m\***: smallest `m` giving power ≥ `1−β` at size `α` for deviation `δ`. Target constant `κ²/(ρ²δ²)` **[EMP]**.
  - **estimation-m\***: smallest `m` with posterior credible half-width ≤ `ε·O_∞` at `1−α`. Target `z²_{1−α/2}/(ρ²ε²)` **[EST, RSI Eq 5.7]**.
- **`selfcalib.py`** (Stage 4) — `ρ̂²`, `R̂_n` from the *same* subsample; predict `m̂*`; store beside truth.

### 2.5 Measurement / diagnostics (`measure/diagnostics.py`)

Per configuration, computed and stored automatically (Phase-4 §8.6):
`O`, `Ô_m`, `ϕ̂`; `Var(Ô_m)` and its split (sampling term, floor); `Var_ν(O)`, `χ_n[O]`, `χ_n/n`, `R_n`; measured `ρ²(ϕ)`; pairwise-correlation / susceptibility proxy; sampling-distribution **skewness** and **excess kurtosis**; **τ_int** (batch-means/windowed ACF) and **R_eff**; **I_eff(m)=1/CV²** and crossover `m_×`; empirical `m*`. Self-calibrated `ρ̂²`, `R̂_n` stored beside ground-truth. Equilibration check (ordered-vs-random initial condition agreement; running-trace stationarity) emitted as pass/flag.

---

## 3. Data schema (per-`(config, replicate)` record)

Tidy columnar (Parquet; CSV fallback). One row per replicate; aggregates derived in analysis.

| field | type | meaning |
|---|---|---|
| `experiment_id` | str | e.g. `stage1_afree` |
| `config_hash` | str | 16-hex of canonicalized config (deterministic filename key) |
| `system` | str | `product_chain`/`cards`/`ising2d`/… |
| `observable` | str | `AD`/`fixed_point_fraction`/`abs_magnetization`/… |
| `ground_truth` | str | `exact_marginal`/`mcmc_cluster`/`transfer_matrix` |
| `n` | int | components (`=L²` for Ising) |
| `L`, `T`, `T_over_Tc` | int/float | Ising only (else null) |
| `regime` | str | `equilibrium`/`off_equilibrium` |
| `phi_target` | float | target mixing level (1.0 at equilibrium) |
| `theta` | float | `t/t_mix` where applicable (else null) |
| `m` | int | sample size (`m=n` census flagged) |
| `delta`, `epsilon` | float | detection deviation / estimation tolerance (null where n/a) |
| `replicate` | int | replicate index |
| `seed_128` | str | 128-bit child seed (hex) |
| `code_version` | str | git commit / version tag |
| `O_full` | float | population value `O(x)` (GT diagnostic) |
| `Ohat_m` | float | estimator `Ô_m` |
| `phi_hat` | float | `Ô_m / O_∞` |
| `Vx` | float | within-config population variance `V_x` |
| `samp_var_raw`, `samp_var_fpc` | float | sampling variance est. (with/without FPC) |
| `skew`, `exkurt` | float | sampling-distribution shape flags |
| `post_mean`, `post_sd`, `ci_lo`, `ci_hi` | float | Bayesian posterior + credible interval |
| `reject_H0` | bool | detection decision at this `(δ,α)` |
| `tau_int`, `R_eff` | float | autocorrelation diagnostics (locally-generated only) |
| `equilibrated_flag` | bool | passed ordered-vs-random + stationarity |
| `access_reads` | int | audited component-read count (must equal `m`) |
| `rho2_hat`, `Rn_hat` | float | self-calibrated (Stage 4; null earlier) |

**Aggregated summaries** (separate table, one row per `(config)`): `Var(Ô_m)` across replicates + its split, `Var_ν(O)`, `χ_n`, `R_n`, `ρ²(ϕ)`, `I_eff(m)`, `m*` (both kinds) — each with **bootstrap 95% CI**. **JSON manifest** sidecar per run: root seed, code version, library versions (numpy/scipy/pandas/pyarrow), host, UTC timestamp, full config grid.

---

## 4. Seed protocol (`core/rng.py`)

- One documented **root seed** per experiment (recorded in manifest).
- `SeedSequence(root).spawn` a fresh child per **`(system, regime, n, m, δ, ϕ, replicate)`** tuple, canonicalized and hashed; the child's 128-bit entropy is stored per row. Independent replicates → non-overlapping streams by construction (spawn guarantees).
- **Config hash** = BLAKE2b of the sorted canonical config dict → deterministic filenames; nothing overwrites silently (write-if-absent + hash check).
- **`distribution_level` mode**: an unseeded path (fresh OS entropy, RSI protocol) for one cross-check that reproducibility holds at the distribution level, per RSI App C.2 and Phase-4 §4.3.

---

## 5. Validation plan

### 5.1 Analytical ground-truth tests (assert against closed forms only — no self-referential fixtures)

| Test | Assertion | Source |
|---|---|---|
| Cards equilibrium | `mean AD over R eq. shuffles → (n²−1)/(3n)` within MC error | Prop 2.6 |
| Cards mean law | `mean AD(t) → (1−(1−2/n)^t)·AD_∞` across `t` | kernel [EST] |
| Cards SNR | measured `ρ²(1) → 2` | RSI §6 |
| Unbiasedness | `E[Ô_m]=O(x)` per `(n,m)`; **`Ô_n=O` identically** (exact, not statistical) | Lemma 2.10 |
| Exact variance | measured `Var(Ô_m)` matches `(n−m)/(m(n−1))·E_σ[V_σ]+Var_σ(AD)` | App E.1/E.3 |
| Chain exact | `ρ²(ϕ)=ϕ/(2−ϕ)`, `Var(ϕ̂_m)=ϕ(2−ϕ)/m`, `ρ²(1)=1`, `a=0` | [EST] |
| Ising GT | Wolff energy/`|M|` vs **Onsager** at test `T` (off-critical exact; near-`T_c` within finite-size) | Onsager |
| Ising dynamics | single-site Glauber detailed-balance sanity (acceptance ratio symmetry) | — |
| Identity | `Var(Ô_m) == sampling term + floor` across `(n,m)`, all systems | ★ |
| Veil | access audit = `m` reads; relabelling invariance; `O(Rm)` memory ceiling honored | §5.5 |

### 5.2 Phase-3 regression suite (the gate)

Reproduce **within Monte Carlo error**, then STOP-gate before Stage 1+:
- Cards `ρ²(ϕ)` curve at `ϕ∈[0.47,1]`.
- `n`-independent precision-law ratios (cards `128→2048`, chain `256→65 536`).
- Detection `m*`: cards `m*δ²≈3.1`, chain `≈6.3`.
- Ising **high-T** subsample variance `∝1/m` (ratios `→1` with `n`); Ising **near-T_c** `m`-independent floor (extra sampling doesn't reduce it).

**⚠ Gate caveat (modelling note M4):** Phase-3's near-`T_c` floor *values* (`0.56–0.58`) used **finite burn-in** and are labelled **qualitatively reliable only** by `Phase3_findings.md §4`. The new instrument uses **Wolff ground truth**, which *corrects* that bias by design. Therefore the regression test for near-`T_c` asserts the **qualitative pattern** (floor dominates; `Var(M̂_m)` flat in `m`; ratio-vs-`Var(site)/m` blows up), **not** the exact Phase-3 numbers. Reproducing the approximate values would actually indicate we reintroduced the burn-in bug. The well-established quantities (cards/chain `m*`, high-T `∝1/m`) are asserted numerically within MC error.

---

## 6. Mapping predictions & primary comparisons → runnable experiments

| Prediction | Experiment (systems, grids) | Measured → scored against |
|---|---|---|
| **P1** floor = whole-system fluctuation | S1–S4 equilibrium; sweep `m` incl. large `m`; measure `Var(Ô_m)` plateau and independent `Var_ν(O)` | plateau within ±10% of `Var_ν(O)` |
| **P2** exponent controls floor | S3 at `T_c`, `L∈{16..128}`; `R_n` vs `n` for `|M|` (a=7/8) and `E` (a≈0) | slope `log R_n` vs `log n` within ±0.1 of `a−1` |
| **P3** size-free `m*` (a=0) | S1,S2,S3-high-T; `m*(δ,n)` vs `n` | `|b|<0.1`, CI excludes 0.25; `m*δ²` within ±10% of `κ²/ρ²` |
| **P4** `m*` breaks (a>0) | S3 at `T_c`, several `δ`; infeasibility boundary | boundary tracks `κ√R_n` within ±15% |
| **P5** `(ρ²,a)` not the system | `m*δ²` across {chain,cards,Ising-hiT}; Ising `M` vs `E` at `T_c` | constants agree; two observables in different regimes |
| **P6** correlation range not exchangeability | Ising high-T (a=0) vs Curie–Weiss `T<T_c` (a=1) | high-T scale-free; CW size-limited |
| **P7** temperature crossover | S3 `T/T_c∈{2.0,1.4,1.15,1.05,1.0}`+`{0.85}`; `m*(δ,n)` | monotone onset tracking `ξ(T)/L` |
| **P8** information saturates | S1–S4; `I_eff(m)=1/CV²` vs `m` | fitted `I_max` matches `1/R_n` within ±15%; `m_×=1/(ρ²R_n)` |

**Primary comparisons** (Phase-4 §4.4): (1) `m*` flat vs `n` for S1/S2/S3-hiT; (2) `m*` vs `n` at `T_c` for `|M|`; (3) `M` vs `E` at `T_c`; (4) **2×2 factorial** {ρ² bounded (S2) vs degenerate (S2′)} × {χ summable (S2) vs extensive (S4<T_c)}; (5) floor plateau `=Var_ν(O)` across S1–S4. Each writes a machine-readable verdict (support/falsify/inconclusive) with measured value + bootstrap CI vs the pre-registered band.

---

## 7. Stage 0 concrete deliverables (what I build first)

Per the handoff staging (build before you run; gated on the regression suite):

1. `core/` (interfaces, rng/seed protocol, types incl. the `Observation` veil boundary).
2. `systems/`: **S1 product_chain**, **S2 cards**, **S3 ising2d** (Wolff GT + Glauber dynamics).
3. `observables/`: `cards_ad`, **`cards_fixedpoint` (S2′ — degenerate-SNR interface proof, exercised from day one)**, `chain_mag`, `ising_obs` (`|M|`, energy, block).
4. `observe/veil.py` (O(Rm) sampler + exchangeability fast-path + access audit).
5. `inference/` estimator + bayes + detect + mstar (both m* kinds).
6. `measure/diagnostics.py`; `storage/`; `analysis/{predictions,bootstrap,score,figures}` (enough to score Stage 0).
7. `tests/`: analytical suite + veil suite + identity suite + Ising-GT suite + **Phase-3 regression suite**.

**Exit gate:** analytical suite green **and** Phase-3 regression green (numeric where [EST]/[EMP], qualitative for near-`T_c` per M4). No Stage-1 runs until then.

Later stages register cleanly onto the same interfaces: **S4 Curie–Weiss** (Stage 3, needed for the 2×2 factorial + CX1), **S5 1D-Ising tunable ξ** and **S6 SSEP** (Stage 4). The interface demonstration: each is a new `systems/*.py` implementing the `System` ABC + an observable; no change to `observe/`, `inference/`, `measure/`, `storage/`, `analysis/` — this is the decoupling proof the handoff asks for.

---

## 8. Modelling subtleties found while cross-reading (flag before coding)

None of these reopen the frozen science; they are design decisions I've made to keep the instrument faithful. Flagged so you can veto.

- **M1 — Cards `AD` is not a "bounded, Θ(1)-intensive" observable.** Phase-4's generic premise is `D_i∈[0,B]` bounded and `O_∞=Θ(1)`. But for cards `D_i=|π(i)−i|∈[0,n−1]` grows with `n`, and `O_∞=AD_∞=(n²−1)/3n=Θ(n)`. The framework survives because **only dimensionless ratios** (`ρ²`, `R_n`, `CV`, `I_eff`) are claimed scale-free, and the exact card laws hold independently (App E). **Decision:** the code works entirely in dimensionless quantities and never assumes bounded-by-constant `D_i` or `O_∞=Θ(1)`. **Consequence for P2:** the exponent `a` is extracted from **`R_n≍n^{a−1}`** (regress `log R_n` vs `log n`), *not* from a raw `χ_n[O]≍n^a` — because the two definitions coincide only when `O_∞=Θ(1)`. All docs already phrase P2 as "slope of `log R_n` vs `log n` = `a−1`", so I follow that; `χ_n` is stored as a secondary diagnostic. This is the principled reading, but I'm flagging it because it's the one place a naïve `χ_n`-based `a` would give nonsense for cards.
- **M2 — Two distinct `m*`.** RSI's `m*(ε)=z²_{1−α/2}/(2ε²)≈1.92/ε²` (cards) is an **estimation / credible-interval** cost (`z_{0.975}=1.96`). Phase-3/4's `m*(δ)=κ²/(ρ²δ²)≈3.1/δ²` (cards) is a **detection / power** cost (`κ=z_{1−α}+z_{1−β}≈2.487`). They differ *by design*. The simulator implements **both**, clearly separated; the Phase-5 headline (P3–P5) is the **detection** one. Flagging so the two constants (`1.92` vs `3.1` for cards) are never treated as a discrepancy.
- **M3 — S2′ fixed-point fraction breaks `O_∞=Θ(1)` too** (`O_∞=E[F]=1/n→0`), and is the intended **degenerate-SNR control** (`ρ²≍1/n`, condition (i) failing). The code must keep `ρ²`, `R_n` well-defined as `O_∞→0` (no division blow-ups; use raw moments, guard/annotate the degenerate regime). Its `a`-label ("0") means "condition (ii) is *not* the binding one here" — the size-limitation comes from `ρ²→0`, the "inferability exponent" axis, not from a correlation floor. Handled in analysis; flagged so the `O_∞→0` edge cases aren't mistaken for bugs.
- **M4 — Phase-3 near-`T_c` numbers are approximate (finite burn-in).** Covered in §5.2 gate caveat; restated here as the one place the regression suite is deliberately qualitative, because the new Wolff ground truth *should* change those values.
- **M5 — `|M|` per-component readout.** `|M|` is a function of the global spin sum, not an average of per-site absolute values, so the natural additive observable is the **signed mean spin** `(1/n)Σs_i` with `|·|` applied at population level; the veil samples signed `s_i`. Below `T_c` / at criticality the sign-symmetry (`±M` degeneracy) is handled in analysis (work with `|Ô_m|` or fold the sampling distribution). Flagged because it's a genuine choice in how the estimator relates to `|M|`; I follow the "signed readout, fold at analysis" convention, consistent with treating magnetization as the intensive additive observable.

**Defaults I've chosen** (state now, easy to change): Python 3.13, NumPy/SciPy/pandas/**pyarrow**(+CSV fallback)/Matplotlib per RSI tooling; colorblind-safe palette; `α=0.05` one-sided detection, `β=0.20` (`κ²≈6.18`); bootstrap `≥1000` resamples; smoke runs at reduced `R` before full `R`.

---

## 9. Genuinely open — do any need your decision?

I don't think any of M1–M5 require a decision (they're the faithful readings), but two packaging/gate points are worth a yes/nod:

1. **Regression gate for near-`T_c` = qualitative (M4).** I'll gate on the *pattern*, not Phase-3's approximate floor values, since Wolff corrects the burn-in bias. If you instead want a numeric near-`T_c` target, it would have to be a *newly* Wolff-generated reference (I can produce and freeze one in Stage 0). Default: qualitative gate.
2. **Package name / location.** Proposed `rsi_partial/` as an installable package inside the project root. Fine, or prefer a different name / a `src/` layout?

Everything else follows the frozen spec directly. On your nod (and answers to the two above if you want to steer them), I proceed to **Stage 0** implementation exactly as scoped in §7, build the analytical + Phase-3 regression suites, and do not launch any Phase-5 run until that gate is green.
