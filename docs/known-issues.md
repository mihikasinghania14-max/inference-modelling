# Known issues, discrepancies and open questions

Every item below was found while assembling this repository. **None has been corrected in the papers, data or code.** Each states what was observed, the evidence, how confident the reading is, and where the relevant files are. "Verified here" means that the check was re-run when the repository was assembled and its output is stored in this repository.

Items are grouped by area. Severity is the assembler's judgement of how much a reader of the papers could be misled: **A** can change a stated number or claim; **B** a documentation or consistency problem; **C** minor.

---

## D. Data and simulation conventions (Stage 1)

### D1 (A). RSI Figure 2 data appear to come from a different chain than the paper describes
* **Paper.** Definition 2.1 and Appendix C.2 say the two positions `a_t, b_t` are drawn *independently*, so `a_t = b_t` occurs with probability `1/n` (self-loops), and that drawing a distinct pair "would not reproduce the chain exactly".
* **Data.** The Figure 2 data (`data/original-research/workbooks/figure2.xlsx`, sheet `G1` of `RSI Data.xlsx`) give for `n = 2` a mean AD of exactly **1.0 with standard deviation 0**. For the independent-pair chain started from the identity with `t_mix = 1` step, the mean is `0.5` with sd `0.5`. A chain that always swaps two distinct positions gives exactly 1.0 with sd 0.
* **Further evidence (verified here).** At `n = 4` (`t_mix = 3`) the sd in the Figure 2 data is 0.518; simulation gives 0.522 for the distinct-pair chain and 0.561 for the independent-pair chain. At `n = 4` the mean is 1.224 (data) against 1.204 (distinct) and 1.093 (independent).
* **Scripts.** Several scripts in `experiments/original-research/legacy-scripts` (`transposition_shuffle.py`, `card1_position.py`, `transposition_sim.py`, `early-drafts/*`) draw distinct pairs with `random.sample` or a rejection loop; `bayesian_simulation.py` and the variance grid use independent pairs (grid at `n = 2`, `τ = 1`: mean 0.479, variance 0.2498, which matches the independent-pair chain).
* **Reading.** The Figure 2 / `RSI Data.xlsx` series were very probably produced with a distinct-pair chain. At `n ≥ 64` the two conventions are indistinguishable at the precision shown, so the paper's qualitative conclusions are not affected; the small-`n` points of Figure 2 and the claim of Appendix C.2 are inconsistent with each other. This is an inference from the numbers; the script that produced the sheet was not identified with certainty.
* **Code.** `rsi_shuffle.shuffle` implements both conventions (`pair_convention="independent" | "distinct"`), tested in `tests/test_rsi_shuffle_simulation.py::test_pair_conventions_differ_at_n2`.

### D2 (B). The folder `INVALIDATED 512` mixes invalid and valid material
* `data with seed, first m.xlsx`, `figure some.xlsx` and `random_transposition_shuffle.py` in that folder use the deterministic "first `m` cards" window (`tracked_cards = np.arange(m)`, `M = 307`). The RSI paper discards this sampling scheme as biased (Section 3.2, Appendix C.4). They are kept under `first-m-window-invalidated/` and are not used by any script.
* `simulation_results.csv.xlsx` in the same folder is **not** of this kind: it is the (n, m, τ) grid generated with uniform sampling without replacement (checked: mean AD at `τ = 1` agrees with `AD_∞`, relative bias −0.4% averaged over `n ≥ 64`). It is stored as `variance-grid/simulation_results.xlsx`. Its name in the original folder should not be read as invalidating it.
* The first-`m` mean at `n = 512`, `m = 32` is 238.9 against `AD_∞ = 170.7`; the bias is large and in the direction the paper predicts.

### D3 (B). The stored grid differs from the paper's grid table and from the "72 configurations" statement
* Appendix C.3 (Table 1) lists 67 `(n, m)` pairs; Appendix E.3 and Figure S3 refer to "all 72 simulated configurations". The stored grid (`variance-grid`) has **77** `(n, m)` pairs at `τ = 1` because its `m` lists (as in `bayesian_simulation.py`) are denser than the table and include, for example, `m = 224` at `n = 512`. Five pairs in the table are absent from the stored grid: (256, 48), (512, 48), (1024, 96), (2048, 192), (4096, 384). Verified here.
* The paper reports the ratio of observed to predicted variance `Var/(n²/18m)` as median 0.95, range [0.76, 1.13]. For the stored grid at `τ = 1` (verified here): median 0.964 over all 77 pairs; for `n ≥ 16` median 0.955 and range [0.725, 1.098]; the largest ratio, 2.25, occurs at `n = 2`. The paper's median is reproduced; its range is not, so the figure was presumably computed on a different set of runs than the stored grid.
* The exact law (`rsi_shuffle.theory.estimator_variance`) explains the departures from `n²/18m`: the ratio tends to `1 − m/(5n)` for large `n` (tested).

### D4 (B). The variance-grid script was run without a seed
* This is stated in the paper (Section 3.1, Appendix C.2), so the stored grid cannot be regenerated bit-for-bit. `rsi_shuffle.grid` adds an optional seed; `scripts/check_variance_grid.py` checks that fresh runs agree with the stored grid and with the exact law within Monte Carlo error (passes).

### D5 (B). An uncaptioned figure differs from the published Figure 9
* `research/original-rsi-research/figures/2.png` shows `m*(ε)` plateaus labelled 710, 177, 82, 47 for ε = 0.05, 0.10, 0.15, 0.20. Figure 9 of the paper labels them 715, 186, 86, 48. The image is probably an earlier rendering from a different run. `1.png` is a higher-resolution version of the four-panel sampling-distribution figure (the paper's Figure 6).

### D6 (C). `t_mix` is the location of the cutoff, not the point at which every statistic is at equilibrium
* `experiments/extension/verification/shuffle_sim.py` (stored output in `results/outputs/extension-verification/shuffle_sim.txt`) finds that at `t = t_mix` and `n = 512`, mean AD (170.198) is within 0.3% of `AD_∞` (170.666) but the number of fixed points is 2.016 against 1 at equilibrium; at `2 t_mix` both agree. The RSI paper sets `t_mix` by rounding `½ n ln n` and states it is a leading-order approximation. Because the extension's certificate concerns one observable only, this is consistent with its own caveat that certification is necessary but not sufficient for equilibrium.

---

## C. Citations and attributions

### C1 (A). The extension attributes the variance of the footrule to a passage that does not contain it
* The extension (Section 4) writes `Var(F) = (n+1)(2n²+7)/45`, "by inclusion–exclusion on distinct index pairs (Singhania, 2026, App. E)". Appendix E of the RSI paper only states `Var = O(n³)` (Diaconis & Graham 1977). The formula itself is correct: it matches exact enumeration for `n ≤ 8` (verified here, `shuffle_algebra.py`, and in `tests/test_rsi_shuffle_theory.py` for `n ≤ 7`) and a symbolic derivation from the pair-moment sums. The integrated manuscript's audit reached the same conclusion and credits Diaconis & Graham with an independent proof.

### C2 (C). Internal numbering in the RSI paper
* The caption of Figure S5 refers to "Theorem 5.3"; the theorem proved in Section 5 is Theorem 5.2.

---

## M. Mathematical statements to read with care

### M1 (B). The RSI posterior uses the leading-order variance
* `rsi_shuffle.bayes` implements the model as written (`σ_m² = n²/(18m)`). That ignores the finite-population factor and the floor `Var_ν(AD)`. The effect on coverage is small for the shuffle but not zero: `shuffle_sim.py` (stored output `results/outputs/extension-verification/shuffle_sim.txt`; `n = 100`, equilibrium decks, nominal 95% interval, prior `σ₀ = 100`) gives coverage of the equilibrium displacement of 0.956 at `m = 4`, rising to 0.971 at `m = n = 100`. The extension's integrated manuscript reports 95.4% to 97.4% from 20,000 decks.

### M2 (B). The Cramér–Rao statement for `ÂD_m`
* "Cramér–Rao, attained by `ÂD_m`" holds inside the Gaussian model with `σ_m² = n²/(18m)`. The integrated manuscript's audit notes that an estimator exploiting the known positional means has per-card variance of order `n²/20` rather than `n²/18`, so the bound is a property of the model together with the sample mean, not of all estimators of the equilibrium displacement.

### M3 (B). Two different sample-size criteria
* RSI's `m*(ε) = z²/(2ε²)` is a two-sided posterior-width criterion. The extension's certification cost uses `κ = z_{1−α} + z_{1−β}` (one-sided test with power). They coincide in form but not in meaning (`rsi_shuffle.theory.m_star_rsi` versus `certification_cost`). At `n = 2048`, `ε = 5%` the exact finite-`n` value is 715 rather than 768 (verified here by `shuffle_algebra.py`: `m*(5%)` for `n = 1024, 2048, 4096` is 669, 715, 741), which explains the plateau in Figure 9.

---

## X. Extension: numbers that do not agree between sources

These concern the critical two-dimensional Ising model. The raw runs behind the extension paper's numbers were **not found** (see `data/extension/README.md`), so the figures below cannot be re-derived from the original data; they are compared with independent computations only.

### X1 (A). Critical `|M|` floor: 0.066 (extension) against about 0.072 (independent)
* Extension: `R_n[|M|] ≈ 0.066`, giving an irreducible relative error `√R_n ≈ 26%`.
* Exact values (verified here, `ising_exact.py`): `R[|M|] = 0.06917, 0.06995, 0.07061` at `L = 4, 6, 8`, increasing with `L`. Independent Wolff runs (verified here, seed 1, 2×10⁵ measurements): `0.0703 ± 0.0004` at `L = 16` and `0.0717 ± 0.0005` at `L = 32`. The audit of the integrated manuscript reports 0.0718 from longer runs and attributes the difference to the original runs being short, without establishing the cause.
* The project notes record an earlier Wolff measurement schedule (a fixed number of spins flipped between measurements) that biased the ensemble through size-biased stopping and gave `R[|M|] ≈ 0.031–0.037`; the audit of the integrated manuscript states that the extension's 0.066 did not come from that schedule (the author's own sampler gave 0.0713 at `L = 8`), and could not establish why the original runs read low.
* The qualitative claim (the order-parameter floor is flat in size and of order 26–27%) is unaffected; the quoted constant is probably 0.066 too low. The integrated manuscript already uses 0.072.

### X2 (A). Energy floor exponent: −0.82 (extension) against −0.864 (exact effective slope)
* The exact energy floor `R_n[E]` from Kaufman's torus partition function (verified here, `ising_exact.py`) has log–log slope **−0.8637** over `16 ≤ L ≤ 128`. The extension reports −0.82 with a 95% interval [−0.84, −0.80], which excludes −0.864. The slope is an effective exponent that approaches its asymptotic value slowly (`R_n[E]` carries a `ln n` correction), so −0.864 is not in conflict with the theory, but the extension's interval does not contain it.

### X3 (B). Susceptibility exponent `a[M]`: 0.889 (extension) against 0.877 (integrated manuscript's later runs)
* Both are consistent with the universal value 7/8 = 0.875.

### X4 (B). No Ising campaign exists in `rsi_partial`
* `src/rsi_partial/experiments/registry.py` lists the Ising grids (`L = 16 … 96`, plus 128 for the critical measurements) and `systems/ising2d.py` provides a Wolff sampler, but `experiments/run.py` implements only `snr_floor` and `detection` for the card shuffle and the product chain. The working notes (`research/research-extension/supplementary-materials/RESULTS.md`, an earlier status log) list the Ising critical-exponent measurement and the full P1–P8 campaign as "deliberately unrun". The extension paper (2 October) reports them as measured. The note `Phase5_criticality_confirmed.md` (6 September) documents a Wolff run "in the cloud container" at `L ∈ {16, 24, 32, 48, 64, 80}` (1600–2500 samples per size; `R_n[M]` slope −0.023, ≈ 0.066; `R_n[E]` slope −0.855; `χ[M]` slope +0.858) and states that a production run to `L ≈ 128–256` with bootstrap-CI slopes was still to be done. The script behind that run, and the later runs to `L = 128` reported in the extension, are not in the accessible material.

### X5 (B). Curie–Weiss values
* The extension's Curie–Weiss slopes (−0.996 above `T_c`, +0.005 below) are recorded in `RESULTS.md` as produced by the instrument, but the stored result files for that run were not found. The independent check `cw.py` verifies the mean-field limits `√2 − 1` and `π/2 − 1` and the exponents `a ≈ 0, 0.5, 1` at `T = 1.5, 1.0, 0.8` of the exact equilibrium distribution; it does not reproduce the slopes of the original run.

### X6 (B). Correlated-field and Algorithm 1 results
* The extension draft reports coverage collapsing from 91% to about 12% on a field with reflecting boundaries (per the project notes). An independent torus computation gives 95% falling to 7% (`field.py`, verified here). The original reflecting-boundary experiment was not found, so the two cannot be compared directly. Algorithm 1 is labelled "proposed, not validated" in the extension; the integrated manuscript reports a validation whose code is `field.py`'s Algorithm 1 section. The first-person descriptions of these computations in the integrated manuscript are the author's to confirm.

---

## R. Regenerated data: what it shows

The files under `data/extension/regenerated` and `data/original-research/regenerated` are new, seeded computations (`docs/provenance.md`). They were not tuned towards any paper. They bear on the issues above as follows; none of them replaces an original, and no paper has been changed.

### R1 (A). The critical-Ising regeneration agrees with the integrated manuscript and with exact values, not with the extension paper's 0.066 and −0.82
* Wolff runs at `L = 8 … 128` (two starts per size, equilibration gate passed at every size): `R_n[|M|]` is flat at 0.0707–0.0725 (slope +0.004, 95% interval [0.000, 0.008]); `R_n[E]` has slope −0.863 [−0.866, −0.860]; `χ_M` has slope +0.878 [0.875, 0.881].
* These sit on the exact values (`R[|M|] = 0.0692, 0.0700, 0.0706` at `L = 4, 6, 8`; Kaufman energy slope −0.8637) and on the integrated manuscript's Table of floors (`|M|` mean 0.0718, slope +0.002; `R_E` at `L = 128` 3.96 × 10⁻⁴ against 3.97 × 10⁻⁴ here). They exclude the extension paper's 0.066 and its energy interval [−0.84, −0.80]. The cause of the original runs' lower values is still not established (X1, X2).

### R2 (B). Stage 1 calibration under both pair conventions
* `data/original-research/regenerated/calibration_pair_conventions.csv` runs the shuffle to `t_mix` under both conventions (`R = 1000` each). At `n = 2` the distinct-pair chain gives exactly 1.0 with sd 0, as in the Figure 2 data, and the independent-pair chain gives 0.50 with sd 0.50; at `n ≥ 64` the two conventions agree within Monte Carlo error. This is the same observation as D1, now reproducible from a seed. It does not identify which script produced the Figure 2 sheet.

### R3 (B). Seeded variance grid
* `variance_grid_seeded.csv` has the same 1617 rows as the stored grid, run with the independent procedure of the original script and a fixed seed. At `τ = 1` and `n ≥ 64` its ratio of simulated variance to the exact law has mean 1.00 (see `data/original-research/regenerated/README.md`). The `m` lists are those of the original script, not the paper's Table 1 (D3).

### R4 (C). Mean law
* `trajectory_n512.csv` includes the closed-form mean law `φ(t) = 1 − (1 − 2/n)^t` next to the simulated `AD(t)/AD_∞`; the largest difference over the trajectory is 0.002 in units of `AD_∞` (at most 2.3 standard errors of a point), so the mean law is not contradicted.

---

## G. Repository hygiene points

* **G1 (resolved by redaction).** The manifests under `results/outputs/extension-cards-smoke` recorded the host name of the machine that generated them and its Windows build number. The host is now `[redacted]` and the platform reads `Windows-11`; nothing else in the manifests was changed, and `scripts/reproduce_extension_cards_smoke.sh` still reports identical results because it compares configuration and seed only. The stored Parquet files contain no host information.
* **G2 (revised copy provided).** The 5 November manuscript source carries the placeholder `AFFILIATION AND E-MAIL TO BE COMPLETED` and `volume/issue TBA`. `research/research-extension/integrated-manuscript/Singhania_Bayesian_Inference_revised.tex` (with `.bbl` and `.pdf`) replaces them for an arXiv-style layout; the original draft is untouched. The e-mail address of the corresponding author is not filled in. See `docs/manuscript-revisions.md`.
* **G3 (revised copy provided).** The original Data Availability statement says code and data are available "on reasonable request". The revised copy points to this repository and says which data are regenerated. The repository URL assumes the name `inference-modelling` under the author's GitHub account and the archive DOI is a placeholder until a deposit exists.
* **G4.** The working notes under `research/research-extension/` contain first-person text written during the project and refer to documents that are not in this repository (see `docs/research-history.md`). Docstrings in `src/rsi_partial` refer to a "design brief" (the build specification given at the start of Phase 5; not included, as it is a task description rather than a research document), to "Phase-N" sections and to "modelling notes M1–M5"; the latter two correspond to the working notes `Phase4_conjecture_and_design.md` and `Phase5_architecture.md`, which are included. In those docstrings "handoff" was replaced by "design brief" by a mechanical text substitution; no code was changed.

* **G5.** The LaTeX source stored with the extension paper is the 15 September version, not the source of the 2 October PDF (see `research/research-extension/README.md`).
