# Research history

Dates are taken from the documents themselves or from file modification times in the author's folders; they are indicative.

| Date | Event |
|---|---|
| June 2026 | RSI programme milestone work (mini-paper, first milestone report and presentation; programme documents, not included) and early experiments: `experiments/original-research/legacy-scripts/early-drafts`, `transposition_shuffle.py`, `card1_position*.py` (undated). |
| 9 July 2026 | RSI paper *Bayesian Inference in Random Transposition Shuffles* (supervisor Prof. A. Ayyer). Variance grid generated with `bayesian_simulation.py` (undated); the first-`m` window data are discarded as biased (Section 3.2). |
| 16 August 2026 | Working notes **Phase 3** pilots on three systems (cards, product chain, 2D Ising) and their interpretation (`Phase3_findings.md`, `Phase3_interpretation.md`); **Phase 4**, the conjecture and experimental design (`Phase4_conjecture_and_design.md`: predictions P1–P8 with support and falsify bands, experimental matrix, controls, counterexamples); **Phase 5** architecture (`Phase5_architecture.md`) and the first version of the instrument `rsi_partial` (file times). |
| 17 August 2026 | Stage 1 smoke results for the card shuffle (manifest timestamps). |
| 31 August 2026 | `RESULTS.md` (74 tests; the Curie–Weiss hard floor recorded as a new result; the Ising critical campaign listed as unrun) and `Phase6_theory.md`. |
| 6 September 2026 | `Phase5_criticality_confirmed.md`, `Phase7_paper_skeleton.md`. |
| 10 September 2026 | Longer draft `EXTENSION_PAPER_v2.md` (and a PDF/HTML rendering). |
| 14–15 September 2026 | First LaTeX versions of the extension (*From Estimation to Certification: An Irreducible Floor in Partial Observation*); figures `fig_decomposition.png`, `fig_criticality.png`. |
| 2 October 2026 | Extension paper PDF *An Irreducible Floor in Partial Observation* (5 pages). |
| October 2026 | Independent verification of the extension's algebra and numerics (`experiments/extension/verification`); an audit of the integrated manuscript records the discrepancies listed in `docs/known-issues.md`. |
| 5 November 2026 | Integrated manuscript *Bayesian Inference in Partially Observed Critical Systems* in the *Bayesian Analysis* format (draft). |
| 14 November 2026 | This repository assembled. |
| 15 November 2026 | Seeded Stage 1 datasets, Wolff campaign for the critical Ising model and exact Curie–Weiss moments added under `regenerated/`; notebooks, figures, Makefile and CI added; repository prepared to accompany the arXiv preprint (arXiv:2610.35973, as given by the author). |

## Documents referred to but not included

The code and working notes refer to a number of documents that are not part of this repository, because they are task descriptions, private application or presentation material, or audits written for the author:

* the Phase-5 build specification (the "design brief" / "handoff"), a task description for writing the instrument;
* adversarial-review, submission-hardening and referee-audit notes and a study guide for defending the work;
* application and grant drafts, interview preparation, a presentation video's graphics and an interactive web page used to present the results ("Probe").

The scientific content that the code depends on is in `research/research-extension/theoretical-framework` and `supplementary-materials`. The discrepancies found by the audits are recorded, with evidence re-checked, in `docs/known-issues.md`. Items excluded for privacy or relevance are listed in `docs/inventory.md`; the author may decide to add any of them.

## Why two stages

The RSI paper's scale-invariance result `CV = 1/√(2m)` was established for one observable of one system. The extension asked whether it is a general property, and found that the exact variance decomposition of the RSI paper (its Appendix E) contains a second term that does not shrink with `m`. The extension identifies that term as the quantity that decides whether the result generalises, evaluates it for the shuffle, and measures it at a critical point where it separates two observables of the same system. Open problems are stated in Section 6 of the extension and in `README.md`.
