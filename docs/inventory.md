# Inventory of source material

Result of the file audit made when the repository was assembled. "Source" is the folder name on the author's machine (the folder `outputs` stands for the author's working folder of generated files); nothing there was modified. Sources searched: all ten connected folders of the author's computer (including nested folders and the `Data Files` folder), the attached project, and the uploaded files. No Jupyter notebooks were found anywhere; the two notebooks in `notebooks/` were written for this repository.

## Included

| Source | Destination | Notes |
|---|---|---|
| `Inference Modelling/RSI_Paper_Mihika_Singhania.pdf` | `research/original-rsi-research/paper/` | The RSI paper (28 pages). Identical to the copy in the project. |
| `RSI figures/*` | `research/original-rsi-research/figures/` | Author's figure files; see its README. About 20 MB in total; Git LFS not needed. |
| `Inference Modelling/Extension_Paper_MihikaSinghania .pdf` | `research/research-extension/paper/Extension_Paper_Mihika_Singhania.pdf` | 2 October extension PDF; identical in bytes to the project's copy. File name normalised (space removed). |
| `extension/{main.tex, references.bib, fig_*.png}` | `research/research-extension/paper/source/` | 15 September version; `supplement-2.tex` in the outputs folder is byte-identical to `main.tex`. |
| `outputs/Singhania_Bayesian_Inference_BA.pdf`, `outputs/JSP_source_BA.zip` (unpacked) | `research/research-extension/integrated-manuscript/` | Without `imsart.cls`, `imsart.sty`, `ba2.bst` (third-party). A revised copy (`…_revised.*`) with the data-availability and front-matter fixes was added on 15 November; see `docs/manuscript-revisions.md`. |
| `RSI IRIS/outputs/{Phase3_findings, Phase3_interpretation, Phase4_conjecture_and_design, Phase5_criticality_confirmed, Phase6_theory, EXTENSION_PAPER_v2}.md` | `research/research-extension/theoretical-framework/` | Working notes, verbatim. |
| `RSI IRIS/outputs/{Phase5_architecture, Phase7_paper_skeleton, RESULTS}.md` | `research/research-extension/supplementary-materials/` | Working notes, verbatim. |
| `data files/*.xlsx` (`RSI PAPER/data files`) | `data/original-research/workbooks/` | Plus CSV exports. `appendix fig1.xlsx` is byte-identical to `figure4.xlsx` and is not duplicated. |
| `INVALIDATED 512/simulation_results.csv.xlsx` | `data/original-research/variance-grid/simulation_results.xlsx` | A valid grid (see `known-issues.md` D2). |
| `INVALIDATED 512/{data with seed, first m.xlsx, figure some.xlsx}` | `data/original-research/first-m-window-invalidated/` | Biased first-`m` window; kept as evidence. |
| `RSI PAPER/python codes here/*.py` (9 files; `hi.py` omitted) and `INVALIDATED 512/random_transposition_shuffle.py` | `experiments/original-research/legacy-scripts/` | Verbatim except two renames. |
| `python codes/{one.py, one.txt, two.py, test.py}` | `legacy-scripts/early-drafts/` | Early exploratory snippets. |
| `RSI IRIS/outputs/rsi_partial/` (package and 74 tests) | `src/rsi_partial/` | Code unchanged except that the word "handoff" in 11 docstrings was replaced by "design brief"; `pyproject.toml` and README rewritten. |
| `RSI IRIS/results/*` (parquet, manifests) | `results/outputs/extension-cards-smoke/` | Unmodified. |
| `outputs/verification_scripts.zip` (unpacked) | `experiments/extension/verification/` | Scripts unchanged; `wolff_analyze.py` and `run_all.sh` added; outputs stored in `results/outputs/extension-verification`. |
| `RSI IRIS/play.py` | not included | An interactive demo script around the instrument; not needed to reproduce anything. |

Written for this repository: `src/rsi_shuffle`, `tests/`, `scripts/`, `notebooks/`, `experiments/extension/campaign`, `data/*/regenerated` (generated, see `docs/provenance.md`), `Makefile`, `.github/workflows`, the documentation, `LICENSE`, `CITATION.cff`, `pyproject.toml`, `requirements.txt`, `environment.yml`, `.gitignore`.

## Found but not included

| Source | Reason |
|---|---|
| `RSI/` (milestone report, milestone presentation, mini-paper with comments, overall feedback, experiment-data notes, a temperature-chart page, programme milestone PDF) | RSI-programme administrative documents with comments from mentors and tutors; not research outputs. The author may add any of them. |
| `Inference Modelling/IRIS_video_graphics/` | Graphics and video clips for a presentation video; not research data. |
| `Inference Modelling/Probe_live_instrument.html`, `Probe_source.zip`, `RSI IRIS/outputs/{probe.html, Probe_netlify.zip, pm_*.png, fig1…fig7 *.png}` | An interactive web page and presentation figures that present the extension's results (the page computes its numbers in the browser from a seed). A presentation tool rather than a research record; the author should decide whether to publish it separately. |
| `outputs/{Change_Log_and_Author_Flags, Referee_Audit_Bayesian_Inference_Partially_Observed_Critical_Systems}.md`, `RSI IRIS/{ADVERSARIAL_REVIEW, SUBMISSION_HARDENING}.md`, `derivations_study_guide.md`, `rsi_vs_extension_and_defense.md`, `interview_defense.html` | Audit, review and study-guide documents written for the author and for competition preparation. The discrepancies they record are re-checked and summarised in `known-issues.md`. |
| `EV_*.md`, `EV_grant*.{md,docx}`, `EV_application_draft_v1.md`, `ev_interview_prep.md` | Grant and application material with personal content. |
| the build-specification and migration files (`*HANDOFF*.md`) | Task and migration descriptions; the latter also contains personal contact details. |
| `supplement.tex`, `supplement-1.tex`, `supplement.pdf`, `Singhania_supplement*.pdf`, `Singhania_revised*` , `EXTENSION_PAPER.{pdf,html}`, `fig_criticality_L128.png` | Earlier versions of the extension and of the integrated manuscript, superseded by the files included. |
| `imsart.cls`, `imsart.sty`, `ba2.bst` | Third-party LaTeX class and style files of the Institute of Mathematical Statistics. |
| `hi.py` | One-line print statement. |
| `__pycache__`, `.pytest_cache`, `.git`-like or editor files | Generated or system artefacts. |

## Checks made

* Searched all included text, code and data files for credentials, tokens, e-mail addresses, telephone numbers and absolute personal paths: none found, with two exceptions that are documented rather than changed: the manifests record the machine host name (`known-issues.md` G1) and the legacy scripts contain paths such as `/mnt/user-data/outputs/…` that do not exist on other machines.
* Spreadsheet properties contain the author's name only.
* No file exceeds 6 MB; Git LFS is not needed.
* No transcripts, prompt text or tool logs are present in the repository (the repository text was searched case-insensitively). The only mentions of the AI assistant are the disclosure in `docs/provenance.md` and the pointers to it. The legacy script whose file name contained a tool name was renamed, and one file-name reference to the build specification inside `Phase5_architecture.md` was replaced by the words "the build specification" (the only edit made to the working notes).
