# Revisions to the integrated manuscript for the data release

`research/research-extension/integrated-manuscript/Singhania_Bayesian_Inference.{tex,bbl,pdf}` is the 5 November draft, kept as found. `Singhania_Bayesian_Inference_revised.{tex,bbl,pdf}` is a copy with the changes below and nothing else; the mathematics, tables and figures are identical. It was compiled with the *Bayesian Analysis* class (`imsart`, not included here) to 17 pages with no undefined references.

| # | Draft | Revised copy |
|---|---|---|
| 1 | `\documentclass[ba]{imsart}` with `\volume{TBA}`, `\issue{TBA}`, `\firstpage{1}`, `\lastpage{1}` (the header read "Number TBA, pp. 1") | `\documentclass[ba,preprint]{imsart}`; the four placeholder lines removed (header without volume or issue) |
| 2 | `\address[A]{AFFILIATION AND E-MAIL TO BE COMPLETED}` | `\address[A]{Independent researcher, India}`. The e-mail line (`\ead`, `\printead`) is not filled in; add the corresponding author's address there |
| 3 | Data Availability: code and data "available from the corresponding author on reasonable request" | Code, exact-computation scripts and data at the repository; archive DOI to be added after deposit. States that the Stage 1 spreadsheets are the original files (generated without a fixed seed, with seeded re-runs alongside), that the raw output of the Ising, Curie–Weiss, correlated-field and Algorithm 1 runs was not recovered, that the repository's Ising and Curie–Weiss data are seeded re-runs agreeing with Table 1 within Monte Carlo error, and that the field and Algorithm 1 results are supported by independent torus computations |

## Still to do before submission

* Fill in the e-mail address (item 2).
* Replace the placeholder archive DOI once the repository is deposited (for example on Zenodo), and check that the GitHub URL matches the repository's final name.
* Decide the AI-use disclosure for the manuscript itself; the repository's statement is in `docs/provenance.md`.
* The Methods paragraph on the Wolff runs describes the manuscript's own runs (two chains of `2 × 10⁵` measurements per size, `L = 4 … 128`). The regenerated runs in `data/extension/regenerated` use a different, smaller schedule at large `L` (see its README, which compares both). Citing the manuscript's numbers and calling the repository's data "regenerated" is consistent; describing the repository's runs by the manuscript's protocol would not be.
