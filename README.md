[README.md](https://github.com/user-attachments/files/34881668/README.md)
# Inference Modelling

Code and data accompanying *Bayesian Inference in Partially Observed Critical Systems*, M. Singhania, arXiv:2610.35973.

## Overview

The paper studies how much of a large random system must be observed, through a uniformly random subset of `m` of its `n` components, before a global statistic can be inferred with controlled uncertainty. For the random transposition shuffle, the relative uncertainty of the mean-displacement estimator is `1/√(2m)`, independent of `n`. The law of total variance separates a sampling term, which decreases as `1/m`, from an `m`-independent term, the floor. At the critical point of the two-dimensional Ising model the floor of the energy density vanishes with system size, while the floor of the magnetisation magnitude does not.

The work extends a project carried out at the Research Science Initiative (RSI), India, under the supervision of Prof. Arvind Ayyer (report dated 9 July 2026).

## Repository structure

| Path | Contents |
|---|---|
| `src/rsi_shuffle` | Random transposition shuffle, estimators, Bayesian posterior, closed-form variance and floor |
| `src/rsi_partial` | Measurement code for partial observation: card shuffle, product chain, 2D Ising, Curie–Weiss |
| `notebooks` | Walkthroughs of the variance law and of the floor |
| `experiments/original-research` | Original RSI scripts |
| `experiments/extension` | Independent verification scripts and the Wolff-sampler campaign for the critical Ising model |
| `data` | RSI spreadsheets and variance grid; regenerated Stage 1 and Stage 2 data |
| `results` | Figures and stored outputs |
| `research` | RSI paper, extension manuscripts and working notes |
| `docs` | Methodology, mathematical background, computational methods, reproducibility, provenance and known issues |

## Installation

Python 3.11 or later.

```bash
pip install -e ".[test]"
```

A C compiler is required for the Wolff sampler. `environment.yml` and `requirements.txt` are provided for conda and pip environments.

## Usage

```bash
pytest                  # 103 tests
make check              # fresh simulation against the stored variance grid and the exact law
make figures            # regenerate all figures
make stage1-data        # seeded Stage 1 datasets (about five minutes)
make ising              # critical Ising campaign (several hours)
```

Example:

```python
from rsi_shuffle import theory

theory.relative_floor(1024)            # 3.91e-04
theory.certification_cost(0.10, 1024)  # approximately 408 observed cards
```

Seeds and expected run times are listed in `docs/reproducibility.md`.

## Data

- `data/original-research` contains the spreadsheets underlying the RSI figures and the variance grid, with CSV exports.
- `data/original-research/regenerated` contains seeded reruns of the Stage 1 experiments: calibration under both pair conventions, the `n = 512` trajectory and the full variance grid.
- `data/extension/regenerated` contains Wolff-sampler block means for the critical Ising model (`L = 8, …, 128`, two starts per size), their analysis, and exact Curie–Weiss moments. At `T_c` the magnetisation floor is approximately 0.072 and the energy-floor exponent is −0.863, in agreement with Table 1 of the integrated manuscript.
- `results/outputs` contains the stored smoke-run output of the measurement code and the output of the verification scripts.

Files under `regenerated/` were produced for this repository with fixed seeds. They are reruns of the experiments reported in the papers and are not the original run outputs. File provenance is documented in `docs/provenance.md`.

## Citation

See `CITATION.cff`.

> M. Singhania. *Bayesian Inference in Partially Observed Critical Systems.* arXiv:2610.35973.

An archive DOI will be added upon deposit.

## Licence

The code is released under the MIT licence (`LICENSE`). The papers, figures and data are not covered by this licence.

## Acknowledgements

The author thanks Prof. Arvind Ayyer, IISc for supervision and guidance and the staff of RSI India, the Centre for Excellence in Education and the Indian Institute of Science. Computing costs of the later work were supported by an Emergent Ventures (Mercatus Center) grant.
