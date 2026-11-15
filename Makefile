# Convenience targets. Everything here is also documented in docs/reproducibility.md.
PY ?= python

.PHONY: install test check figures stage1-data ising cw notebooks clean

install:
	$(PY) -m pip install -e ".[test]"

test:
	$(PY) -m pytest

check:                 ## fresh Monte Carlo against the stored variance grid and the exact law
	$(PY) scripts/check_variance_grid.py

figures:               ## figures from stored and regenerated data
	$(PY) scripts/make_original_figures.py
	$(PY) scripts/make_extension_figures.py

stage1-data:           ## seeded Stage 1 datasets (about a minute)
	$(PY) scripts/regenerate_stage1_data.py

cw:                    ## exact Curie-Weiss moments
	$(PY) experiments/extension/campaign/curie_weiss_exact.py data/extension/regenerated/curie_weiss_exact.csv

ising:                 ## Wolff campaign (hours) and its analysis
	bash experiments/extension/campaign/run_ising_campaign.sh data/extension/regenerated/ising_critical_blocks 2
	$(PY) experiments/extension/campaign/analyse_ising_campaign.py \
	    data/extension/regenerated/ising_critical_blocks data/extension/regenerated/ising_critical

notebooks:             ## re-execute the notebooks in place
	cd notebooks && jupyter nbconvert --to notebook --execute --inplace *.ipynb

clean:
	rm -rf .pytest_cache build *.egg-info
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
