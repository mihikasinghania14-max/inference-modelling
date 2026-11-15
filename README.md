[README.md](https://github.com/user-attachments/files/34881558/README.md)
# inference-modelling

Code and data for *Bayesian Inference in Partially Observed Critical Systems* (arXiv:2610.35973), by Mihika Singhania, independent researcher, India.

The paper starts from a small question. You have a system of `n` pieces and you look at `m` of them, picked at random. How big does `m` have to be before you can say something about the whole system? For a deck of cards shuffled by random transpositions the answer doesn't depend on `n`. For the 2D Ising model at its critical point it does stop working, past a certain precision, and looking at more of the system doesn't help. The paper calls that leftover error the floor. The arguments are in the paper, so I haven't repeated them here.

It began as a project at the Research Science Initiative in India (supervised by Arvind Ayyer; that report is dated 9 July 2026 and is in `research/original-rsi-research`). The Ising, Curie–Weiss and correlated-field parts came afterwards, on my own.

## What's where

- `src/rsi_shuffle` – the shuffle, the estimators, the closed-form variance and floor. This is the RSI half.
- `src/rsi_partial` – the measurement code for the second half (cards, a product chain, Ising, Curie–Weiss). It has its own tests, 74 of them.
- `notebooks` – two short ones, on the variance law and on the floor.
- `experiments/original-research` – my old scripts, left as I wrote them. Some are rough.
- `experiments/extension/verification` – scripts that re-derive the algebra and a few numbers independently, with their output in `results/outputs/extension-verification`.
- `experiments/extension/campaign` – the Wolff sampler (C) and the analysis for the Ising runs. Kaufman's torus formula and Wolff's 1989 cluster update are what it leans on.
- `data` – spreadsheets from the RSI paper and the newer runs. Read the next section first.
- `research` – both papers, the manuscript drafts and my working notes from the time. The notes are untidy and written in the first person, as they were.
- `docs` – notes on things I found wrong or odd. `known-issues.md` is the one to read.

## Running it

```
pip install -e ".[test]"
pytest
```

103 tests, a minute or two. The `Makefile` has shortcuts for the rest. The Ising campaign takes hours (L = 96 and 128 are the slow ones); regenerating the Stage 1 data takes about five minutes. A C compiler is needed for the Wolff code.

## About the data

If you're checking numbers, this is the part that matters.

The Stage 1 spreadsheets are the originals. They were produced without a fixed seed, so they can't be rebuilt exactly. `data/original-research/regenerated` has the same experiments rerun with seeds. Same method, different random numbers.

For Stage 2 I could not find the raw output behind the extension figures: the Ising runs up to L = 128, the Curie–Weiss runs, the correlated fields, the Algorithm 1 runs. What's in `data/extension/regenerated` is new. Wolff runs for L = 8 to 128 and exact Curie–Weiss sums, made this month to replace the ones that went missing. They are not the runs the paper used. They agree with Table 1 of the manuscript within Monte Carlo error (the `|M|` floor is about 0.072, the energy slope about −0.863). They do not agree with my older extension paper, which has 0.066 and −0.82. I don't know why those came out lower.

The correlated-field results and Algorithm 1 were only checked on a torus, independently (`experiments/extension/verification/field.py`). Algorithm 1 is a proposal. I wouldn't call it validated.

## Things that don't line up

- The RSI Figure 2 data at n = 2 and n = 4 look like a chain that swaps two different positions, while the paper defines the chain with independent positions. The difference is gone by n = 64.
- The stored variance grid has 77 (n, m) pairs. The paper's table has 67 and its text says 72.
- `docs/known-issues.md` has the rest, with the evidence.

## Help

Part of the code, the tests, the regenerated data and the docs were written with an AI assistant (Claude, from Anthropic). `docs/provenance.md` says what and where. The RSI paper and the scripts in `experiments/original-research` are my own work.

## Citing, licence, thanks

`CITATION.cff` has the details. There is no archive DOI yet; it goes here once I've deposited the repository. The code is MIT. The papers, figures and data aren't covered by that licence.

Thanks to Prof. Ayyer for the mentoring and for leaving room to explore, to Keerthana Senthil, to Tarun Prasad for reading the work, to Amy Sillman, Max Bee Lindgren and the RSI staff, and to the Centre for Excellence in Education and IISc. The later computing was paid for by an Emergent Ventures (Mercatus Center) grant.
