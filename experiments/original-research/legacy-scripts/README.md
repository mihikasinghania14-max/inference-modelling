# Legacy scripts (preserved as written)

These are the scripts found in the author's working folders. They have **not** been edited, apart from two file renames (listed below). Several have hard-coded paths, interactive plots or no seed; use `src/rsi_shuffle` for new work.

| File | Contents | Notes |
|---|---|---|
| `bayesian_simulation.py` | Generates the `(N, M, τ)` grid CSV (1617 rows) with independent pair draws and uniform `m`-subsets without replacement | No seed. Writes `simulation_results.csv` in the working directory. Its `N` and `M` lists are identical to those in `rsi_shuffle.grid`. |
| `random_transposition_shuffle.py` | One `(N, M)` configuration (`N = 512`, `M = 1` as saved) with histogram and normal fit; uniform sampling without replacement | Calls `plt.show()`. `M` and `N` are edited by hand. |
| `code_runner.py` | Same procedure with a parametrised number of steps `T` (`N = 512`, `M = 32`, `T = 01` as saved) | Renamed from a file name that contained the name of the tool used to draft it. Content unchanged. |
| `transposition_shuffle.py` | Calibration experiment: mean `AD` after `t_mix` against `(n²−1)/(3n)`; `random.sample` picks **distinct** positions | Writes to `/mnt/user-data/outputs/…`; see `docs/known-issues.md` D1. |
| `card1_position.py` | Final position of card 1, `N = 1024` | Distinct pairs. |
| `card1_position_v2.py` | Final position of card 1, `N = 512`, with histogram | Renamed from `card1_position (1).py`. Writes to a hard-coded path and calls `plt.show()`. |
| `transposition_sim.py` | Small deck (`N = 50`), `T = 98` | Distinct pairs. |
| `ad_plot.py` | Redraws RSI Figure 2 from a spreadsheet `plot.xlsx` | Reads an absolute path that does not exist here. Superseded by `scripts/make_original_figures.py`. |
| `first-m-window-invalidated/random_transposition_shuffle_first_m.py` | As `random_transposition_shuffle.py` but tracks cards `0…m−1` (`M = 307`) | **Biased sampling scheme, discarded in the paper.** Kept to document the invalidated data in `data/original-research/first-m-window-invalidated`. |
| `early-drafts/` | Early exploratory snippets (`one.py`, `one.txt`, `two.py`, `test.py`) | Incomplete (parts of `two.py` are commented out); not used for any reported result. Distinct pairs. |

Omitted: `hi.py`, a one-line print statement.
