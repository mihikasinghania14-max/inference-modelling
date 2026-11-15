"""
Random Transposition Shuffle — Average Displacement Experiment
==============================================================
Research goal: Measure how close a deck gets to equilibrium after t_mix
random transpositions, comparing simulated average displacement against
the theoretical value AD_theory = (n^2 - 1) / (3n).

Key improvements over the starter code
---------------------------------------
* Uses random.sample() instead of a while-loop to pick two distinct positions.
* Tracks card positions via a companion `positions` array so we never call
  list.index() (O(n)) inside the hot loop; all displacement lookups are O(1).
* Separates logic into five focused functions with clear docstrings.
"""

import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


# ---------------------------------------------------------------------------
# 1. Core shuffle primitive
# ---------------------------------------------------------------------------

def random_transposition_shuffle(deck: list[int], positions: list[int]) -> None:
    """
    Perform one random transposition in-place.

    Picks two distinct positions i, j uniformly at random (using
    random.sample for efficiency) and swaps the cards there.
    The companion `positions` array is updated simultaneously so that
    positions[card] always gives the current 0-based index of that card.

    Parameters
    ----------
    deck      : list of card values (mutated in place)
    positions : list where positions[card_value] = current index (mutated in place)
                Cards are 0-indexed here (card value - 1) for direct array access.
    """
    n = len(deck)
    # random.sample draws without replacement — gives two *distinct* indices
    i, j = random.sample(range(n), 2)

    # Update the position lookup for both cards before the swap
    positions[deck[i]] = j
    positions[deck[j]] = i

    # Swap the cards in the deck
    deck[i], deck[j] = deck[j], deck[i]


# ---------------------------------------------------------------------------
# 2. Average displacement metric
# ---------------------------------------------------------------------------

def average_displacement(positions: list[int]) -> float:
    """
    Compute average displacement of all cards from their original positions.

    AD = (1/n) * sum( |current_position[card] - original_position[card]| )

    Because the deck starts as [0, 1, 2, ..., n-1] (0-indexed internally),
    the original position of card k is simply k.  So we just need
    |positions[k] - k| for each k.

    Parameters
    ----------
    positions : current position lookup (positions[card] = index in deck)

    Returns
    -------
    float  Average displacement value.
    """
    n = len(positions)
    # Vectorised absolute difference via numpy for speed at large n
    idx = np.arange(n, dtype=np.int64)
    return float(np.mean(np.abs(np.asarray(positions, dtype=np.int64) - idx)))


# ---------------------------------------------------------------------------
# 3. Single trial
# ---------------------------------------------------------------------------

def run_single_trial(n: int, t_mix: int) -> float:
    """
    Run one trial: create an ordered deck, apply t_mix random transpositions,
    return the resulting average displacement.

    Internal representation uses 0-based card values for direct array indexing:
        deck[i]      — which card is at position i
        positions[k] — which position card k is at

    Parameters
    ----------
    n     : deck size
    t_mix : number of random transpositions to apply

    Returns
    -------
    float  Average displacement after t_mix transpositions.
    """
    # Ordered deck: card at position i has value i (0-based)
    deck = list(range(n))
    # Position lookup: card k starts at position k
    positions = list(range(n))

    for _ in range(t_mix):
        random_transposition_shuffle(deck, positions)

    return average_displacement(positions)


# ---------------------------------------------------------------------------
# 4. Full experiment across all deck sizes
# ---------------------------------------------------------------------------

def run_experiment(
    deck_sizes: list[int],
    R: int = 500,
) -> pd.DataFrame:
    """
    Run the shuffling experiment for every deck size and collect statistics.

    For each n:
      - Compute t_mix = ceil( (1/2) * n * ln(n) )
      - Repeat the single-trial experiment R times
      - Record mean AD, std AD, and theoretical AD_theory = (n^2 - 1) / (3n)

    Parameters
    ----------
    deck_sizes : list of deck sizes to test
    R          : number of independent trials per deck size

    Returns
    -------
    pd.DataFrame with columns:
        n, t_mix, mean_AD, std_AD, theoretical_AD
    """
    records = []

    for n in deck_sizes:
        # Theoretical mixing time (number of transpositions)
        t_mix = math.ceil(0.5 * n * math.log(n))

        # Theoretical equilibrium average displacement
        ad_theory = (n ** 2 - 1) / (3 * n)

        print(f"  n={n:>5d} | t_mix={t_mix:>8d} | running {R} trials ...", end="", flush=True)

        # Run R independent trials and collect ADs
        ads = [run_single_trial(n, t_mix) for _ in range(R)]

        mean_ad = float(np.mean(ads))
        std_ad  = float(np.std(ads, ddof=1))   # sample std dev

        print(f"  mean_AD={mean_ad:.4f}  std_AD={std_ad:.4f}")

        records.append({
            "n":              n,
            "t_mix":          t_mix,
            "mean_AD":        mean_ad,
            "std_AD":         std_ad,
            "theoretical_AD": ad_theory,
        })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# 5. Plotting
# ---------------------------------------------------------------------------

def plot_results(df: pd.DataFrame, output_path: str = "shuffle_results.png") -> None:
    """
    Produce a publication-quality plot comparing simulated and theoretical AD.

    - x-axis: deck size n on a log-2 scale
    - y-axis: average displacement
    - Blue points with error bars: simulation mean ± std deviation
    - Orange line: theoretical equilibrium AD_theory = (n^2 - 1) / (3n)

    Parameters
    ----------
    df          : results DataFrame from run_experiment()
    output_path : file path to save the figure
    """
    fig, ax = plt.subplots(figsize=(9, 5))

    # --- Simulation results as scatter + error bars ---
    ax.errorbar(
        df["n"], df["mean_AD"],
        yerr=df["std_AD"],
        fmt="o",
        color="#2563EB",          # blue points
        ecolor="#93C5FD",         # light-blue error bars
        elinewidth=1.5,
        capsize=4,
        capthick=1.5,
        markersize=6,
        label="Simulation mean ± std",
        zorder=3,
    )

    # --- Theoretical curve as a smooth line ---
    ax.plot(
        df["n"], df["theoretical_AD"],
        color="#F97316",           # orange line
        linewidth=2,
        marker="",
        label=r"Theory: $\frac{n^2-1}{3n}$",
        zorder=2,
    )

    # --- Axis scaling and labels ---
    ax.set_xscale("log", base=2)

    # Show actual n values as tick labels (not powers)
    ax.set_xticks(df["n"].tolist())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.tick_params(axis="x", rotation=45)

    ax.set_xlabel("Deck size $n$  (log₂ scale)", fontsize=12)
    ax.set_ylabel("Average displacement (AD)", fontsize=12)
    ax.set_title(
        r"Average Displacement after $t_{\mathrm{mix}}$ Random Transpositions"
        "\n"
        r"$t_{\mathrm{mix}} = \lceil \frac{1}{2} n \ln n \rceil$,  $R = 500$ trials per $n$",
        fontsize=13,
    )

    ax.legend(fontsize=11)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.6)
    fig.tight_layout()

    fig.savefig(output_path, dpi=150)
    print(f"\nPlot saved to: {output_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    DECK_SIZES = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
    R = 500   # trials per deck size

    print("=" * 65)
    print("Random Transposition Shuffle — Average Displacement Experiment")
    print("=" * 65)
    print(f"Deck sizes : {DECK_SIZES}")
    print(f"Trials (R) : {R}\n")

    df = run_experiment(DECK_SIZES, R=R)

    # --- Print formatted results table ---
    print("\n" + "=" * 65)
    print("Results Table")
    print("=" * 65)
    print(
        df.to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}",
            columns=["n", "t_mix", "mean_AD", "std_AD", "theoretical_AD"],
        )
    )

    plot_results(df, output_path="/mnt/user-data/outputs/shuffle_results.png")


if __name__ == "__main__":
    main()
