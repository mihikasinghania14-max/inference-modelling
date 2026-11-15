"""
Partial Observation – Random Transposition Shuffle
Bayesian Inference Data Generator

Generates a CSV with one row per (N, M, τ) combination.
For each combination, R=1000 independent Monte-Carlo simulations are run
and summary statistics are computed.

Output columns:
    N, M, sampling_fraction, t_mix, tau, T,
    mean_AD, variance, standard_deviation, standard_error,
    CI_lower, CI_upper, skewness, kurtosis
"""

# =============================================================================
# Imports
# =============================================================================

import csv
import math
import numpy as np
from scipy.stats import kurtosis, skew

# =============================================================================
# Experiment parameters
# =============================================================================

# Deck sizes and their corresponding M (observed-card) values.
# Each list of M values corresponds to ~5%, 10%, 20%, 40% sampling fractions.
N_VALUES = [
    2,
    4,
    8,
    16,
    32,
    64,
    128,
    256,
    512,
    1024,
    2048,
    4096,
]

M_VALUES = {
    2:    [1, 2],
    4:    [1, 2, 3, 4],
    8:    [1, 2, 4, 6, 8],
    16:   [2, 4, 6, 8, 12, 16],
    32:   [2, 4, 8, 12, 16, 24, 32],
    64:   [4, 8, 16, 24, 32, 48, 64],
    128:  [8, 16, 32, 48, 64, 96, 128],
    256:  [16, 32, 64, 96, 128, 192, 256],
    512:  [32, 64, 96, 128, 192, 256, 384, 512],
    1024: [64, 128, 192, 256, 384, 512, 768, 1024],
    2048: [128, 256, 384, 512, 768, 1024, 1536, 2048],
    4096: [256, 512, 768, 1024, 1536, 2048, 3072, 4096],
}


# 21 equally-spaced normalised time points: 0.0, 0.1, ..., 2.0
TAU_VALUES = [round(i * 0.1, 1) for i in range(21)]

R = 1000          # independent simulations per (N, M, T) combination

OUTPUT_FILE = "simulation_results.csv"

# =============================================================================
# Mixing time
# =============================================================================

def compute_mixing_time(n: int) -> int:
    """
    Return the theoretical RT-shuffle mixing time.

        t_mix = round( N/2 * ln(N) )

    Parameters
    ----------
    n : deck size
    """
    return round((n / 2) * math.log(n))


# =============================================================================
# Single simulation  (unchanged computational logic from original)
# =============================================================================

def run_single_simulation(n: int, m: int, t: int) -> float:
    """
    Run one simulation: start from an ordered deck, apply t random
    transpositions, sample m cards without replacement, and return
    the average displacement of those sampled cards.

    Parameters
    ----------
    n : deck size
    m : number of cards to observe
    t : number of transpositions to apply (whole-number approximation of T = τ·t_mix)
    """
    # Step 1: identity deck
    deck = np.arange(n, dtype=np.int32)

    # Step 2: apply t random transpositions
    if t > 0:
        positions = np.random.randint(0, n, size=(t, 2))
        for i, j in positions:
            deck[i], deck[j] = deck[j], deck[i]

    # Step 3: randomly sample m distinct cards uniformly without replacement
    tracked_cards = np.random.choice(n, size=m, replace=False)
    original_indices = tracked_cards  # card k starts at position k

    # Step 4: build inverse permutation to find final positions
    inv = np.empty(n, dtype=np.int32)
    inv[deck] = np.arange(n, dtype=np.int32)

    # Final positions of the tracked cards
    X = inv[tracked_cards]

    # Average displacement of the sampled cards only
    ad = np.mean(np.abs(X - original_indices))
    return float(ad)


# =============================================================================
# Batch simulations  (unchanged from original)
# =============================================================================

def run_all_simulations(n: int, m: int, t: int, r: int) -> np.ndarray:
    """
    Run r independent simulations and return an array of Average Displacement
    values, one per simulation.

    Parameters
    ----------
    n : deck size
    m : number of observed cards
    t : transpositions per simulation
    r : number of simulations
    """
    ad_values = np.empty(r, dtype=np.float64)
    for k in range(r):
        ad_values[k] = run_single_simulation(n, m, t)
    return ad_values


# =============================================================================
# Statistics  (unchanged from original)
# =============================================================================

def compute_statistics(ad_values: np.ndarray) -> dict:
    """
    Compute summary statistics for an array of simulation results.

    Returns mean, variance, std, stderr, 95%-CI bounds, CV, skewness, kurtosis.
    The 95% CI uses the normal approximation: x̄ ± 1.96 · s / √R.

    Parameters
    ----------
    ad_values : 1-D array of average displacement values
    """
    r      = len(ad_values)
    mean   = float(np.mean(ad_values))
    var    = float(np.var(ad_values, ddof=1))   # unbiased sample variance
    std    = math.sqrt(var)
    stderr = std / math.sqrt(r)
    margin = 1.96 * stderr

    return {
        "mean"              : mean,
        "variance"          : var,
        "std"               : std,
        "stderr"            : stderr,
        "ci_lower"          : mean - margin,
        "ci_upper"          : mean + margin,
        "skewness"          : float(skew(ad_values)),
        "kurtosis"          : float(kurtosis(ad_values, fisher=True)),
    }


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    """
    Iterate over every (N, M, τ) combination, run R simulations, compute
    statistics, and write one CSV row per combination.

    Rows are sorted by N → M → τ (the natural iteration order below).
    """
    csv_columns = [
        "N", "M", "sampling_fraction", "t_mix",
        "tau", "T",
        "mean_AD", "variance", "standard_deviation", "standard_error",
        "CI_lower", "CI_upper",
        "skewness", "kurtosis",
    ]

    total = sum(len(M_VALUES[n]) for n in N_VALUES) * len(TAU_VALUES)
    done  = 0

    with open(OUTPUT_FILE, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
        writer.writeheader()

        for n in N_VALUES:
            t_mix = compute_mixing_time(n)

            for m in M_VALUES[n]:
                sampling_fraction = m / n

                for tau in TAU_VALUES:
                    # Whole-number transpositions for this τ
                    T = round(tau * t_mix)

                    # Progress indicator
                    done += 1
                    print(
                        f"  [{done:>4}/{total}]  N={n:>4}  M={m:>4}  "
                        f"τ={tau:.1f}  T={T:>6}  …",
                        end="\r",
                    )

                    # Run simulations
                    ad_values = run_all_simulations(n, m, T, R)

                    # Compute statistics
                    stats = compute_statistics(ad_values)

                    # Write one CSV row
                    writer.writerow({
                        "N"                  : n,
                        "M"                  : m,
                        "sampling_fraction"  : round(sampling_fraction, 6),
                        "t_mix"              : t_mix,
                        "tau"                : tau,
                        "T"                  : T,
                        "mean_AD"            : round(stats["mean"],     8),
                        "variance"           : round(stats["variance"],  8),
                        "standard_deviation" : round(stats["std"],       8),
                        "standard_error"     : round(stats["stderr"],    8),
                        "CI_lower"           : round(stats["ci_lower"],  8),
                        "CI_upper"           : round(stats["ci_upper"],  8),
                        "skewness"           : round(stats["skewness"],  8),
                        "kurtosis"           : round(stats["kurtosis"],  8),
                    })

    print(f"\n\n  Done. Results written to '{OUTPUT_FILE}'.")
    print(f"  Total rows written: {total}  (= {len(N_VALUES)} N-values "
          f"× 4 M-values × {len(TAU_VALUES)} τ-values)")


if __name__ == "__main__":
    main()
