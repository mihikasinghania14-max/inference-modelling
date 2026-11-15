N = 512        # deck size
M = 1
R = 1000      # number of independent Monte-Carlo simulations

# =============================================================================
# Imports
# =============================================================================

import math
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, kurtosis, skew


# =============================================================================
# Helper functions
# =============================================================================

def compute_mixing_time(n: int) -> int:
    """
    Return the theoretical RT-shuffle mixing time.

    The standard result (Diaconis & Shahshahani 1981) gives

        t_mix = ceil( n * ln(n) / 2 )

    Parameters
    ----------
    n : deck size
    """
    return math.ceil(0.5 * n * math.log(n))   # math.log uses the natural log


def run_single_simulation(n: int, m: int, t_mix: int) -> float:
    # --- Step 1: identity deck
    deck = np.arange(n, dtype=np.int32)

    # --- Step 2: apply t_mix random transpositions
    positions = np.random.randint(0, n, size=(t_mix, 2))
    for i, j in positions:
        deck[i], deck[j] = deck[j], deck[i]

    # --- Step 3: randomly sample m cards to track
    tracked_cards = np.random.choice(n, size=m, replace=False)
    original_indices = tracked_cards  # card i starts at position i

    # Build inverse permutation
    inv = np.empty(n, dtype=np.int32)
    inv[deck] = np.arange(n, dtype=np.int32)

    # Final positions of the tracked cards
    X = inv[tracked_cards]

    # --- Step 4: compute Average Displacement
    ad = np.mean(np.abs(X - original_indices))
    return float(ad)


def run_all_simulations(n: int, m: int, t_mix: int, r: int) -> np.ndarray:
    """
    Run r independent simulations and return an array of Average Displacement
    values, one per simulation.

    Parameters
    ----------
    n     : deck size
    m     : number of observed cards
    t_mix : mixing time (number of transpositions per simulation)
    r     : number of simulations
    """
    ad_values = np.empty(r, dtype=np.float64)
    for k in range(r):
        ad_values[k] = run_single_simulation(n, m, t_mix)
    return ad_values


# =============================================================================
# Statistics
# =============================================================================

def compute_statistics(ad_values: np.ndarray) -> dict:
    """
    Compute summary statistics for an array of simulation results.

    Returns a dictionary containing:
        mean, variance, std, stderr, ci_lower, ci_upper, cv
    where the 95 % CI uses the normal approximation  xbar ± 1.96 * s / sqrt(R).

    Parameters
    ----------
    ad_values : 1-D array of average displacement values
    """
    r       = len(ad_values)
    mean    = float(np.mean(ad_values))
    var     = float(np.var(ad_values, ddof=1))   # sample variance (unbiased)
    std     = math.sqrt(var)
    stderr  = std / math.sqrt(r)
    margin  = 1.96 * stderr

    return {
        "mean"    : mean,
        "variance": var,
        "std"     : std,
        "stderr"  : stderr,
        "ci_lower": mean - margin,
        "ci_upper": mean + margin,
        "cv"      : std / mean,   
        "skewness"  : float(skew(ad_values)),
        "kurtosis"  : float(kurtosis(ad_values, fisher=True)),
    }


# =============================================================================
# Output: printing
# =============================================================================

def print_statistics(stats: dict, n: int, m: int, r: int, t_mix: int) -> None:
    """
    Print a formatted summary of the experiment to the console.

    Parameters
    ----------
    stats  : dict returned by compute_statistics
    n, m, r: experimental parameters
    t_mix  : theoretical mixing time
    """
    sep = "=" * 58
    print(sep)
    print("  Partial Observation – Random Transposition Shuffle")
    print(sep)
    print(f"  Deck size              N       = {n}")
    print(f"  Observed cards         M       = {m}")
    print(f"  Simulations            R       = {r}")
    print(f"  Theoretical t_mix              = {t_mix}")
    print(sep)
    print(f"  Sample mean                    = {stats['mean']:.6f}")
    print(f"  Sample variance                = {stats['variance']:.6f}")
    print(f"  Sample std deviation           = {stats['std']:.6f}")
    print(f"  Standard error                 = {stats['stderr']:.6f}")
    print(f"  95 % CI lower bound            = {stats['ci_lower']:.6f}")
    print(f"  95 % CI upper bound            = {stats['ci_upper']:.6f}")
    print(sep)
    print(f"  Coefficient of variation (CV)  = {stats['cv']:.6f}")
    print(sep)

def plot_histogram(
    ad_values: np.ndarray,
    stats: dict,
    n: int,
    m: int,
    r: int,
    filepath: str = "histogram.png",
) -> None:
    """
    Plot the sampling distribution of Average Displacement with:
      – probability-density histogram (bin count chosen by Freedman–Diaconis)
      – fitted normal density overlay
      – vertical lines for the sample mean and 95 % CI limits
    and save to a PNG file.

    Parameters
    ----------
    ad_values : array of per-simulation AD values
    stats     : dict returned by compute_statistics
    n, m, r   : experimental parameters
    filepath  : output file path
    """
    mean      = stats["mean"]
    std       = stats["std"]
    ci_lower  = stats["ci_lower"]
    ci_upper  = stats["ci_upper"]

    # --- Publication-quality style -----------------------------------------
    plt.rcParams.update({
        "font.family"       : "serif",
        "font.size"         : 12,
        "axes.titlesize"    : 13,
        "axes.labelsize"    : 12,
        "legend.fontsize"   : 10,
        "figure.dpi"        : 150,
    })

    fig, ax = plt.subplots(figsize=(9, 5))

    # --- Histogram (density) with automatic bin count ----------------------
    # np.histogram_bin_edges with 'fd' uses the Freedman–Diaconis estimator,
    # which adapts to the data spread and is robust to outliers.
    bins = np.histogram_bin_edges(ad_values, bins="fd")
    ax.hist(
        ad_values,
        bins=bins,
        density=True,
        color="#4878CF",
        alpha=0.65,
        edgecolor="white",
        linewidth=0.5,
        label="Simulated AD (density)",
    )

    # --- Fitted normal density overlay -------------------------------------
    x_grid = np.linspace(ad_values.min(), ad_values.max(), 400)
    ax.plot(
        x_grid,
        norm.pdf(x_grid, loc=mean, scale=std),
        color="#E84646",
        linewidth=2.0,
        label=f"Fitted normal  (μ={mean:.3f}, σ={std:.3f})",
    )

    # --- Vertical reference lines ------------------------------------------
    ax.axvline(
        mean,
        color="#222222",
        linewidth=1.8,
        linestyle="-",
        label=f"Sample mean = {mean:.4f}",
    )
    ax.axvline(
        ci_lower,
        color="#F39C12",
        linewidth=1.4,
        linestyle="--",
        label=f"95 % CI  [{ci_lower:.4f},  {ci_upper:.4f}]",
    )
    ax.axvline(
        ci_upper,
        color="#F39C12",
        linewidth=1.4,
        linestyle="--",
    )
    ax.axvline(
        170.7,
        color="#39FF14",
        linewidth=1.8,
        linestyle="--",
        label="Reference AD = 170.7",
    )

    # --- Labels, title, legend ---------------------------------------------
    ax.set_xlabel("Average Displacement  $\\widehat{AD}_M$")
    ax.set_ylabel("Probability Density")
    ax.set_title(
        f"Sampling Distribution of Average Displacement "
        f"(N={n},  M={m},  R={r})"
    )
    ax.legend(framealpha=0.9, loc="upper right")

    # Light grid on the y-axis only for readability
    ax.yaxis.grid(True, linestyle=":", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)

    fig.tight_layout()
    plt.show()
    print(f"  Histogram saved           →  {filepath}")


# # =============================================================================
# # Main
# # =============================================================================

def main() -> None:
    """
    Orchestrates the full experiment pipeline:
      1. Validate parameters.
      2. Compute the theoretical mixing time.
      3. Run R independent simulations.
      4. Compute summary statistics.
      5. Print results and excel-ready output.
      6. Display histogram on screen.
    """
    # --- Parameter validation -----------------------------------------------
    if not (1 <= M <= N):
        raise ValueError(f"M must satisfy 1 ≤ M ≤ N, but got M={M}, N={N}.")
    if R < 2:
        raise ValueError(f"R must be at least 2 to compute sample variance, got R={R}.")

    # --- Theoretical mixing time --------------------------------------------
    t_mix = compute_mixing_time(N)
    print(f"\n  Theoretical mixing time  t_mix = ceil(N·ln(N)/2) = {t_mix}\n")

    # --- Monte-Carlo simulations --------------------------------------------
    print(f"  Running {R} simulations (N={N}, M={M}, t_mix={t_mix}) …")
    ad_values = run_all_simulations(N, M, t_mix, R)
    print("  Simulations complete.\n")

    # --- Statistics ---------------------------------------------------------
    stats = compute_statistics(ad_values)

    # --- Console output -----------------------------------------------------
    print_statistics(stats, N, M, R, t_mix)

    # --- Excel-ready two-column output --------------------------------------
# --- Excel-ready output -------------------------------------------------
    print("\n--- COPY INTO EXCEL ---")
    print(M)
    print(f"{stats['mean']}")
    print(f"{stats['variance']}")
    print(f"{stats['std']}")
    print(f"{stats['stderr']}")
    print(f"{stats['ci_lower']}")
    print(f"{stats['ci_upper']}")
    print(f"{stats['cv']}")
    print(f"{stats['skewness']}")
    print(f"{stats['kurtosis']}")
    print("-----------------------\n")

    # --- Show histogram on screen (no file saved) ---------------------------
    plot_histogram(ad_values, stats, N, M, R)
     


if __name__ == "__main__":
    main()


