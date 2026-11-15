"""
Random Transposition Shuffle — Tracking Card 1
===============================================
Simulates 1000 shuffles of a 512-card deck and records
where card 1 ends up after t_mix random transpositions.
"""

import random
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Parameters ───────────────────────────────────────────────────────────────
N     = 512
T_MIX = round((N / 2) * math.log(N))   # natural log (ln)
R     = 1000
# ─────────────────────────────────────────────────────────────────────────────

print(f"# n={N}  t_mix={T_MIX}  R={R}")


def run_single_trial(n: int, t_mix: int) -> int:
    """
    Shuffle an ordered deck and return the final 1-based position of card 1.

    Uses a companion `positions` array so every lookup is O(1) —
    no list.index() calls needed.

    Returns an int in the range [1, n].
    """
    # Deck stored 0-based internally: card value k lives at index k initially
    deck      = list(range(n))   # deck[i]      = card at position i
    positions = list(range(n))   # positions[k] = current index of card k

    for _ in range(t_mix):
        # Pick two distinct positions in one call (no retry loop)
        i, j = random.sample(range(n), 2)

        # Update position lookup before the swap
        positions[deck[i]] = j
        positions[deck[j]] = i

        # Swap the cards
        deck[i], deck[j] = deck[j], deck[i]

    # Card 1 is stored internally as card 0 (0-based); convert to 1-based
    return positions[0] + 1


# ── Run all simulations ───────────────────────────────────────────────────────
final_positions = [run_single_trial(N, T_MIX) for _ in range(R)]


# ── Mean and standard deviation ───────────────────────────────────────────────
mean = sum(final_positions) / R
std  = (sum((x - mean) ** 2 for x in final_positions) / (R - 1)) ** 0.5

print(f"# mean position = {mean:.4f}")
print(f"# std deviation  = {std:.4f}")


# ── Binned frequencies (bins of width 16) ────────────────────────────────────
BIN_WIDTH = 16
num_bins  = N // BIN_WIDTH   # = 32 bins for n=512

# Count how many trials fell into each bin
bin_counts = [0] * num_bins
for pos in final_positions:
    bin_index = (pos - 1) // BIN_WIDTH   # pos 1-16 -> bin 0, 17-32 -> bin 1, ...
    bin_counts[bin_index] += 1

print(f"\n# Binned frequencies (bin width = {BIN_WIDTH})")
print(f"# {'Bin range':<12}  {'Frequency':>10}")
print(f"# {'-'*24}")
for b in range(num_bins):
    lo = b * BIN_WIDTH + 1
    hi = lo + BIN_WIDTH - 1
    print(f"# {lo:>3}–{hi:<6}     {bin_counts[b]:>6}")


# ── Plot: individual (width-1) and binned (width-16) histograms ───────────────
fig, (ax_ind, ax_bin) = plt.subplots(
    2, 1,
    figsize=(10, 7),
    sharex=True,                  # both plots share the same x-axis
    gridspec_kw={"hspace": 0.08}, # small gap between the two panels
)

# Expected uniform frequency per bin, shown as a reference line on both panels
uniform_individual = R / N                     # ≈ 1000/512 ≈ 1.95 per position
uniform_binned     = R / num_bins              # ≈ 1000/32  = 31.25 per bin

# --- Top panel: individual positions (bin width = 1) -------------------------
# Build per-position counts (index 1..N; index 0 unused)
pos_counts = [0] * (N + 1)
for p in final_positions:
    pos_counts[p] += 1

ax_ind.bar(
    range(1, N + 1),
    pos_counts[1:],
    width=1,
    color="#2563EB",
    alpha=0.7,
    label="Individual positions (bin width = 1)",
)
ax_ind.axhline(
    uniform_individual,
    color="#F97316", linewidth=1.5, linestyle="--",
    label=f"Uniform expected ({uniform_individual:.2f})",
)
ax_ind.set_ylabel("Frequency", fontsize=11)
ax_ind.set_title(
    f"Final position of card 1 after {T_MIX} random transpositions "
    f"(n={N}, R={R})",
    fontsize=12,
)
ax_ind.legend(fontsize=9, loc="upper right")
ax_ind.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

# --- Bottom panel: binned positions (bin width = 16) -------------------------
bin_centres = [b * BIN_WIDTH + BIN_WIDTH / 2 + 1 for b in range(num_bins)]

ax_bin.bar(
    bin_centres,
    bin_counts,
    width=BIN_WIDTH * 0.9,        # slight gap between bars for readability
    color="#7C3AED",
    alpha=0.7,
    label=f"Binned positions (bin width = {BIN_WIDTH})",
)
ax_bin.axhline(
    uniform_binned,
    color="#F97316", linewidth=1.5, linestyle="--",
    label=f"Uniform expected ({uniform_binned:.1f})",
)
ax_bin.set_xlabel("Position of card 1", fontsize=11)
ax_bin.set_ylabel("Frequency", fontsize=11)
ax_bin.legend(fontsize=9, loc="upper right")
ax_bin.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

# Shared x-axis ticks every 32 positions
ax_bin.set_xlim(0, N + 1)
ax_bin.xaxis.set_major_locator(ticker.MultipleLocator(32))

# Annotation box with mean and std
textstr = f"mean = {mean:.1f}\nstd  = {std:.1f}"
ax_bin.text(
    0.02, 0.92, textstr,
    transform=ax_bin.transAxes,
    fontsize=9,
    verticalalignment="top",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="grey", alpha=0.8),
)

fig.savefig("/mnt/user-data/outputs/card1_histogram.png", dpi=150, bbox_inches="tight")
print("\n# Plot saved to card1_histogram.png")
plt.show()
