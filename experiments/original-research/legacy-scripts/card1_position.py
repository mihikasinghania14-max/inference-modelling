import random
import math

N     = 1024
T_MIX = round((N / 2) * math.log(N))   
R     = 1000

print(f"# n={N}  t_mix={T_MIX}  R={R}")


def run_single_trial(n: int, t_mix: int) -> int:
    deck      = list(range(n))    
    positions = list(range(n))   

    for _ in range(t_mix):
        
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
