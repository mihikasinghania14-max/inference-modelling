
import random
import math


N      = 50 # deck size
T_MIX  = 98 # number of random transpositions  (= ceil(0.5 * n * ln(n)))
R      = 1000    # number of trials (no need to change this)


def random_transposition_shuffle(deck: list, positions: list) -> None:
    i, j = random.sample(range(len(deck)), 2)

    # Keep the position lookup consistent before the swap
    positions[deck[i]] = j
    positions[deck[j]] = i

    deck[i], deck[j] = deck[j], deck[i]


def average_displacement(positions: list, n: int) -> float:
    total = 0
    for k in range(n):
        total += abs(positions[k] - k)
    return total / n


def run_single_trial(n: int, t_mix: int) -> float:
    deck      = list(range(n))   # card k starts at position k (0-based)
    positions = list(range(n))   # positions[k] = current index of card k

    for _ in range(t_mix):
        random_transposition_shuffle(deck, positions)

    return average_displacement(positions, n)


def main():
    print(f"# n={N}  t_mix={T_MIX}  R={R}")
    print(f"# theoretical_AD = {(N**2 - 1) / (3 * N):.4f}")
    print("# ---1000 AD values---")

    values = [run_single_trial(N, T_MIX) for _ in range(R)]

    for ad in values:
        print(f"{ad:.6f}")

    mean = sum(values) / R
    std  = (sum((x - mean) ** 2 for x in values) / (R - 1)) ** 0.5
    print(f"\n# mean={mean:.6f}  std={std:.6f}")

if __name__ == "__main__":
    main()
