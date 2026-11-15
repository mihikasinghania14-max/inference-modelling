import random

def random_transposition_shuffle(deck):
    n = len(deck)

    # Pick two different random positions
    position_1 = random.randint(0, n - 1)

    position_2 = position_1
    while position_2 == position_1:
        position_2 = random.randint(0, n - 1)

    # Swap the cards at those two positions
    temp_card = deck[position_1]
    deck[position_1] = deck[position_2]
    deck[position_2] = temp_card

    return deck


# --- Run it ---
n = 10
deck = list(range(1, n + 1))  # Build the deck once, outside the loop

print("Starting deck:", deck)

for i in range(100):
    deck = random_transposition_shuffle(deck)

print("Final deck:  ", deck)