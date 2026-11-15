import random

def check_fp(deck, original):
    n = len(deck)
    fixed_points = 0
    for i in range(n): 
        if deck[i] == original[i]:
            fixed_points = fixed_points + 1
    
    return fixed_points

#def check_ad(deck):
    n = len(deck)
    dis = 0
    for i in range(0,n):
        if deck[i] != i + 1:
            dis = dis + abs(deck[i] - (i + 1))
    avg_dis = dis / n
    return avg_dis 

def check_inversions(deck):
    n = len(deck)
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if deck[i] > deck[j]:
                inversions = inversions + 1
    return inversions

# ONLY FOR AVG DISPLACEMENT- REDONE ON 10/06

def random_transposition_shuffle(deck):
    n = len(deck)

    position_1 = random.randint(0, n - 1)

    position_2 = position_1
    while position_2 == position_1:
        position_2 = random.randint(0, n - 1)

    temp_card = deck[position_1]
    deck[position_1] = deck[position_2]
    deck[position_2] = temp_card

    return deck

def check_ad(deck, original):
    n = len(deck)
    dis = 0
    for i in original:
        original_position = original.index(i)
        current_position = deck.index(i)
        dis = dis + abs(current_position - original_position)
    avg_dis = dis / n
    return avg_dis

def main():
    n = 4096
    deck = list(range(1, n + 1))  
    original = list(range(1, n + 1))

    print("Starting deck:", deck)

    for i in range(17035):
        deck = random_transposition_shuffle(deck)

    print("Final deck:  ", deck)
    print("Avg displacement:", check_ad(deck, original))
    print("Fixed points:", check_fp(deck, original))
    print("Inversions:", check_inversions(deck))


# def main():
#     n = 10
#     deck = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
#     original = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

#     print("Starting deck:", deck)

#     for i in range(100):
#         deck = random_transposition_shuffle(deck)

#     print("Final deck:  ", deck)
#     print("Avg displacement:", check_ad(deck, original))
#     print("Fixed points:", check_fp(deck, original))
#     print("Inversions:", check_inversions(deck))

if __name__ == "__main__":
    main()


 