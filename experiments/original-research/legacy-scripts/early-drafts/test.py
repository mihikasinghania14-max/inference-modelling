import random

def average_displacement(deck):
    n = len(deck)
    total = 0
    for i, card in enumerate(deck):
        total += abs(card - i)
    return total / n


def run_trial(n, tolerance=0.05):
    deck = list(range(n))

    target = (n**2 - 1)/(3*n)

    lower = target*(1-tolerance)
    upper = target*(1+tolerance)

    swaps = 0

    while True:

        i = random.randint(0,n-1)
        j = random.randint(0,n-1)

        deck[i], deck[j] = deck[j], deck[i]

        swaps += 1

        ad = average_displacement(deck)

        if lower <= ad <= upper:
            return swaps


def many_trials(n,trials=1000):

    results=[]

    for _ in range(trials):
        results.append(run_trial(n))

    mean=sum(results)/len(results)

    variance=sum((x-mean)**2 for x in results)/(len(results)-1)

    sd=variance**0.5

    return mean,sd,results

if __name__ == "__main__":
    for n in [50,100,200,400,800]:
        print(many_trials(n))


