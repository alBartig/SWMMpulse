import random

def get_weighted_population(pop, fraction):
    n = 0
    for i in range(pop):
        p = random.random()
        if p <= fraction:
            n += 1
    return n