from datetime import datetime, timedelta
import numpy as np
from pobject import PObject
import random
from pconstants import Loading
from helpers import get_weighted_population

class Units:
    GRAM = 'g'
    MILLIGRAM = 'mg'
    MICROGRAM = 'ug'
    COUNT = '#'


class Constituent:
    def __init__(self, name, specific_load, unit, degradation_coefficient=0, probability=1):
        self.name = name
        self.specific_load = specific_load
        self.unit = unit
        self.degradation_coefficient = degradation_coefficient
        self.probability = probability

    def __repr__(self):
        return f'{self.name}'


class DefaultConstituents:
    FECAL = Constituent(Loading.FECAL, 200, Units.GRAM)
    COV = Constituent(Loading.COV, 10000, Units.COUNT)
    PEP = Constituent(Loading.PEP, 10000, Units.COUNT)

class DefaultPatterns:
    _BRISTOLmen = [1.4, 0.3, 0.1, 0.0, 0.3, 1.7, 9.1, 21, 13, 9, 6.9, 4.9, 1.9, 3.6, 2.5, 2, 2.9, 2.3, 4.1, 4.0, 2.7,
                   2.1, 2.2, 2.0]
    _BRISTOLwomen = [2.7, 0.1, 0.1, 0.08, 0.02, 0.2, 3.3, 16.9, 20.1, 12.4, 8.8, 5.0, 2.5, 2.6, 3.1, 2.5, 3.0, 2.2, 4.3,
                     3.3, 2.1, 1.5, 2.0, 1.2]
    BRISTOL = [round((w * 0.5 + v * 0.5), 3) for w, v in zip(_BRISTOLmen, _BRISTOLwomen)]


class Group:
    def __init__(self, name, constituents=None, weight=0.5, pattern=DefaultPatterns.BRISTOL, dailypoops=1):
        self.name = name
        self.constituents = constituents
        self.weight = weight #as fraction of overall population
        self.pattern = pattern
        self.dailypoops = dailypoops

    def __repr__(self):
        return  f'{self.name}, Constituents: {" ".join(self.constituents)}, Weight: {self.weight}, Pattern: {self.pattern}'

    def set(self, **kwargs):
        [setattr(self, k, v) for k, v in kwargs.items()]


    def plist(self, origin, population, nowd=None, startid=0):
        """
        Creates packet list
        Args:
            origin (str): name of node of origin
            population (int): number of overall population at node
            nowd (datetime): date

        Returns:
            list: list of packets
        """
        if nowd is None:
            nowd = datetime(day=2,month=3,year=2021)

        nt = 8640 #number of timeslots, a step every 10 seconds

        #prepare arrays
        index = [i for i in range(nt)]
        pattern_ratio = round(nt / len(self.pattern),2)
        pattern = [self.pattern[int(np.floor(i/pattern_ratio))] for i in index]
        packetlist = []

        #prepare packets
        weighted_population = get_weighted_population(population, self.weight)
        npoops = int(round(weighted_population*self.dailypoops))
        occs = random.choices(index, weights=pattern, k=npoops)
        for i, occ in enumerate(occs, start=startid):
            h,r = divmod(occ,360)
            m,s = divmod(r,6)
            nowh = timedelta(hours=h,minutes=m,seconds=10*s)
            packetlist.append(PObject(origin, self.name, nowd+nowh, self.constituents, pid=f"P{i}"))
        return packetlist


class DefaultGroups:
    HEALTHY = Group(name='Healthy', weight=0.8, constituents=[DefaultConstituents.FECAL, DefaultConstituents.PEP])
    INFECTED = Group(name='Infected', weight=0.2, constituents=[DefaultConstituents.FECAL,\
                                                                DefaultConstituents.PEP, DefaultConstituents.COV])


class Environment:
    def __init__(self, groups=[DefaultGroups.HEALTHY, DefaultGroups.INFECTED], dispersion=0.16, date=datetime(day=17,month=8,year=2020)):
        self.dispersion = dispersion
        self.groups = groups
        self.date = date

    @property
    def constituents(self):
        constituents = set()
        for group in self.groups:
            [constituents.add(c) for c in group.constituents]
        return list(constituents)

def test_fractions():
    env = Environment()
    population = 1759
    fractions = [0.00030,0.00010,0.00005,0.00003]
    test_arr = []
    for fraction in fractions:
        env.groups[0].weight = (1 - fraction)
        env.groups[1].weight = fraction
        temp = [[],[]]
        for i in range(1000):
            plist = []
            [plist.extend(group.plist('sample_node',population)) for group in env.groups]
            n = len([p for p in plist if p.classification == 'Infected'])
            temp[0].append(n)
            temp[1].append(len(plist))
        print(f"{fraction}: Packets: min: {min(temp[1])}, max: {max(temp[1])}, Overall: {np.mean(temp[1])}")
    #print(plist)

def random_correlation():
    from collections import Counter
    env = Environment()
    population = 2000
    plist = []
    [plist.extend(group.plist('sample_node',population)) for group in env.groups]
    originminutes = [p.t0.minute for p in plist]
    c = Counter(originminutes)
    print(f"Minuten: {list(c.keys())}")
    print("Function finished")

def random_pattern():
    from collections import Counter
    import matplotlib.pyplot as plt

    env = Environment()
    population = 20000
    plist = []
    [plist.extend(group.plist("sample_node",population)) for group in env.groups]
    ots = [p.t0 for p in plist]
    c = Counter(ots)
    x, heights = list(c.keys()), list(c.values())
    plt.bar(x, heights, width = 0.002)
    plt.show()

    print("finished")

def main():
    test_fractions()

if __name__ == '__main__':
    main()