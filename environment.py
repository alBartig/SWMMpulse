from datetime import datetime, timedelta
import numpy as np
from pobject import PObject
import random

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
    FECAL = Constituent('Fecal Matter', 200, Units.GRAM)
    COV = Constituent('Cov RNA', 10000, Units.COUNT)
    PEP = Constituent('Pepper virus', 10000, Units.COUNT)

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


    def plist(self, origin, population, nowd=None):
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
        npoops = int(round(population*self.weight*self.dailypoops))
        occs = random.choices(index, weights=pattern, k=npoops)
        for occ in occs:
            h,r = divmod(occ,360)
            m,s = divmod(r,60)
            nowh = timedelta(hours=h,minutes=m,seconds=s)
            packetlist.append(PObject(origin, self.name, nowd+nowh, self.constituents))
        return packetlist


class DefaultGroups:
    HEALTHY = Group(name='Healthy', weight=0.8, constituents=[DefaultConstituents.FECAL, DefaultConstituents.PEP])
    INFECTED = Group(name='Infected', weight=0.2, constituents=[DefaultConstituents.FECAL,\
                                                                DefaultConstituents.PEP, DefaultConstituents.COV])


class Environment:
    def __init__(self, groups=[DefaultGroups.HEALTHY, DefaultGroups.INFECTED], dispersion=1.6, date=datetime(day=2,month=3,year=2021)):
        self.dispersion = dispersion
        self.groups = groups
        self.date = date

if __name__ == '__main__':
    env = Environment()
    population = 100
    plist = []
    [plist.extend(group.plist('sample_node',population)) for group in env.groups]
    print(plist)