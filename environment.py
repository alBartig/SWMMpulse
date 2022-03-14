from typing import Dict, Any

from datetime import datetime, timedelta
import numpy as np
from pobject import PObject
import random
from pconstants import Loading
from helpers import get_weighted_population
import pandas as pd
from swmm_api import read_out_file

class Units:
    GRAM = 'g'
    MILLIGRAM = 'mg'
    MICROGRAM = 'ug'
    COUNT = '#'


class PACKET:
    CLASSIFICATION = "classification"
    PACKETID = "pid"
    CONSTITUENTS = "constituents"
    ORIGIN = "origin"
    T0 = "t0"
    AGE = "age"
    ARRIVAL_TIME = "arrival_time"
    LINK = "link"
    NODE = "node"


class CONSTITUENT:
    DECAY_RATE = "decay_rate"
    SPECIFIC_LOAD = "specific_load"
    DISPERSION_RATE = "dispersion_rate"


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
    COV = Constituent(Loading.COV, 10000, Units.COUNT, degradation_coefficient=0.114)
    PEP = Constituent(Loading.PEP, 10000, Units.COUNT)

class DefaultPatterns:
    _BRISTOLmen = [1.4, 0.3, 0.1, 0.0, 0.3, 1.7, 9.1, 21, 13, 9, 6.9, 4.9, 1.9, 3.6, 2.5, 2, 2.9, 2.3, 4.1, 4.0, 2.7,
                   2.1, 2.2, 2.0]
    _BRISTOLwomen = [2.7, 0.1, 0.1, 0.08, 0.02, 0.2, 3.3, 16.9, 20.1, 12.4, 8.8, 5.0, 2.5, 2.6, 3.1, 2.5, 3.0, 2.2, 4.3,
                     3.3, 2.1, 1.5, 2.0, 1.2]
    BRISTOL = [round((w * 0.5 + v * 0.5), 3) for w, v in zip(_BRISTOLmen, _BRISTOLwomen)]


class Group:
    time_range = pd.date_range(datetime.datetime.today().date(), periods=24*60*6, freq="10S")

    def __init__(self, name, constituents=None, weight=0.5, pattern=DefaultPatterns.BRISTOL, dailypoops=1):
        self.name = name
        self.constituents = constituents
        self.weight = weight #as fraction of overall population
        self.pattern = self._standardize_pattern(pattern)
        self.dailypoops = dailypoops

    def __repr__(self):
        return  f'{self.name}, Constituents: {" ".join(self.constituents)}, Weight: {self.weight}, Pattern: {self.pattern}'

    def set(self, **kwargs):
        [setattr(self, k, v) for k, v in kwargs.items()]

    def _standardize_pattern(self, pattern):
        """
        Takes the input weights and interpolates weights for a frequency of 10S
        Args:
            pattern (list): list of weights

        Returns:
            list
        """
        dest_range = Group.time_range #date_range in destination frequency
        orig_range = pd.date_range(datetime.date.today(), datetime.date.today() + datetime.timedelta(days=1),
                                  periods=len(pattern)+1, inclusive="left") #date_range with period of input pattern
        s_weights = pd.Series(pattern, index=orig_range)
        s_weights = s_weights.reindex(dest_range)
        s_weights[s_weights.index[-1]] = s_weights[s_weights.index[0]]
        s_weights.interpolate(inplace=True)
        return s_weights.values

    def _packet_gen(self, i, node=None):
        """
        creates a packet dictionary with an id and an origin time. The time is assigned randomly
        according to the time pattern in the environment
        Args:
            i (int): index number for packet, should be unique in the list

        Returns:
            dict
        """
        packet = {PACKET.PACKETID: f"P{int(i):d}",
                  PACKET.T0: random.choices(Group.time_range.values, weights=self.pattern, k=1),
                  PACKET.CLASSIFICATION: self.name,
                  PACKET.CONSTITUENTS: self.constituents,
                  PACKET.ORIGIN: node}
        return packet

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
    INFECTED = Group(name='Infected', weight=0.2, constituents=[DefaultConstituents.FECAL,
                                                                DefaultConstituents.PEP, DefaultConstituents.COV])


class Environment:
    def __init__(self, groups=None, dispersion=0.16, date=datetime(day=1, month=1, year=2000)):
        if groups is None:
            groups = [DefaultGroups.HEALTHY, DefaultGroups.INFECTED]
        self.dispersion = dispersion
        self.groups = groups
        self.date = date

    def get_packets(self, population):
        """
        Returns a dictionary with number of packets of the given environment.
        A pool of individuals is created by creating a list with the node-strings in it, as often as there are people at
        that node. This pool is shuffled. For each group, the weighted slice is taken from the pool and then extended/shrunken
        according to the dailypoops of that group. Then a dictionary with all packets and their metadata is created for each poop.
        Args:
            population (list): list of tuples (node, population)

        Returns:
            dict
        """
        #create each individual in environment as a string in order to pick randomly later
        people = []
        for node in population:
            people += node[0] * node[1]
        random.shuffle(people) #shuffle list of people
        pop_tot = len(people)

        weights_tot = sum(g.weight for g in self.groups) #calculate sum of weights of groups
        packets = {} #prepare dictionary for packets
        j = 0 #counter for packets
        for group in self.groups:
            n = int(np.around(group.weight * pop_tot / weights_tot, 0)) #calculate number of people to pick
            picks = people[j, j+n] #slice the chunk of people
            d, r = (n * group.dailypoops / n), (n * group.dailypoops / n) #calculate modulo operator to multiply and slice to extend pickinglist
            picks = picks*d + picks[:r] #extend picks to make up for extra daily poops
            packets.update({i: group._packet_gen(i, node) for i, node in zip(range(j,j+n),picks)}) #generate n packets and add to dictionary
            j += n #counter for packets
        return packets

    @property
    def constituents(self):
        constituents = set()
        for group in self.groups:
            [constituents.add(c) for c in group.constituents]
        return list(constituents)

    def read_swmmoutfile(self, outfile_path):
        """
        Reads swmm-outfile, interpolates to 10s frequency and reindexes to chosen environment date
        Args:
            outfile_path (str): path to swmm-outfile

        Returns:
            None
        """
        with read_out_file(outfile_path) as out:
            #import outfile as dataframe
            df = out.to_frame()
            df = df.loc[:, ("link", slice(None), ["Flow_rate", "Flow_velocity"])].droplevel(0, axis=1)
            #reindex and interpolate dataframe
            _date = df.index[0].date()
            _dtindex = pd.date_range(_date, periods=8641, freq="10S")
            df = df.reindex(_dtindex)
            df.loc[df.index[-1], :] = df.loc[df.index[0], :]
            df.interpolate(direction="both", inplace=True)
            df.drop(df.index[-1], inplace=True)
            df.index = df.index.map(lambda d: d.replace(year=self.date.year, month=self.date.month, day=self.date.day))
            self.flow_rates = df.xs("Flow-rate", axis=1, level=1).to_dict(orient="list")
            self.flow_velocities = df.xs("Flow-velocities", axis=1, level=1).to_dict(orient="list")
            self.TIMESTEPS = len(df.index)
        return None

    def _calc_index(self, time):
        hr, mn, sk = time.hour, time.minute, time.second #calculate index of tm in timeseries
        i0 = (hr * 360 + mn * 6 + sk) % self.TIMESTEPS
        return i0

    def get_velocity(self, time, link):
        return self.flow_velocities.get(link)[self._calc_index(time)]

    def get_flow_rate(self, time, link):
        return self.flow_rates.get(link)[self._calc_index(time)]


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