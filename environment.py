import datetime as dt
import numpy as np
import random
import pandas as pd
from swmm_api import read_out_file


class UNITS:
    GRAM = 'g'
    MILLIGRAM = 'mg'
    MICROGRAM = 'ug'
    COUNT = '#'
    UNIT = "unit"
    DATE = "date"


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
    NAME = "name"
    FECAL = 'Fecal-Matter'
    COV = 'Cov-RNA'
    PEP = 'Pepper-virus'


class GROUP:
    GROUPS = "groups"
    NAME = "name"
    DAILYPOOPS = "dailypoops"
    WEIGHT = "weight"
    CONSTITUENTS = "constituents"
    PATTERN = "pattern"


class DEFAULT:
    CONST_FECAL = {CONSTITUENT.NAME: CONSTITUENT.FECAL, CONSTITUENT.SPECIFIC_LOAD: 200,
                   UNITS.UNIT: UNITS.GRAM, CONSTITUENT.DECAY_RATE: 0.114}
    CONST_COV = {CONSTITUENT.NAME: CONSTITUENT.COV, CONSTITUENT.SPECIFIC_LOAD: 1000, UNITS.UNIT: UNITS.COUNT}
    CONST_PEP = {CONSTITUENT.NAME: CONSTITUENT.PEP, CONSTITUENT.SPECIFIC_LOAD: 1000, UNITS.UNIT: UNITS.COUNT}

    DEFAULT_CONSTITUENTS = [CONST_FECAL, CONST_COV]

    PATTERN_BRISTOL_MEN = [1.4, 0.3, 0.1, 0.0, 0.3, 1.7, 9.1, 21, 13, 9, 6.9, 4.9, 1.9, 3.6, 2.5, 2, 2.9, 2.3, 4.1, 4.0,
                           2.7,
                           2.1, 2.2, 2.0]
    PATTERN_BRISTOL_WOMEN = [2.7, 0.1, 0.1, 0.08, 0.02, 0.2, 3.3, 16.9, 20.1, 12.4, 8.8, 5.0, 2.5, 2.6, 3.1, 2.5, 3.0,
                             2.2, 4.3,
                             3.3, 2.1, 1.5, 2.0, 1.2]
    PATTERN_BRISTOL = [round((w * 0.5 + v * 0.5), 3) for w, v in zip(PATTERN_BRISTOL_MEN, PATTERN_BRISTOL_WOMEN)]

    GROUP_HEALTHY = {GROUP.NAME: "Healthy", GROUP.WEIGHT: 0.8,
                     GROUP.CONSTITUENTS: [CONST_FECAL], GROUP.PATTERN: PATTERN_BRISTOL}
    GROUP_INFECTED = {GROUP.NAME: "Infected", GROUP.WEIGHT: 0.2,
                      GROUP.CONSTITUENTS: [CONST_FECAL, CONST_COV], GROUP.PATTERN: PATTERN_BRISTOL}

    ENVIRONMENT = {
        UNITS.DATE: dt.date(day=1, month=1, year=2000),
        GROUP.GROUPS: [GROUP_HEALTHY, GROUP_INFECTED],
        CONSTITUENT.DISPERSION_RATE: 0.16
    }


class Environment:
    def __init__(self, information=None):
        if information is None:
            information = DEFAULT.ENVIRONMENT
        self.information = information
        self.time_range = pd.date_range(information.get(UNITS.DATE), periods=24 * 60 * 6, freq="10S")
        # standardize patterns in groups
        for group in self.information.get(GROUP.GROUPS):
            group[GROUP.PATTERN] = self._standardize_pattern(group.get(GROUP.PATTERN))

    def _standardize_pattern(self, pattern):
        """
        Takes the input weights and interpolates weights for a frequency of 10S
        Args:
            pattern (list): list of weights

        Returns:
            list
        """
        date = self.information.get(UNITS.DATE)
        pattern.append(pattern[0])  # end pattern as it started
        dest_range = pd.date_range(date, date + dt.timedelta(days=1),
                                   periods=len(self.time_range) + 1,
                                   inclusive="both")  # date_range in destination frequency
        orig_range = pd.date_range(date, date + dt.timedelta(days=1),
                                   periods=len(pattern),
                                   inclusive="both")  # date_range with period of input pattern
        s_weights = pd.Series(pattern, index=orig_range)
        s_weights = s_weights.reindex(dest_range)
        s_weights.interpolate(limit_direction="both", inplace=True)
        s_weights.drop(s_weights.index[-1], inplace=True)
        return s_weights.values

    def _packet_gen(self, i, group, node=None):
        """
        creates a packet dictionary with an id and an origin time. The time is assigned randomly
        according to the time pattern in the environment
        Args:
            i (int): index number for packet, should be unique in the list

        Returns:
            dict
        """
        packet = {PACKET.PACKETID: f"P{int(i):d}",
                  PACKET.T0: random.choices(self.time_range.values, weights=group.get(GROUP.PATTERN), k=1),
                  PACKET.CLASSIFICATION: group.get(GROUP.NAME),
                  PACKET.CONSTITUENTS: group.get(GROUP.CONSTITUENTS),
                  PACKET.ORIGIN: node}
        return packet

    def get_packets(self, population):
        """
        Returns a dictionary with number of packets of the given environment.
        A pool of individuals is created by creating a list with the node-strings in it, as often as there are people at
        that node. This pool is shuffled. For each group, the weighted slice is taken from the pool and then
        extended/shrunken according to the dailypoops of that group. Then a dictionary with all packets and their
        metadata is created for each poop.
        Args:
            population (list): list of tuples (node, population)

        Returns:
            dict
        """
        # create each individual in environment as a string in order to pick randomly later
        people = []
        for node in population:
            people += node[0] * node[1]
        random.shuffle(people)  # shuffle list of people
        pop_tot = len(people)

        weights_tot = sum(group.get(GROUP.WEIGHT) for group in
                          self.information.get(GROUP.GROUPS))  # calculate sum of weights of groups
        packets = {}  # prepare dictionary for packets
        j = 0  # counter for packets
        for group in self.information.get(GROUP.GROUPS):
            n = int(np.around(group.get(GROUP.WEIGHT) * pop_tot / weights_tot, 0))  # calculate number of people to pick
            picks = people[j: j + n]  # slice the chunk of people
            dailypoops = group.get(GROUP.DAILYPOOPS)
            quotient, rest = (n * dailypoops / n), (
                        n * dailypoops / n)  # calculate modulo operator to multiply and slice to extend pickinglist
            # extend picks to make up for extra daily poops from non integer dailypoops
            picks = picks * quotient + picks[:rest]
            packets.update({i: self._packet_gen(i, node) for i, node in
                            zip(range(j, j + n), picks)})  # generate n packets and add to dictionary
            j += n  # counter for packets
        return packets

    def read_swmmoutfile(self, outfile_path):
        """
        Reads swmm-outfile, interpolates to 10s frequency and reindexes to chosen environment date
        Args:
            outfile_path (str): path to swmm-outfile

        Returns:
            None
        """
        with read_out_file(outfile_path) as out:
            # import outfile as dataframe
            df = out.to_frame()
            df = df.loc[:, ("link", slice(None), ["Flow_rate", "Flow_velocity"])].droplevel(0, axis=1)
            # reindex and interpolate dataframe
            _date = df.index[0].date()
            _dtindex = pd.date_range(_date, periods=8641, freq="10S")
            df = df.reindex(_dtindex)
            df.loc[df.index[-1], :] = df.loc[df.index[0], :]
            df.interpolate(direction="both", inplace=True)
            df.drop(df.index[-1], inplace=True)
            date = self.information.get(UNITS.DATE)
            df.index = df.index.map(lambda d: d.replace(year=date.year, month=date.month, day=date.day))
            self.flow_rates = df.xs("Flow-rate", axis=1, level=1).to_dict(orient="list")
            self.flow_velocities = df.xs("Flow-velocities", axis=1, level=1).to_dict(orient="list")
            self.TIMESTEPS = len(df.index)
        return None

    def _calc_index(self, time):
        hr, mn, sk = time.hour, time.minute, time.second  # calculate index of tm in timeseries
        i0 = (hr * 360 + mn * 6 + sk) % self.TIMESTEPS
        return i0

    def get_velocity(self, time, link):
        return self.flow_velocities.get(link)[self._calc_index(time)]

    def get_flow_rate(self, time, link):
        return self.flow_rates.get(link)[self._calc_index(time)]


def test_fractions():
    env = Environment()
    population = 1759
    fractions = [0.00030, 0.00010, 0.00005, 0.00003]
    test_arr = []
    for fraction in fractions:
        env.groups[0].weight = (1 - fraction)
        env.groups[1].weight = fraction
        temp = [[], []]
        for i in range(1000):
            plist = []
            [plist.extend(group.plist('sample_node', population)) for group in env.groups]
            n = len([p for p in plist if p.classification == 'Infected'])
            temp[0].append(n)
            temp[1].append(len(plist))
        print(f"{fraction}: Packets: min: {min(temp[1])}, max: {max(temp[1])}, Overall: {np.mean(temp[1])}")
    # print(plist)


def random_correlation():
    from collections import Counter
    env = Environment()
    population = 2000
    plist = []
    [plist.extend(group.plist('sample_node', population)) for group in env.groups]
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
    [plist.extend(group.plist("sample_node", population)) for group in env.groups]
    ots = [p.t0 for p in plist]
    c = Counter(ots)
    x, heights = list(c.keys()), list(c.values())
    plt.bar(x, heights, width=0.002)
    plt.show()

    print("finished")


def main():
    test_fractions()


if __name__ == '__main__':
    main()
