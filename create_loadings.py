import numpy as np
import random
import matplotlib.pyplot as plt
import datetime
from pconstants import Distributions
from packet_router import Packet

class Group:
    def __init__(self, community, population, poops_per_person, classification, distribution=Distributions.BRISTOL, **kwargs):
        self.community = community
        self.population = population
        self.dailypoops = population * poops_per_person
        self.distribution = distribution
        self.classification = classification
        self.nmen = kwargs.get('nmen')
        if self.nmen == None:
            self.nmen = 686
        self.nwomen = kwargs.get('nwomen')
        if self.nwomen == None:
            self.nwomen = 888

        self.generate_weights()

    def generate_weights(self):
        '''
        kwargs:
        nmen - number of men, in case of "Bristol"-distribution
        nwomen - number of women, in case of "Bristol"-distribution
        '''
        weights = []
        if self.distribution == Distributions.BRISTOL:
            weightsmen = [1.4,0.3,0.1,0.0,0.3,1.7,9.1,21,13,9,6.9,4.9,1.9,3.6,2.5,2,2.9,2.3,4.1,4.0,2.7,2.1,2.2,2.0]
            weightswomen = [2.7,0.1,0.1,0.08,0.02,0.2,3.3,16.9,20.1,12.4,8.8,5.0,2.5,2.6,3.1,2.5,3.0,2.2,4.3,3.3,2.1,1.5,2.0,1.2]
            for w,v in zip(weightsmen,weightswomen):
                weights.append(round((w*self.nmen + v*self.nwomen)/(self.nmen+self.nwomen),3))
            self.weights = weights

    def generate_packets(self, **kwargs):
        '''
        kwargs:
        day - optional specification of day as datetime object
        '''
        nowd = kwargs.get('day')
        if nowd == None:
            nowd = datetime.datetime(day=2,month=3,year=2021)

        nt = 8640 #number of timeslots, a step every 10 seconds

        #prepare arrays
        index = [i for i in range(nt)]
        weights_ratio = round(nt / len(self.weights),2)
        weights = [self.weights[int(np.floor(i/weights_ratio))] for i in index]
        packets = []

        #prepare packets
        occs = random.choices(index, weights=weights, k=self.dailypoops)
        for occ in occs:
            h,r = divmod(occ,360)
            m,s = divmod(r,60)
            nowh = datetime.timedelta(hours=h,minutes=m,seconds=s)
            packets.append(Packet(self.community,self.classification,nowd+nowh))
        return packets

class Community:
    def __init__(self,name):
        self.name = name
        self.groups = {}
        self.packets = []

    def group_from_dict(self, dict):
        pop = dict.get('population')
        ppp = dict.get('ppp')
        cls = dict.get('classification')
        dis = dict.get('distribution')
        self.add_group(pop, ppp, cls, dis)

    def add_group(self, subcommunity):
        if type(subcommunity) == Group:
            self.groups[subcommunity.classification] = subcommunity
        if type(subcommunity) == list:
            for subcom in subcommunity:
                self.groups[subcom.classification] = subcom

    def community_timeseries(self):
        self.packets = []
        for subcom in self.groups.values():
            grouppackets = subcom.generate_packets()
            [self.packets.append(packet) for packet in grouppackets]
        return self.packets

if __name__ == '__main__':
    location = 'sample_city'
    sample_com = Community(location)
    group1 = Group(location,95,1,'healthy')
    group2 = Group(location,5,1,'infected')
    sample_com.add_group([group1,group2])
    packets = sample_com.community_timeseries()
    print('PACKETS IN SAMPLE_CITY')
    [print(packet) for packet in packets]