import numpy as np
import datetime
from pconstants import Discretization

class PObject:
    def __init__(self, origin, classification, t0, contents):
        self.origin = origin
        self.classification = classification
        self.t0 = t0
        self.contents = contents

    def __repr__(self):
        return f'[{self.origin}, {self.classification}, {self.t0}]\n'

    def as_seconds(self, time):
        return time.hour*3600 + time.minute*60 + time.second

    def set_lookup(self, qlookup):
        try:
            self.qlookup = qlookup
            return True
        except Exception:
            print('Unknown Error, could not assign lookup table')
            return False

    def set_arrival(self, tm, dispersion):
        """
        Assigns arrival-time at node of interest for packet. Calculates subsequent age and derived packet properties
        like the dispersion, degradation and maximum loads of constituents
        Args:
            tm (datetime.datetime): time of arrival at node of interest
            dispersion (float): dispersion coefficient [m2/s]
        """
        self.tm = tm
        self.age = tm - self.t0
        self.max_loads = {}
        for constituent,properties in self.contents.items():
            self.max_loads[constituent] = self.calculate_maximum_load(properties, dispersion)

    def calculate_maximum_load(self, properties, dispersion):
        '''
        :param spec_load: Load contained in packet, g, mg, 1
        :param spec_degr: Degradation coefficient: 1/s
        :param dispersion: Dispersion coefficient: m2/s
        :return: tuple (peak load [unit/m], spread [+- m])
        '''
        spec_load = properties['spec_load']
        spec_degr = properties['spec_degr']
        load_max = spec_load * self.calculate_degradation(spec_degr) * self.calculate_dispersion(dispersion)
        return load_max

    def calculate_dispersion(self, dispersion):
        '''
        :param dispersion: Dispersion coefficient: m2/s
        :return: multiplier to reduce max load
        '''
        if self.age != None:
            return 1/(2*(np.pi*dispersion*self.age)**0.5)
        else:
            print('packet does not have age information')

    def calculate_spread(self, dispersion):
        '''
        :param dispersion: Dispersion coefficient: m2/s
        :return: spread in front of and behind expected mean, 2SD for 95% of load in m
        '''
        if self.age != None:
            return 2*(2*dispersion*self.age)**0.5
        else:
            print('packet does not have age information')

    def calculate_degradation(self, spec_degr):
        '''
        :param spec_degr: degradation coefficient: 1/s
        :return: multiplier to reduce max load
        '''
        if self.age != None:
            return np.e**(-spec_degr*self.age)
        else:
            print('packet does not have age information')

    def calculate_time_boundaries(self, dispersion,link):
        """
        Calculates first and last time parts of the load can be measured at node of interest.
        Args:
            dispersion (float): dispersion coefficient [m2/s]
            link (str): name of link upstream of node
        """
        #calculating the spread, start end end times of the dispersed packet
        spread_distance = self.calculate_spread(dispersion)
        try:
            spread_time = datetime.timedelta(self.qlookup.m_to_s(link,time=self.tm))
        except Exception:
            print('Unknown Error: Could not convert spread distance [m] to time [s]')
            exit()
        self.ta = self.ta - datetime.timedelta(seconds=spread_time)
        self.te = self.ta + datetime.timedelta(seconds=spread_time)

    def interpolate(self, t, constituent):
        """
        Returns an interpolated load/m for the given timestamp
        Args:
            t (datetime.datetime): time for which specific load is to be calculated
            constituent (str): name of constituent. Must be contained in packet
        Returns:
            float: specific load at time t
        """
        h = self.max_loads[constituent]
        if self.ts < t <= self.ta:
            return (t-self.ts).seconds*h/(self.ta-self.ts).seconds
        elif self.ta < t < self.te:
            return (self.te-t).seconds*h/(self.te-self.ta).seconds
        else:
            return 0.0

    @property
    def prepare_timesteps(self):
        """
        prepares an array with all timesteps for a load-distribution
        Args:
            ta (datetime.datetime): timestamp of beginning of load arrival
            te (datetime.datetime): timestamp of end of load
        Returns:
            timesteps (list)
        """
        first = self.tlookup.find_smaller(self.ta)
        last = self.tlookup.find_larger(self.te)
        diff = last - first
        size = int(diff.seconds/Discretization.TIMESTEPLENGTH + 1)
        return [first+datetime.timedelta(seconds=step*Discretization.TIMESTEPLENGTH) for step in range(size)]

    def prepare_loadheights(self, timesteps, constituent):
        """
        Prepares an array with a corresponding specific load for each timestep in input array and constituent
        Args:
            timesteps (list): array with timesteps (datetime.datetime)
            constituent (str): name of the constituent
        Returns:
            list: list of specific loads, starts and ends with 0
        """
        heights = []
        for step in timesteps:
            heights.append(self.interpolate(step,constituent))
        return heights

    def calculate_trapezoid(self, h1, h2, dt=Discretization.TIMESTEPLENGTH):
        """
        integrates the values h1 and h2 across timesteplength dt, which is assumed
         to be the assigned TIMESTEPLENGTH if no other value is given.
         A trapezoidal shape is assumed.
        Args:
            h1 (float): first value
            h2 (float): last value
            dt (integer): timestep
        Returns:
            float
        """
        if type(dt) == datetime.timedelta:
            dt = dt.seconds
        return 0.5 * (h1+h2) * dt

    def integrate_timestep_loads(self, heights,timesteps):
        """
        integrates the load for each timestep over the timesteplength
        Args:
            heights (list): specific load at the beginning of each timestep
            timesteps (list): list of timesteps
        Returns:
            list
        """
        loads = []
        for i,time in enumerate(timesteps):
            if i == 0:
                #first load
                loads.append(self.calculate_trapezoid(0,heights[1],timesteps[0]-self.ta))
            elif i == len(timesteps-1):
                #last load
                loads.append(self.calculate_trapezoid(0,heights[i],self.te-timesteps[i]))
            else:
                loads.append(self.calculate_trapezoid(heights[i],heights[i+1]))
        return loads

    def get_loads(self, constituent):
        """
        Calculates the load distribution of given constituent and returns loads per timestep and timesteps
        Args:
            constituent (str): name of constituent, must be contained in Packet.contents
        Returns:
            (tuple): Tuple of lists (loads,timesteps)
        """
        timesteps = self.prepare_timesteps()
        heights = self.prepare_loadheights(timesteps,constituent)
        loads = self.integrate_timestep_loads(heights,timesteps)
        return {'values':loads,'timestamps':timesteps,'origin':self.origin,'traveltime':self.age,'class':self.classification}