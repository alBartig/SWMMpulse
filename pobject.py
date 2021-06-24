import pconstants

import numpy as np
import datetime
from pconstants import Discretization, Loading
from Exceptions import MassbalanceError
from timeseries import TDict

class PObject:
    def __init__(self, origin, classification, t0, contents):
        self.origin = origin
        self.classification = classification
        self.t0 = t0
        self.contents = contents

    def __repr__(self):
        return f'PObject[{self.origin}, {self.classification}, {self.t0}]\n'

    def __contains__(self, constituent):
        for c in self.contents:
            if c.name == constituent:
                return True
        return False

    def as_seconds(self, time):
        return time.hour*3600 + time.minute*60 + time.second

    def _calc_arriving_load(self, constituent):
        '''
        :param spec_load: Load contained in packet, g, mg, 1
        :param spec_degr: Degradation coefficient: 1/s
        :return: peak load [unit/m]
        '''
        spec_load = constituent.specific_load
        spec_degr = constituent.degradation_coefficient
        load_red = spec_load * self._factor_degradation(spec_degr)
        return load_red

    def _factor_dispersion_peak(self, dispersion):
        '''
        :param dispersion: Dispersion coefficient: m2/s
        :return: multiplier to reduce max load
        '''
        if self.age != 0:
            return 1/(2*(np.pi*dispersion*self.age)**0.5)
        else:
            return 1

    def _calculate_dispersion_spread(self, dispersion):
        '''
        :param dispersion: Dispersion coefficient: m2/s
        :return: spread in front of and behind expected mean, 2SD for 95% of load in m
        '''
        if self.age != None:
            return 2*(2*dispersion*self.age)**0.5
        else:
            print('packet does not have age information')

    def _factor_degradation(self, spec_degr):
        '''
        :param spec_degr: degradation coefficient: 1/s
        :return: multiplier to reduce max load
        '''
        if self.age != None:
            return np.e**(-spec_degr*self.age)
        else:
            print('packet does not have age information')

    def _peak_triangular(self, dist):
        """
        returns peak factor of a triangular distribution with a given spread
        Args:
            dist (float): spread distance (one side of the triangle)

        Returns:
            peak (float): peak factor
        """
        return 1/dist

    def interpolate(self, xg, x2, y2, x1=0, y1=0):
        """
        interpolation formula
        Args:
            xg (float):
            x2 (float):
            y2 (float):
            x1 (float):
            y1 (float):

        Returns:
            float
        """
        return y1 + (y2 - y1) * (xg - x1) / (x2 - x1)

    def _calculate_time_boundaries(self, dispersion, link, lookup):
        """
        Calculates first and last time parts of the load can be measured at node of interest.
        Args:
            dispersion (float): dispersion coefficient [m2/s]
            link (str): name of link upstream of node
        """
        #calculating the spread, start end end times of the dispersed packet
        spread_distance = self._calculate_dispersion_spread(dispersion)
        try:
            calc_spread = lookup.m_to_s(link,self.tm,spread_distance)
            spread_time =datetime.timedelta(seconds=max(calc_spread,6))
        except Exception:
            print('Unknown Error: Could not convert spread distance [m] to time [s]')
            exit()
        self.ta = self.tm - spread_time
        self.te = self.tm + spread_time

    def _interpolate(self, xg, x2, y2, x1=0, y1=0):
        return y1 + (y2 - y1) * (xg - x1) / (x2 - x1)

    def interpolate(self, h_max, t, constituent):
        """
        Returns an interpolated load/m for the given timestamp
        Args:
            t (datetime.datetime): time for which specific load is to be calculated
            constituent (str): name of constituent. Must be contained in packet
        Returns:
            float: specific load at time t
        """
        if self.ta < t <= self.tm:
            return (t-self.ta).seconds*h_max/(self.tm-self.ta).seconds
        elif self.tm < t < self.te:
            return (self.te-t).seconds*h_max/(self.te-self.tm).seconds
        else:
            return 0.0

    def prepare_timesteps(self, lookup):
        """
        prepares an array with all timesteps for a load-distribution
        Args:
            ta (datetime.datetime): timestamp of beginning of load arrival
            te (datetime.datetime): timestamp of end of load
        Returns:
            timesteps (list)
        """
        first,shift_front = lookup.find_smaller_dt(self.ta)
        last,shift_back = lookup.find_smaller_dt(self.te)
        diff = (last + shift_back) - (first + shift_front)
        size = int(diff/Discretization.TIMESTEPLENGTH + 1)
        return [first+step*Discretization.TIMESTEPLENGTH for step in range(size)]

    def _prepare_loadheights(self, h_max, timesteps, constituent):
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
            heights.append(self.interpolate(h_max, step,constituent))
        return heights

    def _calculate_trapezoid(self, h1, h2, dt=Discretization.TIMESTEPLENGTH):
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
        if type(dt) is not float:
            dt = dt.total_seconds()
        trapezoid = 0.5 * (h1+h2) * dt
        if type(trapezoid) is not float:
            print('Implausible load calculated')
            print(f'Load: {trapezoid}')
            raise BaseException
        return trapezoid

    def _integrate_timestep_loads(self, heights, timesteps):
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
            try:
                if i == 0:
                    #first load
                    loads.append(self._calculate_trapezoid(0, heights[i+1], timesteps[i+1] - self.ta))
                elif i == len(timesteps)-1:
                    #last load
                    loads.append(self._calculate_trapezoid(heights[i], 0, self.te - timesteps[i]))
                else:
                    loads.append(self._calculate_trapezoid(heights[i], heights[i + 1]))
            except:
                print(f'Error calculating trapezoids: at index {i}')
                print(f'Heights ({len(heights)}): {heights}')
                print(f'Timesteps ({len(timesteps)}): {timesteps}')
                raise Exception
        return loads

    def _check_mass_balance(self, loads, constituent, spec_degr):
        expected_total = self.max_loads[constituent] * self._factor_degradation(spec_degr)
        if abs(1 - sum(loads)/expected_total) > 0.01:
            return False
        else:
            return True

    def get_loads(self, constituent, link, qlookup, tlookup):
        """
        Calculates the load distribution of given constituent and returns loads per timestep and timesteps
        Args:
            constituent (str): name of constituent, must be contained in Packet.contents
        Returns:
            (tuple): Tuple of lists (loads,timesteps)
        """
        mred = self.max_loads[constituent]
        #calculate dispersion spread
        dist = self._calculate_dispersion_spread(self.dispersion)
        #looking up flow velocity and calculate spatial resolution corresponding to timestep and flowvelocity
        v = qlookup.lookup_v(link,self.tm)[0]
        dx = round(v,1) * Discretization.TIMESTEPLENGTH.total_seconds() / Discretization.REFINE
        #creating linear spatial array with distance to load peak
        d = np.arange(-2*dist,2*dist,dx)
        #calculate peak loading [unit/m]
        peak_m = self._peak_triangular(dist) * mred
        #interpolate loading for each section
        cx = [max(0, self._interpolate(abs(x), dist, 0, y1=peak_m)) for x in d]
        #integrate loading for each section and aggregate to timesteps
        mx = [c*dx for c in cx]
        mt = [sum(mx[i:(i+Discretization.REFINE)]) for i in range(0,len(mx),Discretization.REFINE)]
        #convert spatial array to temporal arry
        ref_time = tlookup.find_closest(self.tm)
        t = [ref_time + datetime.timedelta(seconds=round(x/v)) for x in d[::Discretization.REFINE]]
        #check mass-balance
        totalm = sum(mt)
        deviation = abs(totalm-mred)/mred
        if not deviation <= Loading.MAX_CONTINUITY_ERROR:
            print(f"Expected load: {mred}\n"
                  f"Calculated load: {totalm}\n"
                  f"Continuity Error: {deviation*100:.2f} %\n"
                  f"Allowerd Error: {Loading.MAX_CONTINUITY_ERROR*100:.2} %")
            raise MassbalanceError
        return {'values':mt,'timestamps':t,'origin':self.origin,'traveltime':self.age,'class':self.classification}

    def set_lookup(self, qlookup):
        try:
            self.qlookup = qlookup
            return True
        except Exception:
            print('Unknown Error, could not assign lookup table')
            return False

    def set_arrival(self, tm):
        """
        Assigns arrival-time at node of interest for packet. Calculates subsequent age and derived packet properties
        like the dispersion, degradation and maximum loads of constituents
        Args:
            tm (datetime.datetime): time of arrival at node of interest
            dispersion (float): dispersion coefficient [m2/s]
        """
        self.tm = tm
        try:
            td = tm - self.t0
            self.age = td.total_seconds()
        except:
            raise ValueError
        self.max_loads = {}
        for constituent in self.contents:
            self.max_loads[constituent.name] = self._calc_arriving_load(constituent)
        return self

    def set_dispersion(self, dispersion):
        self.dispersion = dispersion

def main():
    from environment import Environment
    import datetime as dt
    from pconstants import Loading
    from application import prepare_environment
    env = Environment()

    t0 = dt.datetime.now()
    tm = t0 + dt.timedelta(seconds=30*60)
    evalnode = 'MH327-088-1'
    graph, qlut = prepare_environment()
    link = graph.get_outletlinks(evalnode)[0]

    start, end = qlut.timestamps[0], qlut.timestamps[0] + datetime.timedelta(days=1)
    timestamps = [start + i * Discretization.TIMESTEPLENGTH for i in
                  range(int((end - start) / Discretization.TIMESTEPLENGTH))]
    tlut = TDict(timestamps)

    x = PObject("Testorigin", "Infected", t0=t0, contents=env.groups[1].constituents)
    x.set_arrival(tm)
    x.set_dispersion(env.dispersion)
    loads = x.get_loads(Loading.FECAL,link,qlut, tlut)
    print(loads)
    print("main finished")

if __name__ == "__main__":
    main()
    print("Module finished")