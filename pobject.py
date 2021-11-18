import pconstants
import pandas as pd
import numpy as np
import datetime
from pconstants import Discretization, Loading
from Exceptions import MassbalanceError
from timeseries import TDict


class PObject:
    def __init__(self, origin, classification, t0, contents, pid=None):
        self.origin = origin
        self.classification = classification
        self.t0 = t0
        self.contents = contents
        self.pid = pid

    def __repr__(self):
        return f'PObject[{self.origin}, {self.classification}, {self.t0}]\n'

    def __contains__(self, constituent):
        for c in self.contents:
            if c.name == constituent:
                return True
        return False

    def as_seconds(self, time):
        return time.hour * 3600 + time.minute * 60 + time.second

    def to_list(self):
        return [self.pid, self.origin, self.classification, self.t0, [c.name for c in self.contents]]

    @property
    def tags(self):
        return ["pid", "origin_location", "classification", "origin_time", "contents"]

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
            return 1 / (2 * (np.pi * dispersion * self.age) ** 0.5)
        else:
            return 1

    def _calculate_dispersion_spread(self, dispersion):
        '''
        returns 2SD of normal distributed dispersion after Fick
        SD = (2*D*t)**0.5
        95% of load are within 2SD of the mean
        Args:
            dispersion (float): Dispersion coefficient: m2/s

        Returns:
            2SD (float): 2SD of normal distributed dispersion
        '''
        if self.age != None:
            return 2 * (2 * dispersion * self.age) ** 0.5
        else:
            print('packet does not have age information')

    def _factor_degradation(self, spec_degr):
        '''
        :param spec_degr: degradation coefficient: 1/d
        :return: multiplier to reduce max load
        '''
        if self.age != None:
            return np.e ** (-spec_degr * self.age / 86400)
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
        try:
            return 1 / dist
        except:
            print(f"Trying to divide 1 by {dist}, Type: {type(dist)}")
            raise ValueError

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
        # calculating the spread, start end end times of the dispersed packet
        spread_distance = self._calculate_dispersion_spread(dispersion)
        try:
            calc_spread = lookup.m_to_s(link, self.tm, spread_distance)
            spread_time = datetime.timedelta(seconds=max(calc_spread, 6))
        except Exception:
            print('Unknown Error: Could not convert spread distance [m] to time [s]')
            raise BaseException
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
            return (t - self.ta).seconds * h_max / (self.tm - self.ta).seconds
        elif self.tm < t < self.te:
            return (self.te - t).seconds * h_max / (self.te - self.tm).seconds
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
        first, shift_front = lookup.find_smaller_dt(self.ta)
        last, shift_back = lookup.find_smaller_dt(self.te)
        diff = (last + shift_back) - (first + shift_front)
        size = int(diff / Discretization.TIMESTEPLENGTH + 1)
        return [first + step * Discretization.TIMESTEPLENGTH for step in range(size)]

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
            heights.append(self.interpolate(h_max, step, constituent))
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
        trapezoid = 0.5 * (h1 + h2) * dt
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
        for i, time in enumerate(timesteps):
            try:
                if i == 0:
                    # first load
                    loads.append(self._calculate_trapezoid(0, heights[i + 1], timesteps[i + 1] - self.ta))
                elif i == len(timesteps) - 1:
                    # last load
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
        if abs(1 - sum(loads) / expected_total) > 0.01:
            return False
        else:
            return True

    def calc_load(self, constituent, link, qlookup):
        """
        Calculates the load distribution of given constituent and returns loads per timestep and timesteps
        Args:
            constituent (str): name of constituent, must be contained in Packet.contents
        Returns:
            (tuple): Tuple of lists (loads,timesteps)
        """
        reduced_load = self.max_loads[constituent]
        # calculate dispersion spread
        dispersion_spread = self._calculate_dispersion_spread(self.dispersion)
        # looking up flow velocity and calculate spatial resolution corresponding to timestep and flowvelocity
        flow_velocity = round(qlookup.lookup_v(link, self.tm)[0], 1)
        dx = flow_velocity * Discretization.TIMESTEPLENGTH.total_seconds()
        # creating linear spatial array with distance to load peak, making sure x=0.0 is contained in d
        borders = (2 * np.ceil(dispersion_spread / dx) - 1) * dx
        dx = dx / Discretization.REFINE
        d = np.arange(-borders, borders, dx)
        # check whether pulse load is spread-out over more than one dx / catching div by 0
        if len(d) == 0:
            spatial_loading_aggregated = [reduced_load / dx]
            d = [0]
        else:
            # calculate peak loading [unit/m]
            peak_loading = self._peak_triangular(dispersion_spread) * reduced_load
            # interpolate loading for each section
            spatial_linear_loading = [max(0, self._interpolate(abs(x), dispersion_spread, 0, y1=peak_loading)) for x in
                                      d]
            # integrate loading for each section and aggregate to timesteps
            spatial_loading_fine = [c * dx for c in spatial_linear_loading]
            spatial_loading_aggregated = [sum(spatial_loading_fine[i:(i + Discretization.REFINE)]) for i in
                                          range(0, len(spatial_loading_fine), Discretization.REFINE)]
        # convert spatial array to temporal array
        ref_time = qlookup.find_seconds(self.tm)
        values = np.zeros(len(qlookup.timestamps))
        mid = qlookup.timestamps.index(qlookup.norm_time(ref_time))
        indexer = [round(di / (dx * Discretization.REFINE)) + mid for di in d[::Discretization.REFINE]]
        try:
            values[indexer] = spatial_loading_aggregated
        except:
            pass
        # t = [ref_time + datetime.timedelta(seconds=round(x/flow_velocity)) for x in d[::Discretization.REFINE]]
        # check mass-balance
        totalm = sum(spatial_loading_aggregated)
        deviation = abs(totalm - reduced_load) / reduced_load
        if not deviation <= Loading.MAX_CONTINUITY_ERROR:
            print(f"Expected load: {reduced_load}\n"
                  f"Calculated load: {totalm}\n"
                  f"Continuity Error: {deviation * 100:.2f} %\n"
                  f"Allowed Error: {Loading.MAX_CONTINUITY_ERROR * 100:.2} %")
            raise MassbalanceError
        return values

    def spread(self, dispersion_spread=2, lred=1, f=1, res=1):
        d = dispersion_spread
        x = np.arange(-f * d, d + res, res)
        peak = 2 * lred / ((f + 1) * d)
        y = np.interp(x, [-f * d, 0, d], [0, peak, 0])

        return x, y

    def superimpose(self, series, res=1):
        xmin = min([s[0][0] for s in series])
        xmax = max([s[0][-1] for s in series])

        xcomb = np.round(np.arange(xmin, xmax + res, res), 8)
        ycomb = np.zeros(len(xcomb))

        for s in series:
            ycomb[np.nonzero(np.isin(xcomb, np.round(s[0], 8)))] += s[1]

        return xcomb, ycomb

    def calc_composite_load(self, constituent, link, qlookup):
        """
        Calculates the load distribution of given constituent and returns loads per timestep and timesteps
        Args:
            constituent (str): name of constituent, must be contained in Packet.contents
        Returns:
            (tuple): Tuple of lists (loads,timesteps)
        """
        reduced_load = self.max_loads[constituent]
        fractions = Loading.FRACTIONS[constituent]
        ref_time = qlookup.norm_time(self.tm)

        series = []
        for fraction in fractions:
            # calculate dispersion spread
            dispersion_spread = self._calculate_dispersion_spread(fraction.get(Loading.DISPERSION))
            # looking up flow velocity and calculate spatial resolution corresponding to timestep and flowvelocity
            flow_velocity = round(qlookup.lookup_v(link, ref_time)[0], 1)
            dx = flow_velocity * Discretization.TIMESTEPLENGTH.total_seconds()
            # creating linear spatial array with distance to load peak, making sure x=0.0 is contained in d
            # set x-array borders to be a multiple of dx, x-max is therefore usually a bit larger than dispersion spread
            border = np.ceil(dispersion_spread / dx) * dx
            # dx = dx / Discretization.REFINE
            # check whether pulse load is spread-out over more than one dx / catching div by 0
            fraction_of_load = fraction.get(Loading.FRACTION)
            fractionized_load = reduced_load * fraction_of_load
            if border == 0:
                x, y = [-dx, 0, dx], [0, fractionized_load / dx, 0]
            else:
                skewedness = fraction.get(Loading.SKEWEDNESS)

                # calculate peak loading [unit/m]
                x, y = self.spread(border, fractionized_load, skewedness, dx)

            series.append((x, y))

        # combine all fractions to one combined timeseries
        xcomb, ycomb = self.superimpose(series, dx)

        # convert loading from unit/m to unit/volume
        flow_rate = qlookup.lookup_q(link, ref_time)[0]
        ccomb = ycomb * flow_velocity / flow_rate

        # convert spatial array to temporal array
        t = ref_time + datetime.timedelta(
            seconds=Discretization.TIMESTEPLENGTH.total_seconds()) * -xcomb / flow_velocity
        t = qlookup.norm_array(t)

        # check mass-balance
        totalm = sum(ccomb) * Discretization.TIMESTEPLENGTH.total_seconds() * flow_rate
        deviation = abs(totalm - reduced_load) / reduced_load
        if not deviation <= Loading.MAX_CONTINUITY_ERROR:
            print(f"Expected load: {reduced_load}\n"
                  f"Calculated load: {totalm}\n"
                  f"Continuity Error: {deviation * 100:.2f} %\n"
                  f"Allowed Error: {Loading.MAX_CONTINUITY_ERROR * 100:.2} %")
            raise MassbalanceError
        return pd.Series(ccomb, index=t, name=self.pid)

    def unpack_loading(self, values, timestamps, qlookup):
        """
        unpacks values to full day series. Pads all new values with 0
        Args:
            values (list):
            timestamps (list):
            qlookup (QSeries):

        Returns:
            list: Extended list of load values
        """
        indexers, values = timestamps, values
        serieslength = len(qlookup.timestamps)
        valcount = len(values)

        unpacked = np.zeros(serieslength)
        first = qlookup.timestamps.index(qlookup.norm_time(indexers[0]))
        try:
            unpacked[first:first + valcount] = values
        except:
            firstsnipped = serieslength - first
            unpacked[first:] = values[:firstsnipped]
            unpacked[:valcount - firstsnipped] = values[firstsnipped:]
        return unpacked

    def get_entry(self, constituent, link, qlookup):
        """
        calculates load and unpacks it to return as an entry dictionary to be used in postprocessing object
        Args:
            constituent:
            link:
            qlookup:

        Returns:
            dict: entry dict to be used in Postprocessing object
        """
        # values = self.calc_load(constituent, link, qlookup)
        values = self.calc_composite_load(constituent, link, qlookup)
        # values = self.unpack_loading(values, timestamps, qlookup)
        tags = {}
        return {"values": values, 'origin': self.origin, 'traveltime': self.age, 'class': self.classification,
                "pid": self.pid}

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
            print(f"tm: {tm}, t0: {self.t0}")
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
    tm = t0 + dt.timedelta(seconds=30 * 60)
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
    loads = x.get_entry(Loading.FECAL, link, qlut)
    print(loads)
    print("main finished")


if __name__ == "__main__":
    main()
    print("Module finished")
