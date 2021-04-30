from Exceptions import *
from swmm_api import read_out_file
import datetime
import numpy as np
import pandas as pd

class Timeseries:
    import pandas as pd
    """
    Timeseries object that offers additional functionality.
    Consinsts of dictionary with following structure:
    {timestamps: timestamps,
    entries: [entry1: {values: values, **kwargs},
              entry2: {values: values, **kwargs},...]
    locations:location, constistuent:constituent, other **kwargs}
    """
    def __init__(self, start, end, step,**kwargs):
        self.count,rest = divmod((start-end).total_seconds(),step)
        if rest != 0:
            raise StepsizeError
        timestamps = np.array([start + datetime.timedelta(seconds=i*step) for i in range(self.count)])
        self.tdict = {'timestamps':timestamps,'entries':[],**kwargs}
        self.tlookup = Tlookup(timestamps)

    @property
    def start(self):
        return self.tdict['timestamps'][0]
    @property
    def end(self):
        return self.tdict['timestamps'][-1]

    def append(self, packetdict):
        """
        appends a time-value - tuple from a Pobject-class object to the tdict.
        Creates an array over entire ts-lengths with 0 values in unused fields
        Args:
            packetseries (dict): dictionary, must contain 'values' and 'timestamps'. Usually returned by Pobject.get_loads()

        Returns:
            bool: true if appending was successful, False if not
        """
        try:
            packetdict['values'] = self.expand(packetdict.pop('values'),packetdict.pop('timestamps'))
            self.tdict['entries'].append([packetdict])
            return True
        except:
            raise PacketdictError

    def _expand(self, values, timestamps):
        """
        creates array with zeros of the size of self.count. Inserts the passed values at the correct indexes
        Args:
            values (list): List of values to add
            timestamps (list): timestamps corresponding to the values

        Returns:
            np.array: Array of the size of the timeseries with passed values at correct fields
        """
        expanded = np.zeros(self.count)
        for time,value in zip(timestamps,values):
            index = np.where(self.timeseries == time)
            expanded[index] = value
        return expanded

    @property
    def timeseries(self, start=None, end=None):
        """
        returns all entries aggregated into one series along with the timesteps
        Args:
            start (datetime.datetime): start of timeseries (optional), default: first value in timestamps
            end (datetime.datetime): end of timeseries (optional), default: last value in timestamps

        Returns:
            tuple: (values,timestamps)
        """
        #set parameters
        if start == None:
            start = self.tlookup.get_index(self.start)
        else:
            start = self.tlookup.get_index(self.tlookup.find_smaller(start))

        if end == None:
            end = self.tlookup.get_index(self.end)
        else:
            end = self.tlookup.get_index(self.tlookup.find_larger(end))

        #aggregate entries
        timeseries = sum([entry['values'][start:end] for entry in self.tdict['entries']])
        return (self.timestamps,timeseries)

    def plot(self, **kwargs):
        """
        creates plot of timeseries
        Args:
            **kwargs: start time, end time
        """
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = fig.add_subplots(111)
        ax1.plot(*self.timeseries(**kwargs))

    def distribution(self, time, key='traveltime'):
        """
        plots distribution of occurences in a given timewindow
        Args:
            time (datetime.datetime or tuple: either a given point in time or a tuple of start and end time
            key (string): Value, the distribution is looked for. Default is 'traveltime'

        Returns:

        """
        #generate values
        if type(time) == tuple or type(time) == list:
            start = self.tlookup.get_index(self.tlookup.find_smaller(time[0]))
            end = self.tlookup.get_index(self.tlookup.find_larger(time[1]))
            values = [(entry[key],np.mean(entry['values'][start:end])) for entry in self.tdict['entries']]
        elif type(time) == datetime.datetime:
            index = self.tlookup.get_index(self.tlookup.find_closest(time))
            values = [(entry[key],entry['values'][index]) for entry in self.tdict['entries']]



class Tlookup:
    def __init__(self, time_values):
        time_values = pd.to_datetime(time_values)
        self.times = time_values
        self.date = self.times[0].date()
        self.timedict = self.register_timeseries(self.times)

    def register_timeseries(self,series):
        timedict = {}
        for step in series:
            h = step.hour
            m = step.minute
            s = step.second
            if timedict.__contains__(h):
                if timedict[h].__contains__(m):
                    timedict[h][m].append(s)
                else:
                    timedict[h][m] = [s]
            else:
                timedict[h] = {m:[s]}

        return timedict

    def find_closest(self, time):
        h = time.hour
        m = time.minute
        s = time.second
        minutes_in_hour = self.timedict.get(h)
        if minutes_in_hour != None:
            seconds_in_minute = minutes_in_hour.get(m)
            if seconds_in_minute != None:
                s2 = min(seconds_in_minute, key = lambda x: abs(x-m))
                time =  datetime.time(hour=h,minute=m,second=s2)
            else:
                minutes = minutes_in_hour.keys()
                m2 = min(minutes, key = lambda x: abs(x-m))
                seconds_in_minute = minutes_in_hour.get(m2)
                if m2 > m:
                    time = datetime.time(hour=h,minute=m2,second=max(seconds_in_minute))
                else:
                    time = datetime.time(hour=h,minute=m2,second=min(seconds_in_minute))
        return datetime.datetime.combine(self.date, time)

    def find_smaller(self, time):
        h = time.hour
        m = time.minute
        s = time.second
        minutes_in_hour = self.timedict.get(h)
        if minutes_in_hour != None:
            seconds_in_minute = minutes_in_hour.get(m)
            if seconds_in_minute != None:
                s2 = min([second for second in seconds_in_minute if second<s], key = lambda x: abs(x-m))
                time =  datetime.time(hour=h,minute=m,second=s2)
            else:
                minutes = minutes_in_hour.keys()
                m2 = min([minute for minute in minutes if minute<m], key = lambda x: abs(x-m))
                seconds_in_minute = minutes_in_hour.get(m2)
                time = datetime.time(hour=h,minute=m2,second=max[seconds_in_minute])
        return datetime.datetime.combine(self.date, time)

    def find_larger(self, time):
        h = time.hour
        m = time.minute
        s = time.second
        minutes_in_hour = self.timedict.get(h)
        if minutes_in_hour != None:
            seconds_in_minute = minutes_in_hour.get(m)
            if seconds_in_minute != None:
                s2 = min([second for second in seconds_in_minute if second>s], key = lambda x: abs(x-m))
                time =  datetime.time(hour=h,minute=m,second=s2)
            else:
                minutes = minutes_in_hour.keys()
                m2 = min([minute for minute in minutes if minute>m], key = lambda x: abs(x-m))
                seconds_in_minute = minutes_in_hour.get(m2)
                time =  datetime.time(hour=h,minute=m2,second=min[seconds_in_minute])
        return datetime.datetime.combine(self.date, time)

    def get_index(self, datetime):
        return self.times.index(datetime)

    def get_datetime(self, index):
        return self.times[index]


class qlookup:
    def __init__(self, fpath):
        out = read_out_file(fpath)
        lookup = out.to_frame()
        lookup = lookup.xs(('link','Flow_velocity'),axis=1,level=[0,2])
        self.table = lookup
        self.date = self.table.index[0].date()
        self.tlookup = Tlookup(self.table.index.values)

    def lookup_v(self, link, time):
        querytime = self.tlookup.find_closest(time)
        v = self.table.at[querytime,link]
        return v

    def m_to_s(self, link, time, distance):
        '''
        :param link: str; Link to be queried for
        :param time: datetime obj, Time to be queried for
        :param distance: float, Distance to be converted into time
        :return: float, flowtime difference in seconds
        '''
        velocity = self.lookup_v(link,time)
        s = distance/velocity
        return s

if __name__ == '__main__':
    lpath = 'C:/Users/alber/Documents/swmm/swmmpulse/HS_calib_120_simp.out'
    lt = qlookup(lpath)
    print()