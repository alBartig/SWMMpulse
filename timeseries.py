from Exceptions import *
import datetime
import numpy as np
import pandas as pd

class TDict:
    def __init__(self, timestamps):
        #time_values = pd.to_datetime(time_values)
        self.timestamps = pd.to_datetime(timestamps)
        self.date = self.timestamps[0].date()
        self.lookup = self.register_timeseries(self.timestamps)

    @classmethod
    def from_bounds(cls, start, end, step):
        """
        Initializes Timeseries
        Args:
            start (datetime.datetime):
            end (datetime.datetime):
            step (datetime.timedelta):
            **kwargs:
        """
        count,rest = divmod((end - start),step)
        if rest.total_seconds() != 0:
            raise StepsizeError
        timestamps = np.array([start + i * step for i in range(count)])
        return cls(timestamps)

    @property
    def start(self):
        return self.timestamps[0]

    @property
    def end(self):
        return self.timestamps[-1]

    @property
    def count(self):
        return len(self.timestamps)

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
        minutes_in_hour = self.lookup.get(h)
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
        minutes_in_hour = self.lookup.get(h)
        if minutes_in_hour != None:
            seconds_in_minute = minutes_in_hour.get(m)
            if seconds_in_minute != None:
                seconds = [second for second in seconds_in_minute if second < s]
                if len(seconds) > 1:
                    s2 = seconds[-1]
                    time =  datetime.time(hour=h,minute=m,second=s2)
                else:
                    time = time - datetime.timedelta(seconds=s+1)
                    time = self.find_smaller(time)
            else:
                time = time - datetime.timedelta(seconds=s+1)
                time = self.find_smaller(time)
        return datetime.datetime.combine(self.date, time)

    def find_larger(self, time):
        h = time.hour
        m = time.minute
        s = time.second
        minutes_in_hour = self.lookup.get(h)
        if minutes_in_hour != None:
            seconds_in_minute = minutes_in_hour.get(m)
            if seconds_in_minute != None:
                seconds = [second for second in seconds_in_minute if second<s]
                if len(seconds) > 1:
                    s2 = seconds[0]
                    time = datetime.time(hour=h,minute=m,second=s2)
                else:
                    time = time + datetime.timedelta(seconds=60-s)
                    time = self.find_larger(time)
            else:
                time = time + datetime.timedelta(seconds=60-s)
                time = self.find_larger(time)
        return datetime.datetime.combine(self.date, time)

    def get_index(self, datetime):
        return self.timestamps.get_loc(datetime)

    def get_datetime(self, index):
        return self.timestamps[index]

class TSeries(TDict):
    """
    Timeseries object that offers additional functionality.
    Consinsts of dictionary with following structure:
    {timestamps: timestamps,
    entries: [{values: values, **kwargs},
              {values: values, **kwargs},...]
    locations:location, constituent:constituent, other **kwargs}
    """
    def __init__(self, timestamps, entries, expand=True, eid=None, **kwargs):
        super().__init__(timestamps)
        if eid is not None:
            try:
                self.eids = [entry[eid] for entry in entries]
            except:
                print(f'Could not find {eid} in entry')
        self.entries = []
        [self.append(entry, expand) for entry in entries]
        [self.__setattr__(key, value) for key, value in kwargs.items()]

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
            index = np.where(self.timestamps == time)
            expanded[index] = value
        return expanded

    def append(self, entrydict, expand=True):
        """
        appends a time-value - tuple from a Pobject-class object to the tdict.
        Creates an array over entire ts-lengths with 0 values in unused fields
        Args:
            packetseries (dict): dictionary, must contain 'values' and 'timestamps'. Usually returned by Pobject.get_loads()

        Returns:
            bool: true if appending was successful, False if not
        """
        if expand is True:
            try:
                entrydict['values'] = self._expand(entrydict.pop('values'),entrydict.pop('timestamps'))
                self.entries.append(entrydict)
                return True
            except:
                print('Timestamps list expected, not found. Set kwarg "expand" to False to append without padding values')
                raise PacketdictError
        else:
            self.entries.append(entrydict)

    def entry(self, eid):
        i = self.eids.index(eid)
        return self.entries[i]

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
            start = self.get_index(self.find_smaller(time[0]))
            end = self.get_index(self.find_larger(time[1]))
            values = [(entry[key],np.mean(entry['values'][start:end])) for entry in self.entries]
        elif type(time) == datetime.datetime:
            index = self.get_index(self.find_closest(time))
            values = [(entry[key],entry['values'][index]) for entry in self.entries]

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
            start = self.get_index(self.start)
        else:
            start = self.get_index(self.tlookup.find_smaller(start))

        if end == None:
            end = self.get_index(self.end)
        else:
            end = self.get_index(self.find_larger(end))

        #aggregate entries
        timeseries = sum([entry['values'][start:end] for entry in self.entries])
        return (self.timestamps,timeseries)

    def plot(self, eid=None, **kwargs):
        """
        creates plot of timeseries
        Args:
            **kwargs: start time, end time
        """
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        if eid is None:
            ax1.plot(self.timestamps, self.timeseries())
        else:
            ax1.plot(self.timestamps, self.entry(eid)['values'])
        plt.savefig(f'/mnt/c/Users/albert/Documents/SWMMpulse/tseries_plot_{eid}.png')


class QSeries(TSeries):
    """
    Creates a lookup table for flow velocity from a swmm out-file
    """
    def __init__(self, out_file, **kwargs):
        from swmm_api import read_out_file
        with read_out_file(out_file) as out:
            df = out.to_frame()
            df = df.xs(('link','Flow_velocity'),axis=1,level=[0,2])
            timestamps = df.index.values
            entries = [{'values':df[col].values,'link':col} for col in df]
            super().__init__(timestamps, entries, expand=False, eid='link')

    def lookup_v(self, link, time):
        """
        Returns velocity in link at time. If v == 0 in link at time, the first v > 0
        is returned together with delay time. If v never rises above 0 after, it returns None,None
        Args:
            link (str):
            time (str):

        Returns:
            tuple: (velocity, delay)
        """
        querytime = self.find_closest(time)
        index = self.get_index(querytime)
        iorg = index
        vvalues = self.entries[self.eids.index(link)]['values']
        v = vvalues[index]
        imax = self.count - 1
        while v == 0:
            if index != imax: #Reset index to 0 if at end of array
                index += 1
            else:
                index = 0
            v = vvalues[index]
            if index == iorg: #Break if full loop was done
                return None,None
        else:
            if index >= iorg: #Return if index has not looped around index > time
                delay = self.get_datetime(index)-time
                return v,delay
            else: #Return if index has looped - index < time
                delay = datetime.timedelta(days=1) - (time-self.get_datetime(index))
                return v,delay

    def m_to_s(self, link, time, distance):
        '''
        :param link: str; Link to be queried for
        :param time: datetime obj, Time to be queried for
        :param distance: float, Distance to be converted into time
        :return: float, flowtime difference in seconds
        '''
        velocity,delay = self.lookup_v(link,time)
        s = distance/velocity
        return s

if __name__ == '__main__':
    print('Test QSeries')
    #of_path = 'C:/Users/alber/Documents/swmm/swmmpulse/HS_calib_120_simp.out'
    of_path = '/mnt/c/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.out'
    qlut = QSeries(of_path)
    print('Test QSeries finished')