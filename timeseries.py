from Exceptions import *
import datetime
import numpy as np
import pandas as pd
from tqdm import tqdm

class TDict:
    def __init__(self, timestamps):
        self.timeframe = pd.DataFrame.from_records([[ts.hour,ts.minute,ts.second,ts] for ts in pd.to_datetime(timestamps)],
                                                   columns=['hours','minutes','seconds','timestamps'])
        self.date = self.timestamps[0].date()
        self.timedict = self.register_timeseries(self.timeframe.timestamps)

    @property
    def timestamps(self):
        return self.timeframe.timestamps

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

    def _index(self, time):
        return self.timedict.get(time.hour).get(time.minute).get(time.second)

    def _time(self, index):
        if index > (self.count - 1):
            index -= self.count
        elif index < 0:
            index += self.count
        return self.timestamps[index]

    def _closer(self, ref, t1, t2):
        times = [t1,t2]
        diffs = [abs(ref - t) for t in times]
        index = diffs.index(min(diffs))
        return times[index]

    def register_timeseries(self,series):
        timedict = {}
        for i,step in enumerate(series):
            h,m,s = step.hour, step.minute, step.second
            if timedict.__contains__(h):
                if timedict[h].__contains__(m):
                    timedict[h][m][s] = i
                else:
                    timedict[h][m] = {s:i}
            else:
                timedict[h] = {m: {s:i}}
        return timedict

    def _strp_time(self, time):
        return (time.hour, time.minute, time.second)

    def check_time(self, time):
        """
        Checks if time is contained in timedict and returns the missing key if not
        Args:
            time (datetime.datetime): querytime

        Returns:
            time (datetime.datetime): if exact time is contained in timedict
            key (string): string of ['hour','minute','second'] where the data failed to be retrieved
        """
        try:
            hours = self.timedict[time.hour]
        except:
            return 'hour'
        try:
            minutes = hours[time.minute]
        except:
            return 'minute'
        try:
            seconds = minutes[time.second]
        except:
            return 'second'
        return True

    def _close_sample(self, time, key):
        h,m,s = self._strp_time(time)
        if key == 'second':
            return self.timeframe.loc[(self.timeframe.hours == h) & (self.timeframe.minutes == m)].timestamps
        elif key == 'minute':
            return self.timeframe.loc[(self.timeframe.hours == h)].timestamps
        elif key == 'hour':
            return self.timeframe.hours.timestamps

    def get_closest(self, time):
        time = datetime.datetime.combine(self.date,time.time())
        key = self.check_time(time)
        if key is True:
            return time
        else:
            ls = list(self._close_sample(time,key))
            greater = self._greater_ls(time, ls)
            if greater:
                i = self._index(greater)
            else:
                i = self._index(ls[-1])+1
                greater = self._time(i)
            lesser = self._time(i-1)
            return self._closer(time,greater,lesser)

    def get_greater(self, time):
        querydate = time.date()
        delay = datetime.timedelta(days=(querydate - self.date).days)
        time = datetime.datetime.combine(self.date,time.time())
        key = self.check_time(time)
        if key is True:
            return time
        else:
            ls = list(self._close_sample(time,key))
            greater = self._greater_ls(time,ls)
            if greater:
                return greater,delay
            else:
                i = self._index(ls[-1])+1
                return self._time(i),delay

    def get_lesser(self, time):
        querydate = time.date()
        delay = datetime.timedelta(days=(querydate - self.date).days)
        time = datetime.datetime.combine(self.date,time.time())
        key = self.check_time(time)
        if key is True:
            return time
        else:
            ls = list(self._close_sample(time,key))
            lesser = self._lesser_ls(time,ls)
            if lesser:
                return lesser,delay
            else:
                i = self._index(ls[0])-1
                return self._time(i),delay

    def _lesser_ls(self, val, ls):
        i = len(ls) - 1
        try:
            while ls[i] > val:
                i -= 1
            return ls[i]
        except:
            return False

    def _greater_ls(self, val, ls):
        i = 0
        try:
            while ls[i] < val:
                i += 1
            return ls[i]
        except:
            return False

    def _closest_ls(self, val, ls):
        i = 0
        imax = len(ls) - 1
        try:
            while ls[i] < val and i < imax:
                i += 1
            if abs(ls[i] - val) < abs(ls[i - 1] - val):
                return ls[i]
            else:
                return ls[i - 1]
        except:
            return False

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
        [self.append(entry, expand) for entry in tqdm(entries)]
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
        try:
            expanded = np.zeros(self.count)
            for time,value in zip(timestamps,values):
                index = self._index(time)
                expanded[index] = value
        except:
            raise Exception
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
            entrydict['values'] = self._expand(entrydict.pop('values'),entrydict.pop('timestamps'))
            self.entries.append(entrydict)
            return True
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
            start = self._index(self.find_smaller(time[0]))
            end = self._index(self.find_larger(time[1]))
            values = [(entry[key],np.mean(entry['values'][start:end])) for entry in self.entries]
        elif type(time) == datetime.datetime:
            index = self._index(self.get_closest(time))
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
            start = self._index(self.start)
        else:
            start = self._index(self.tlookup.find_smaller(start))

        if end == None:
            end = min(self._index(self.end) + 1, len(self.timestamps))
        else:
            end = min(self._index(self.find_larger(end)) + 1, len(self.timestamps))

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
        if time.date() != self.date:
            time = datetime.datetime.combine(self.date, time.time())
        querytime = self.get_closest(time)
        index = self._index(querytime)
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
        if index > iorg: #Return if index has not looped around index > time
            delay = self._time(index) - time
        elif index < iorg: #Return if index has looped - index < time
            delay = datetime.timedelta(days=1) - (time - self._time(index))
        else:
            delay = datetime.timedelta(seconds=0)
        if delay < datetime.timedelta(seconds=0):
            print(f'delay: {delay}, time: {time}, ')
            raise PlausibilityError()
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