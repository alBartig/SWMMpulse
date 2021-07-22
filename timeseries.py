from Exceptions import *
import datetime as dt
import numpy as np
import pandas as pd
from tqdm import tqdm

# class TStamps:
#     def __init__(self, timestamps, freq):
#         timestamps = pd.to_datetime(timestamps)
#         self.date = timestamps[0].date()
#         self.timestamps = pd.date_range(dt.datetime.combine(self.date,dt.time(hour=0,minute=0,second=0)),
#                                         dt.datetime.combine(self.date,dt.time(hour=23,minute=59,second=59)),
#                                         freq=freq)
#
#
# class TData:
#     def __init__(self, timestamps, entries=None, freq="10S", expand=True, eid=None, **kwargs):
#         self._reg_datetimeindex(timestamps, freq)
#         self.tags = []
#
#     def set_tags(self, tags):
#         self.tags = tags
#         return None
#
#     def _reg_datetimeindex(self, timestamps, freq):
#         """
#         expands entered timestamps to specified frequency
#         Args:
#             timestamps (list): list of datetime values or timestamps
#             freq: Desired frequency, default = "10S"
#
#         Returns:
#             None
#         """
#         timestamps = pd.to_datetime(timestamps)
#         self.date = timestamps[0].date()
#         self.timestamps = pd.date_range(dt.datetime.combine(self.date, dt.time(hour=0, minute=0, second=0)),
#                                         dt.datetime.combine(self.date, dt.time(hour=23, minute=59, second=59)),
#                                         freq=freq)
#         return None
#
#     def _expand_series(self, entry):
#         """
#         Maps an input series on the stored datetimeindex. Missing values are interpolated.
#         Entry must have datetime index
#         Args:
#             entry (pd.Series):
#
#         Returns:
#             pd.Series, expanded to stored datetimeindex
#         """
#         expanded = pd.Series(entry, index=self.timestamps).interpolate(method="time", limit_direction="both")
#         return expanded
#
#     def _append_entry(self, entry, tags=None):
#         s = pd.DataFrame(self._expand_series(entry))
#         if tags:
#             try:
#                 tagnames = self.tags
#             except:
#                 print("Cannot tag series before assigning tags to object.\nAssign tags with 'set_tags()'")
#                 raise BaseException
#             if type(tags) is dict:
#                 mi = pd.MultiIndex.from_tuples([tuple(tags.get(name, "na") for name in tagnames)], names=tagnames)
#             elif type(tags) is list:
#                 mi = pd.MultiIndex.from_tuples([tuple(tags[:len(tagnames)])], names=tagnames)
#
#             s.columns = mi
#
#         try:
#             self.df = self.df.join(s)
#         except:
#             self.df = pd.DataFrame(s, index=self.timestamps)
#         return None
#
#     def _append_df(self, df):
#         for column in df:
#             self._append_entry(df[column])
#         return None
#
#     def _append_entries(entries, tags):
#         """
#         """
#         pass
#
# class TSeries:
#     def __init__(self, timestamps, entries, freq="10S", expand=True, eid=None, **kwargs):
#         self._reg_datetimeindex(timestamps, freq)
#
#     def _reg_datetimeindex(self, timestamps, freq):
#         """
#         expands entered timestamps to specified frequency
#         Args:
#             timestamps (list): list of datetime values or timestamps
#             freq: Desired frequency, default = "10S"
#
#         Returns:
#             None
#         """
#         timestamps = pd.to_datetime(timestamps)
#         self.date = timestamps[0].date()
#         self.timestamps = pd.date_range(dt.datetime.combine(self.date,dt.time(hour=0,minute=0,second=0)),
#                                         dt.datetime.combine(self.date,dt.time(hour=23,minute=59,second=59)),
#                                         freq=freq)
#         return None
#
#     def _expand_series(self, entry):
#         """
#         Maps an input series on the stored datetimeindex. Missing values are interpolated.
#         Entry must have datetime index
#         Args:
#             entry (pd.Series):
#
#         Returns:
#             pd.Series, expanded to stored datetimeindex
#         """
#         expanded = pd.Series(entry, index=self.timestamps).interpolate(method="time", limit_direction="both")
#         return expanded
#
class TDict:
    def __init__(self, timestamps):
        #time_values = pd.to_datetime(time_values)
        timestamps = pd.to_datetime(timestamps)
        self.date = timestamps[0].date()
        self.timestamps, self.sortby = self._normalize(timestamps)
        self.lookup = self.register_timeseries(self.timestamps)

    def _normalize(self, timestamps):
        temp = []
        for t in timestamps:
            if t.date() == self.date:
                temp.append(t)
            else:
                t = datetime.datetime.combine(self.date,t.time())
                if not temp.__contains__(t):
                    temp.append(t)
                else:
                    print("Series spans more than 24 h! Duplicate time-values found.\nPlease remove duplicate time-values regardless of date")
                    raise TimevalueError
        sortby = np.argsort(temp)
        return ([temp[i] for i in sortby],sortby)

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

    def _get_secs(self,time):
        return self.lookup.get(time.hour,{}).get(time.minute,[])

    def _get_secs_padded(self, time):
        """
        same as self._get_secs, but also queries seconds from 60s ahead to cover min + 60s
        Args:
            time:

        Returns:
            list: seconds in minute and if available following minute first value
        """
        secs1 = self._get_secs(time)
        time2 = time+datetime.timedelta(seconds=60)
        secs2 = [s+60 for s in self._get_secs(time2)[:1]]
        return  secs1 + secs2

    def _get_time_from_sec(self, time, secs):
        return datetime.datetime.combine(self.date, datetime.time(time.hour,time.minute)) + datetime.timedelta(seconds=secs)

    def _get_time_from_min(self, time, min, first=True):
        time = datetime.time(time.hour,min)
        secs = abs(first-1)*max(self._get_secs(time)) + first*min(self._get_secs(time))
        time.replace(second=secs)
        return datetime.datetime.combine(self.date,time)

    def find_closest(self, time):
        h = time.hour
        m = time.minute
        s = time.second
        minutes_in_hour = self.lookup.get(h)
        if minutes_in_hour != None:
            seconds_in_minute = self._get_secs_padded(time)
            if seconds_in_minute != None:
                s2 = min(seconds_in_minute, key = lambda x: abs(x-s))
                time = self._get_time_from_sec(time, s2)
            else:
                minutes = minutes_in_hour.keys()
                m2 = min(minutes, key = lambda x: abs(x-m))
                if m2 > m:
                    time = self._get_time_from_min(time, m2, first=False)
                    #time = datetime.time(hour=h,minute=m2,second=max(seconds_in_minute))
                else:
                    time = self._get_time_from_min(time, m2, first=True)
                    #time = datetime.time(hour=h,minute=m2,second=min(seconds_in_minute))
        return time

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

    def find_larger_dt(self, time):
        larger = [ts for ts in self.timestamps if ts >= time]
        if larger == []:
            return self.timestamps[0],datetime.timedelta(days=1)
        else:
            try:
                return larger[0],datetime.timedelta(days=0)
            except:
                print(f'larger: {larger}, time: {time}')


    def find_smaller_dt(self, time):
        querydate = time.date()
        querytime = time.time()
        compare_dt = datetime.datetime.combine(self.date, querytime)
        delay = datetime.timedelta(days=(querydate - self.date).days)
        smaller = [ts for ts in self.timestamps if ts <= compare_dt]
        return smaller[-1],delay

    def find_closest_dt(self, time):
        smaller,shift_front = self.find_smaller_dt(time)
        ismaller = self.get_index(smaller)
        if ismaller < len(self.timestamps)-1:
            ilarger = ismaller + 1
            larger = self.timestamps[ilarger]
            if larger-time < time-smaller:
                return larger
            else:
                return smaller
        else:
            return smaller

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
                    time = dt.time(hour=h,minute=m,second=s2)
                else:
                    time = time + dt.timedelta(seconds=60-s)
                    time = self.find_larger(time)
            else:
                time = time + dt.timedelta(seconds=60-s)
                time = self.find_larger(time)
        return datetime.datetime.combine(self.date, time)

    def get_index(self, querydate):
        querytime = querydate.time()
        compare_dt = dt.datetime.combine(self.date, querytime)
        try:
            return self.timestamps.index(compare_dt)
        except:
            raise BaseException

    def get_datetime(self, index):
        return self.timestamps[index]

    def _to_csv(self, content, opath):
        """

        Args:
            content (dict):
            path (str):

        Returns:

        """
        import csv
        import os
        import ntpath

        def path_leaf(path):
            head, tail = ntpath.split(path)
            return tail or ntpath.basename(head)

        def path_leaf_inverse(path):
            return path.strip('/\\').strip(path_leaf(path))

        if not os.path.exists(path_leaf_inverse(opath)):
            os.makedirs(path_leaf_inverse(opath))

        with open(opath, "w+", newline="") as csvfile:
            writer = csv.writer(csvfile, delimiter=';')
            writer.writerow(list(content.keys()))
            for row in zip(*list(content.values())):
                writer.writerow(row)
        print(f'csv-file saved to "{opath}"')
#
# class TSeries(TDict):
#     """
#     Timeseries object that offers additional functionality.
#     Consinsts of dictionary with following structure:
#     {timestamps: timestamps,
#     entries: [{values: values, **kwargs},
#               {values: values, **kwargs},...]
#     locations:location, constituent:constituent, other **kwargs}
#     """
#     def __init__(self, timestamps, entries, expand=True, eid=None, **kwargs):
#         super().__init__(timestamps)
#         if eid is not None:
#             try:
#                 self.eids = [entry[eid] for entry in entries]
#             except:
#                 print(f'Could not find {eid} in entry')
#         self.entries = []
#         if expand == True:
#             [self.append(entry,expand) for entry in entries]
#         else:
#             for entry in entries:
#                 try:
#                     vals = [entry['values'][i] for i in self.sortby]
#                 except:
#                     pass
#                 entry['values'] = vals
#                 self.append(entry, expand)
#         [self.__setattr__(key, value) for key, value in kwargs.items()]
#
#     def _expand(self, values, timestamps):
#         """
#         creates array with zeros of the size of self.count. Inserts the passed values at the correct indexes
#         Args:
#             values (list): List of values to add
#             timestamps (list): timestamps corresponding to the values
#
#         Returns:
#             np.array: Array of the size of the timeseries with passed values at correct fields
#         """
#         try:
#             expanded = np.zeros(self.count)
#             i0 = self.get_index(timestamps[0])
#             ie = self.get_index(timestamps[-1])
#             if ie-i0 == len(timestamps)-1:
#                 for i in range(len(timestamps)):
#                     expanded[i0+i] = values[i]
#             else:
#                 for time,value in zip(timestamps,values):
#                     index = self.get_index(time)
#                     expanded[index] = value
#         except:
#             raise Exception
#         return expanded
#
#     def _explode_eid(self, eid):
#         #extract eid as series
#         series = self._extract_eid(eid)
#         #create high-res daterange
#         dr = pd.date_range(self.date, periods=8640, freq="10S")
#         #create empty DataFrame with high-res dr as index
#         df = pd.DataFrame(index=dr)
#         #join with old series
#         df = df.join(series)
#         #interpolate missing values
#         df.interpolate(inplace=True)
#         return df[series.name]
#
#     def _extract_eid(self, eid):
#         return pd.Series(self.entry(eid)['values'], index=self.timestamps, name=eid)
#
#     def append(self, entrydict, expand=True):
#         """
#         appends a time-value - tuple from a Pobject-class object to the tdict.
#         Creates an array over entire ts-lengths with 0 values in unused fields
#         Args:
#             packetseries (dict): dictionary, must contain 'values' and 'timestamps'. Usually returned by Pobject.get_loads()
#
#         Returns:
#             bool: true if appending was successful, False if not
#         """
#         if expand is True:
#             entrydict['values'] = self._expand(entrydict.pop('values'),entrydict.pop('timestamps'))
#             self.entries.append(entrydict)
#             return True
#         else:
#             self.entries.append(entrydict)
#
#     def entry(self, eid):
#         i = self.eids.index(eid)
#         return self.entries[i]
#
#     def distribution(self, time, key='traveltime'):
#         """
#         plots distribution of occurences in a given timewindow
#         Args:
#             time (datetime.datetime or tuple: either a given point in time or a tuple of start and end time
#             key (string): Value, the distribution is looked for. Default is 'traveltime'
#
#         Returns:
#
#         """
#         #generate values
#         if type(time) == tuple or type(time) == list:
#             start = self.get_index(self.find_smaller(time[0]))
#             end = self.get_index(self.find_larger(time[1]))
#             values = [(entry[key],np.mean(entry['values'][start:end])) for entry in self.entries]
#         elif type(time) == datetime.datetime:
#             index = self.get_index(self.find_closest(time))
#             values = [(entry[key],entry['values'][index]) for entry in self.entries]
#
#     def timeseries(self, start=None, end=None, wpath=None):
#         """
#         returns all entries aggregated into one series along with the timesteps
#         Args:
#             start (datetime.datetime): start of timeseries (optional), default: first value in timestamps
#             end (datetime.datetime): end of timeseries (optional), default: last value in timestamps
#
#         Returns:
#             tuple: (values,timestamps)
#         """
#         #set parameters
#         if start == None:
#             start = self.get_index(self.start)
#         else:
#             start = self.get_index(self.find_smaller(start))
#
#         if end == None:
#             end = min(self.get_index(self.end)+1,len(self.timestamps))
#         else:
#             end = min(self.get_index(self.find_larger(end))+1,len(self.timestamps))
#
#         #aggregate entries
#         timestamps = self.timestamps[start:end]
#         timeseries = [0] * len(timestamps)
#         timeseries = sum([entry['values'][start:end] for entry in self.entries], timeseries)
#
#         if wpath is not None:
#             self._to_csv({'timestamps':self.timestamps,'values':timeseries},wpath)
#
#         return (timestamps,timeseries)
#
#     def traveltime_distribution(self, start=None, end=None):
#         #creating bins:
#         nbins = 10
#         traveltimes = set([entry["Traveltime"] for entry in self.entries])
#         tmin, tmax = min(traveltimes), max(traveltimes)
#         binwidth = round((tmax-tmin)/nbins)
#         binspecs = [(tmin + i * binwidth, tmin + (i+1) * binwidth) for i in range(nbins)]
#
#
#     def plot(self, eid=None, **kwargs):
#         """
#         creates plot of timeseries
#         Args:
#             **kwargs: start time, end time
#         """
#         import matplotlib.pyplot as plt
#         fig = plt.figure()
#         ax1 = fig.add_subplot(111)
#         x,y = self.timeseries()
#         if eid is None:
#             ax1.plot(x,y)
#         else:
#             ax1.plot(self.timestamps[0], self.entry(eid)['values'])
#         plt.savefig(f'/mnt/c/Users/albert/Documents/SWMMpulse/tseries_plot_{eid}.png')

class TSeries:
    freqmult = {"S": 1,
                "M": 60,
                "H": 3600}

    def __init__(self, timestamps, entries=None, freq="10S", tagnames=[]):
        self._reg_datetimeindex(timestamps, freq)
        self.tags = tagnames
        self.freq = int(freq[:-1] * self.freqmult[freq[-1:]])

        if entries:
            self._append_entries(entries, timestamps)

    def set_tags(self, tags):
        self.tags = tags
        return None

    def norm_time(self, time):
        return dt.datetime.combine(self.date, time.time())

    def _reg_datetimeindex(self, timestamps, freq):
        """
        expands entered timestamps to specified frequency
        Args:
            timestamps (list): list of datetime values or timestamps
            freq: Desired frequency, default = "10S"

        Returns:
            None
        """
        timestamps = pd.to_datetime(timestamps)
        self.date = timestamps[0].date()
        self.timestamps = pd.date_range(dt.datetime.combine(self.date, dt.time(hour=0, minute=0, second=0)),
                                        dt.datetime.combine(self.date, dt.time(hour=23, minute=59, second=59)),
                                        freq=freq)
        return None

    def find_closest_dt(self, time):
        secs = time.second

        m = secs % self.freq
        if m < self.freq / 2:
            delt = dt.timedelta(seconds=(-m))
        else:
            delt = dt.timedelta(seconds=(self.freq - m))
        return time + delt

    def _expand_series(self, entry):
        """
        Maps an input series on the stored datetimeindex. Missing values are interpolated.
        Entry must have datetime index
        Args:
            entry (pd.Series):

        Returns:
            pd.Series, expanded to stored datetimeindex
        """
        expanded = pd.Series(entry, index=self.timestamps).interpolate(method="time", limit_direction="both")
        return expanded

    def _append_series(self, series, tags=None):
        s = pd.DataFrame(self._expand_series(series))
        if tags:
            try:
                tagnames = self.tags
            except:
                print("Cannot tag series before assigning tags to object.\nAssign tags with 'set_tags()'")
                raise BaseException
            if type(tags) is dict:
                mi = pd.MultiIndex.from_tuples([tuple(tags.get(name, "na") for name in tagnames)], names=tagnames)
            elif type(tags) is list:
                mi = pd.MultiIndex.from_tuples([tuple(tags[:len(tagnames)])], names=tagnames)

            s.columns = mi

        try:
            self.df = self.df.join(s)
        except:
            self.df = pd.DataFrame(s, index=self.timestamps)
        return None

    def _append_df(self, df):
        for column in df:
            self._append_entry(df[column])
        return None

    def _append_entry(self, entry, timestamps=None):
        if timestamps is None:
            s = pd.Series(entry["values"], index=entry["timestamps"])
        else:
            s = pd.Series(entry["values"], index=timestamps)

        self._append_series(s, tags=entry["tags"])
        return None

    def _append_entries(self, entries, timestamps):
        """
        """
        for entry in entries:
            self._append_entry(entry, timestamps=timestamps)

# class QSeries(TSeries):
#     """
#     Creates a lookup table for flow velocity from a swmm out-file
#     """
#     def __init__(self, out_file, **kwargs):
#         from swmm_api import read_out_file
#         with read_out_file(out_file) as out:
#             df = out.to_frame()
#             df = df.xs(('link','Flow_velocity'),axis=1,level=[0,2])
#             timestamps = df.index.values
#             entries = [{'values':df[col].values,'link':col} for col in df]
#             super().__init__(timestamps, entries, expand=False, eid='link')
#
#     def lookup_v(self, link, time):
#         """
#         Returns velocity in link at time. If v == 0 in link at time, the first v > 0
#         is returned together with delay time. If v never rises above 0 after, it returns None,None
#         Args:
#             link (str):
#             time (str):
#
#         Returns:
#             tuple: (velocity, delay)
#         """
#         if time.date() != self.date:
#             time = datetime.datetime.combine(self.date, time.time())
#         querytime = self.find_closest_dt(time)
#         index = self.get_index(querytime)
#         iorg = index
#         try:
#             vvalues = self.entries[self.eids.index(link)]['values']
#         except:
#             pass
#         v = vvalues[index]
#         imax = self.count - 1
#         while v == 0:
#             if index != imax: #Reset index to 0 if at end of array
#                 index += 1
#             else:
#                 index = 0
#             v = vvalues[index]
#             if index == iorg: #Break if full loop was done
#                 return None,None
#         if index > iorg: #Return if index has not looped around index > time
#             delay = self.get_datetime(index)-time
#         elif index < iorg: #Return if index has looped - index < time
#             delay = datetime.timedelta(days=1) - (time-self.get_datetime(index))
#         else:
#             delay = datetime.timedelta(seconds=0)
#         if delay < datetime.timedelta(seconds=0):
#             print(f'delay: {delay}, time: {time}, ')
#             raise PlausibilityError()
#         return v,delay
#
#     def m_to_s(self, link, time, distance):
#         '''
#         :param link: str; Link to be queried for
#         :param time: datetime obj, Time to be queried for
#         :param distance: float, Distance to be converted into time
#         :return: float, flowtime difference in seconds
#         '''
#         velocity,delay = self.lookup_v(link,time)
#         s = distance/velocity
#         return s

class QSeries():
    freqmult = {"S": 1,
                "M": 60,
                "H": 3600}

    def __init__(self, out_file, freq="10S"):
        dftemp = self.load_ofile(out_file)
        self.date, timestamps = self._reg_datetimeindex(dftemp.index, freq)
        self.freq = int(freq[:-1] * self.freqmult[freq[-1:]])
        dfq = pd.DataFrame(index=timestamps)
        self.dfq = dfq.join(dftemp).interpolate(limit_direction="both")

    def norm_time(self, time):
        return dt.datetime.combine(self.date, time.time())

    @staticmethod
    def find_seconds(time, freq):
        secs = time.second

        m = secs % freq
        if m < freq / 2:
            delt = dt.timedelta(seconds=(-m))
        else:
            delt = dt.timedelta(seconds=(freq - m))
        return time + delt

    @staticmethod
    def _reg_datetimeindex(timestamps, freq):
        """
        expands entered timestamps to specified frequency
        Args:
            timestamps (list): list of datetime values or timestamps
            freq: Desired frequency, default = "10S"

        Returns:
            None
        """
        timestamps = pd.to_datetime(timestamps)
        date = timestamps[0]
        timestamps = pd.date_range(date.floor("D"), date.ceil("D"), freq=freq)[:-1]
        return date, timestamps

    def load_ofile(self, out_file, **kwargs):
        from swmm_api import read_out_file
        with read_out_file(out_file) as out:
            df = out.to_frame()
            df = df.xs(('link', 'Flow_velocity'), axis=1, level=[0, 2])
            timestamps = df.index.values
            entries = [{'values': df[col].values, 'link': col} for col in df]
            return df

    def lookup_v(self, link, time):
        time = self.norm_time(time)
        s = self.dfq[link]
        try:
            v = s[time]
        except:
            v = s[self.norm_time(self.find_seconds(time, self.freq))]
        if v != 0.0:
            return v, dt.timedelta(seconds=0)
        else:
            s = s[s > 0]
            try:
                temp = s[time:]
                v, t = temp[0], temp.index[0]
                return v, t-time
            except:
                try:
                    temp = s[:time]
                    v, t = temp[0], temp.index[0]
                    return v, dt.timedelta(days=1) - (time-t)
                except:
                    return None,None

    def meter_to_seconds(self, link, time, distance):
        velocity, delay = self.lookup_v(link, time)
        s = distance / velocity
        return s

def test_qseries():
    print('Test QSeries')
    #of_path = 'C:/Users/alber/Documents/swmm/swmmpulse/HS_calib_120_simp.out'
    of_path = '/mnt/c/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.out'
    qlut = QSeries(of_path)
    print('Test QSeries finished')

def test_diffs():
    start, end = qlut.timestamps[0], qlut.timestamps[0] + datetime.timedelta(days=1)
    timestamps = [start + i * Discretization.TIMESTEPLENGTH for i in
                  range(int((end - start) / Discretization.TIMESTEPLENGTH))]
    tlut = TDict(timestamps)

    t0 = dt.datetime.now()
    deltas = [i for i in range(86400)]
    diffs = []
    times = []
    for delta in deltas:
        time = t0 + dt.timedelta(seconds=delta)
        closest = tlut.find_closest(time)
        diffs.append(min((closest - time).seconds, abs((closest - time).seconds - 86400)))
        times.append(time)
    df = pd.DataFrame.from_records(zip(deltas, diffs, times), columns=['Deltas', 'Differenzen', 'Timestamps'])
    print(f"Mean: {np.mean(diffs)} s\n"
          f"Min: {np.min(diffs)} s\n"
          f"Max: {np.max(diffs)} s")
    print("Finished Main")

def load_qlut():
    lpath = '/mnt/c/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.out'
    qlut = QSeries(lpath)
    return qlut

def load_graph():
    from network_lib import DirectedTree as ntwk, GraphConstants as gc
    gpath = '/mnt/c/Users/albert/documents/SWMMpulse/HS_calib_120_simp/'
    graph = ntwk.from_directory(gpath)
    return graph

def test_pproc(graph, qlut, evalnode):
    from prouter import _Route_table, Postprocessing
    import os
    from environment import Environment

    fpath = f'/mnt/c/Users/albert/documents/SWMMpulse/'
    #fpath = f'C:/Users/alber/documents/swmm/swmmpulse/'
    routetable = _Route_table.from_file(os.path.join(fpath,'route_tables','rtb_testdrive.pickle'))
    location = [evalnode,graph.get_outletlinks(evalnode)[0]]
    rt_slice = routetable._extract_node(location, Environment().dispersion)

    pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)
    eval_constituent = Loading.FECAL
    entry_loc = os.path.join(fpath, 'entries', "entries_testdrive.pickle")
    pproc.process_constituent(eval_constituent, entry_loc, load=False)
    pproc.Fecal_Matter.timeseries(wpath=os.path.join(fpath, 'timeseries', 'timeseries_testdrive.pickle'))
    print("Function finished")

def main():
    from application import testdrive

    evalnode = 'MH327-088-1'
    graph = load_graph()
    qlut = load_qlut()

    test_pproc(graph, qlut, evalnode)


if __name__ == '__main__':
    main()