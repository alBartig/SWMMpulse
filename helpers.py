import random
import time
import datetime

def get_weighted_population(pop, fraction):
    n = 0
    for i in range(pop):
        p = random.random()
        if p <= fraction:
            n += 1
    return n

def test_finders():
    ls = [i for i in range(100)]
    testvals = [0,10.5,50,90,110]
    print("Smaller")
    start = time.time()
    for val in testvals:
        print(get_smaller_ls(val,ls))
    print(time.time()-start," s")
    print("Greater")
    start = time.time()
    for val in testvals:
        print(get_greater_ls(val,ls))
    print(time.time()-start," s")
    print("Closest")
    start = time.time()
    for val in testvals:
        print(get_closest_ls(val,ls))
    print(time.time()-start," s")

class Timetree:

    TIMEOUT = 10

    def __init__(self, timestamps):
        self.timeframe = pd.DataFrame.from_records([[ts.hour,ts.minute,ts.second,ts] for ts in timestamps],
                                                   columns=['hours','minutes','seconds','timestamps'])
        self.date = timestamps[0].date()
        self.timedict = self.register_timeseries(self.timeframe.timestamps)

    @property
    def timestamps(self):
        return self.timeframe.timestamps

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

    def check_hour(self, time):
        if self.timeframe.hours.unique().values.__contains__(time.hour):
            return True
        else:
            return False

    def check_minute(self, time):
        if self.timeframe.loc[self.timedict.hour == time.hour].minute.unique().values.__contains__(time.minute):
            return True
        else:
            return False

    def check_second(self, time):
        if self.timeframe.loc[(self.timedict.hour == time.hour) & (self.timedict.minute == time.minute)].second.unique().values.__contains__(time.second):
            return True
        else:
            return False

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
        date = time.date()
        time = datetime.datetime.combine(self.date,time.time())
        key = self.check_time(time)
        if key is True:
            return time
        else:
            ls = list(self._close_sample(time,key))
            greater = self.get_greater_ls(time,ls)
            if greater:
                i = self._index(greater)
            else:
                i = self._index(ls[-1])+1
                greater = self.timeframe.timestamps.loc[i]
            lesser = self.timeframe.loc[i-1].timestamps
            return self._closer(time,greater,lesser)

    def get_lesser_ls(self, val, ls):
        i = len(ls) - 1
        try:
            while ls[i] > val:
                i -= 1
            return ls[i]
        except:
            return False

    def get_greater_ls(self, val, ls):
        i = 0
        try:
            while ls[i] < val:
                i += 1
            return ls[i]
        except:
            return False

    def get_closest_ls(self, val, ls):
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

    def get_surrounding_ls(self, val, ls):
        i = 0
        imax = len(ls) - 1
        while ls[i] < val and i < imax:
            i += 1
        if ls[i] >= val:
            s = ls[i-1] if i > 0 else False
            g = ls[i]
            return s,g
        else:
            return ls[i],False

    def _closer(self, ref, t1, t2):
        times = [t1,t2]
        diffs = [abs(ref - t) for t in times]
        index = diffs.index(min(diffs))
        return times[index]

    def _closest_time(self, time):
        """
        Checks if time is in dict and if not, returns closest time available
        Args:
            time (datetime.datetime):

        Returns:
            time (datetime.datetime)
        """
        key = self.check_time(time)
        if key:
            return time
        elif key == 'second':
            return self._closest_second(time)
        elif key == 'minute':
            return self._closest_minute(time)
        elif key == 'hour':
            return self._closest_hour(time)

    def _closest_second(self, time):
        seconds = self.timedict.loc[(self.timedict.hour == time.hour) & (self.timedict.minute == time.minute)].second.unique()
        s,g = self.get_surrounding_ls(time.second,seconds)
        if s == False:
            s = self._get_smaller_second(time)
        if g == False:
            g = self._get_larger_second(time)
        return self._closer(time,s,g)

    def _closest_minute(self, time):
        minutes = self.timedict.loc[self.timedict.hour == time.hour].minute.unique()
        s,g = self.get_surrounding_ls(time.minute,minutes)
        ds,dg = [abs(x-time.minute) for x in [s,g]]
        return [s,g][[ds,dg].index(min(ds,dg))]

    def _closest_hour(self, time):
        hours = self.timedict.hours.unique()
        s,g = self.get_surrounding_ls(time.hour,hours)
        ds,dg = [abs(x-time.hour) for x in [s,g]]
        return [s,g][[ds,dg].index(min(ds,dg))]


    def _get_smaller_second(self, time):
        if type(time) is datetime.time:
            time = datetime.datetime.combine(datetime.datetime.date(),time)
        day = time.day

        def _get_greatest(time):
            seconds = self.timedict.loc[(self.timedict.hour == time.hour) & (self.timedict.minute == time.minute)].second.unique()
            try:
                return seconds[-1]
            except IndexError:
                return False

        seconds = self.timedict.loc[(self.timedict.hour == time.hour) & (self.timedict.minute == time.minute)].second.unique()
        s = self.get_lesser_ls(time.second, seconds)
        i = 0
        while not s:
            if i > self.TIMEOUT:
                print('TIMOUT FINDING SMALLER SECOND')
                exit()
            time = time - datetime.timedelta(seconds=time.second+1)
            if time.day != day:
                time.day = day
            s = _get_greatest(time)
        return s

    def _get_larger_second(self, time):
        if type(time) is datetime.time:
            time = datetime.datetime.combine(datetime.datetime.date(),time)
        day = time.day

        def _get_least(time):
            seconds = self.timedict.loc[(self.timedict.hour == time.hour) & (self.timedict.minute == time.minute)].second.unique()
            try:
                return seconds[0]
            except IndexError:
                return False

        seconds = self.timedict.loc[(self.timedict.hour == time.hour) & (self.timedict.minute == time.minute)].second.unique()
        s = self.get_greater_ls(time.second, seconds)
        i = 0
        while not s:
            if i > self.TIMEOUT:
                print('TIMOUT FINDING SMALLER SECOND')
                exit()
            time = time + datetime.timedelta(seconds=60 - time.second)
            if time.day != day:
                time.day = day
            s = _get_least(time)
        return s

def test_dict(tt):
    t1 = datetime.datetime.strptime('2020-08-17 12:58:00', "%Y-%m-%d %H:%M:%S")
    closest = tt.get_closest(t1)
    print(closest)

if __name__ == "__main__":
    of_path = '/mnt/c/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.out'
    from swmm_api import read_out_file
    import pandas as pd
    with read_out_file(of_path) as out:
        df = out.to_frame()
        df = df.xs(('link', 'Flow_velocity'), axis=1, level=[0, 2])
        timestamps = pd.to_datetime(df.index.values)
        tt = Timetree(timestamps)

    test_dict(tt)
    print('Module finished')