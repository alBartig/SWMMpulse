import pandas as pd
import datetime as dt

class Sampler:
    @staticmethod
    def get_sampling_hours(freq):
        start = dt.datetime.today().date()
        end = pd.to_datetime(dt.datetime.today()).ceil("D")
        return pd.date_range(start, end, freq=freq, closed="left")

    @staticmethod
    def get_sampling_times(starttime, duration, interval):
        nsteps = duration/interval
        interval = f"{interval}S"
        return pd.date_range(starttime, periods=nsteps, freq=interval)

    @staticmethod
    def build_sampling_index(samplingfreq, duration=120, interval=10):
        starttimes = Sampler.get_sampling_hours(samplingfreq)
        sampling_index = []
        for time in starttimes:
            sampling_index += Sampler.get_sampling_times(time, duration, interval).to_list()
        return pd.DatetimeIndex(sampling_index)