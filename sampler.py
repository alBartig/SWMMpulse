import pandas as pd
import datetime as dt
import numpy as np

strategies = {"A":{"samplingfreq":"1H",
                   "samplingtime":60,
                   "volume":250,
                   "details":"standard"},
              "B":{"samplingfreq":f"{5/3}H",
                   "samplingtime":60,
                   "volume":250},
              "C":{"samplingfreq":f"{1/4}H",
                   "samplingtime":60,
                   "volume":200},
              "D":{"samplingfreq":f"{1/3}H",
                   "samplingtime":180,
                   "volume":50},
              "E":{"samplingfreq":f"{5/3}H",
                   "samplingtime":60,
                   "volume":250},
              "F":{"samplingfreq":f"{5/3}H",
                   "samplingtime":60,
                   "volume":250},
              "G":{"samplingfreq":f"{5/3}H",
                   "samplingtime":60,
                   "volume":250}}

class STRATEGY:
    SAMPLINGTIME = "samplingtime"
    SAMPLINGFREQ = "samplingfreq"
    SAMPLINGDURATION = "samplingduration"
    SAMPLINGVOLUME = "samplingvolume"
    SAMPLECOUNT = "samplecount"
    WEIGHTING = "weighting"
    WEIGHTS = "weights"
    GRAB = "grab"
    COMPOSITE = "composite"

class Sampler:
    def __init__(self, strategies=strategies):
        self.strategies = strategies

    def add_strategies(self, strategies):
        """
        Adds strategies to the sampling strategy dict.
        "strategies" is a dictionary containing the strategies and
        the names as key value pairs:
        {key1: strategy1, key2: strategy2}

        strategies must contain:
            "kind":             "time"/"flow"/"volume"/"grab"
            "samplingtime":     dt.time, time the grab sample is taken
            "samplingduration": int, duration in seconds, that the sample is taken over
            "samplingfreq":     float, frequency in hours, the samples are taken
            "samplingvolume":   In volume weighted sampling, the volume of flow delimiting each sample
            "flow":             pd.Series, flow

        Args:
            strategies (dict): dictionary of dictionaries containing the strategies

        Returns: None

        """
        for key, strategy in strategies:
            self.strategies[key] = strategy

    def add_flows(self, flow):
        """
        Adds flows to the sampler for flow- and volume-weighted sampling
        Args:
            flow (pd.Series): Series containing the flows with same index as concentrations

        Returns: None

        """
        self.flows = flow

    @staticmethod
    def get_sampling_hours(freq):
        start = dt.datetime.today().date()
        end = pd.to_datetime(dt.datetime.today()).ceil("D")
        return pd.date_range(start, end, freq=freq, closed="left")

    @staticmethod
    def get_sampling_times(starttime, duration, interval=10):
        nsteps = duration/interval
        interval = f"{interval}S"
        return pd.date_range(starttime, periods=nsteps, freq=interval)

    @staticmethod
    def sampling_index_time(samplingfreq, duration=120, interval=10):
        starttimes = Sampler.get_sampling_hours(samplingfreq)
        sampling_index = []
        for time in starttimes:
            sampling_index += Sampler.get_sampling_times(time, duration, interval).to_list()
        return pd.DatetimeIndex(sampling_index)

    def sampling_index_volume(self, samplecount=24, duration=120, interval=10):
        flows = self.flows
        date = flows.index[0].date()
        starttimes = (flows.cumsum()//(flows.sum()/(samplecount-1))).drop_duplicates().index
        sampling_index = []
        for time in starttimes:
            sampling_index += Sampler.get_sampling_times(time, duration, interval).to_list()
        sampling_index = pd.DatetimeIndex(sampling_index)
        return sampling_index[sampling_index.date == date]

    def sample(self, df_timeseries, strategy):
        """
        Takes a dataframe with timeseries (index=datetimeindex, columes=measurements) and samples each column
        according to strategy, returns series with mean sample concentrations
        Args:
            df_timeseries (pd.DataFrame): DataFrame of timeseries
            strategy (dict): dictionary containing {"samplingfreq", "samplingtime", "volume"}
                            a dict containing sample strategies can be loaded with "from sampler import strategies

        Returns:
            pd.Series: with mean sample concentrations for each timeseries

        """
        kind = strategy.get("kind", "time")
        if kind == "time":
            samplingindex = Sampler.sampling_index_time(samplingfreq=strategy.get(STRATEGY.SAMPLINGFREQ),
                                                        duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.loc[samplingindex,:].mean()

        elif kind == "flow":
            samplingindex = Sampler.sampling_index_time(samplingfreq=strategy.get(STRATEGY.SAMPLINGFREQ),
                                                        duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.multiply(self.flows, axis="index")
            sample_concs = sample_concs.div(self.flows.mean(), axis="index").loc[samplingindex,:].mean()

        elif kind == "volume":
            samplingindex = self.sampling_index_volume(samplecount=strategy.get(STRATEGY.SAMPLECOUNT),
                                                        duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.loc[samplingindex,:].mean()

        elif kind == "grab":
            samplingtime = dt.datetime.combine(df_timeseries.index[0].date(), strategy.get(STRATEGY.SAMPLINGTIME))
            samplingindex = Sampler.get_sampling_times(starttime=samplingtime,
                                                       duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.loc[samplingindex,:].mean()

        sample_concs.rename("concentration", inplace=True)
        return sample_concs

def main():
    import matplotlib.pyplot as plt
    ncols = 100

    # creating test timeseries
    dtindex = pd.date_range(dt.datetime.today().date(), periods=8640, freq="10S")

    dfdict = {i: np.cumsum(0.005*(np.random.normal(loc=-2+i/(ncols/4), scale=100, size=8640) - 0.5)) + 100\
              for i in range(ncols)}
    df_means = pd.DataFrame([(-2+i/(ncols/4)) for i in range(ncols)], columns=["means"])
    df_timeseries = pd.DataFrame(dfdict, index=dtindex)

    fig, ax = plt.subplots(ncols=2, nrows=3, figsize=[12,12], facecolor="white")

    for k, (name, strategy) in enumerate(strategies.items()):
        if k < 6:
            i, j = k//2, k%2
            samples = Sampler.sample(df_timeseries, strategy)
            temp = df_means.join(samples)

            temp.plot(x="means",y="concentration", kind="scatter", ax=ax[i, j])

            textstr = '\n'.join([f'freq={float(strategy.get("samplingfreq").strip("H")):.2f} H',
                                 f'dur={strategy.get("samplingtime")} s',
                                 f'corr={temp["means"].corr(temp["concentration"]):.2f}'])

            # these are matplotlib.patch.Patch properties
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

            # place a text box in upper left in axes coords
            ax[i, j].text(0.05, 0.95, textstr, transform=ax[i, j].transAxes, fontsize=14,
                          verticalalignment='top', bbox=props)
            ax[i, j].set(title=f"strategy {name}", xlabel="fraction of infected")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()