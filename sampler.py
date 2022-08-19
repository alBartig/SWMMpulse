import pandas as pd
import datetime as dt
import numpy as np
from environment import Environment, DirectedTree, PACKET, CONSTITUENT, DirectedTree, HYDRAULICS, GROUP, DEFAULT
from prouter import Router

STRATEGIES = {"A":{"kind":"time",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=23, minute=59)},
              "B":{"kind":"time",
                   "samplecount":72,
                   "samplingduration":60,
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=23, minute=59)},
              "C":{"kind":"time",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":250,
                   "start-time":dt.time(hour=6, minute=0),
                   "end-time":dt.time(hour=12, minute=59)},
              "D":{"kind":"flow",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":200,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=23, minute=59)},
              "E":{"kind":"flow",
                   "samplecount":72,
                   "samplingduration":60,
                   "volume":200,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=23, minute=59)},
              "F":{"kind":"flow",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":200,
                   "start-time":dt.time(hour=6, minute=0),
                   "end-time":dt.time(hour=12, minute=59)},
              "G":{"kind":"volume",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":50,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=23, minute=59)},
              "H":{"kind":"volume",
                   "samplecount":72,
                   "samplingduration":60,
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=23, minute=59)},
              "I":{"kind":"volume",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":250,
                   "start-time":dt.time(hour=6, minute=0),
                   "end-time":dt.time(hour=12, minute=59)},
              "J":{"kind":"grab",
                   "samplingtime":dt.time(hour=9),
                   "samplingduration":120,
                   "volume":1000},
              "K":{"kind":"grab",
                   "samplingtime":dt.time(hour=12),
                   "samplingduration":120,
                   "volume":1000}}


class STRATEGY:
    SAMPLINGTIME = "samplingtime"
    START = "start-time"
    END = "end-time"
    SAMPLINGFREQ = "samplingfreq"
    SAMPLINGDURATION = "samplingduration"
    SAMPLINGVOLUME = "samplingvolume"
    SAMPLECOUNT = "samplecount"
    WEIGHTING = "weighting"
    WEIGHTS = "weights"
    GRAB = "grab"
    COMPOSITE = "composite"

class Sampler:
    def __init__(self, strategies=None):
        if strategies is None:
            strategies = strategies
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
        dt_index = pd.date_range("01-01-2000", periods=8640, freq="10S")
        self.flows = pd.Series(flow, index=dt_index)

    @staticmethod
    def get_sampling_hours(start, end, freq):
        """
        Creates a date_range with the points in time that are being sampled
        Args:
            freq (str): Frequency-string

        Returns: pd.date_range

        """
        start = dt.datetime.combine(dt.date(year=2000, month=1, day=1), start)
        end = dt.datetime.combine(dt.date(year=2000, month=1, day=1), end)
        return pd.date_range(start, end, freq=freq)

    @staticmethod
    def get_sampling_times(starttime, duration, freq="10S"):
        """
        returns the sampling interval for a given starttime
        Args:
            starttime (timestamp): 
            duration (dt.timedelta): 
            interval (str): 

        Returns:

        """
        if type(duration) == pd._libs.tslibs.timedeltas.Timedelta or type(duration) == dt.timedelta:
            nsteps = duration / pd.to_timedelta(freq)
        elif type(duration) == int:
            nsteps = dt.timedelta(seconds=duration) / pd.to_timedelta(freq)
        else:
            print(f"duration: {pd.to_timedelta(freq)}, type: {pd.to_timedelta(freq)}\n, should be: Timedelta"
                  f"freq: {freq}, type: {type(freq)}, should be: Timedelta")
            raise TypeError()
        return pd.date_range(starttime, periods=nsteps, freq=freq)

    @staticmethod
    def sampling_index_time(start, end, samplingfreq, duration=dt.timedelta(seconds=60), freq="10S"):
        """
        Returns a indexer for all sampled timesteps so that the correct timesteps can be picked from a given timeseries
        Args:
            samplingfreq (str): Frequency by which the samples is to sample 
            duration (dt.timedelta): Timedelta, over which the sample is taken 
            freq (str): Frequency of the timeseries 

        Returns:
            pd.DatetimeIndex
        """
        starttimes = Sampler.get_sampling_hours(start, end, samplingfreq)
        sampling_index = []

        if type(duration) is not dt.timedelta:
            duration = dt.timedelta(seconds=duration)

        for time in starttimes:
            sampling_index += Sampler.get_sampling_times(time, duration, freq).to_list()
        return pd.DatetimeIndex(sampling_index)

    def sampling_index_volume(self, start, end, samplecount=24, duration=dt.timedelta(seconds=120), freq="10S"):
        """
        Creates a volumne-weighted sampling index
        Args:
            samplecount (int): Number of samples to be taken over 24 hours
            duration (dt.timedelta): timedelta, over which sample is to be taken 
            freq (str): Frequency of the timeseries

        Returns:
            pd.DatetimeIndex
        """
        start = dt.datetime.combine(dt.date(year=2000, month=1, day=1), start)
        end = dt.datetime.combine(dt.date(year=2000, month=1, day=1), end)

        flows = self.flows[start:end]

        date = flows.index[0].date()
        starttimes = (flows.cumsum()//(flows.sum()/(samplecount-1))).drop_duplicates().index
        sampling_index = []
        for time in starttimes:
            sampling_index += Sampler.get_sampling_times(time, duration, freq).to_list()
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
        dt_index = pd.date_range("01-01-2000", periods=8640, freq="10S")
        df_timeseries.set_index(dt_index)

        kind = strategy.get("kind", "")
        if kind == "time":
            samplingfreq = f"{24 / strategy.get(STRATEGY.SAMPLECOUNT)}H"
            samplingindex = Sampler.sampling_index_time(start = strategy.get(STRATEGY.START),
                                                        end = strategy.get(STRATEGY.END),
                                                        samplingfreq=samplingfreq,
                                                        duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.loc[samplingindex,:].mean()

        elif kind == "flow":
            samplingfreq = f"{24 / strategy.get(STRATEGY.SAMPLECOUNT)}H"
            samplingindex = Sampler.sampling_index_time(start = strategy.get(STRATEGY.START),
                                                        end = strategy.get(STRATEGY.END),
                                                        samplingfreq=samplingfreq,
                                                        duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.multiply(self.flows, axis="index")
            sample_concs = sample_concs.div(self.flows.mean(), axis="index").loc[samplingindex,:].mean()

        elif kind == "volume":
            samplingindex = self.sampling_index_volume(start = strategy.get(STRATEGY.START),
                                                        end = strategy.get(STRATEGY.END),
                                                        samplecount=strategy.get(STRATEGY.SAMPLECOUNT),
                                                        duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.loc[samplingindex,:].mean()

        elif kind == "grab":
            samplingtime = dt.datetime.combine(df_timeseries.index[0].date(), strategy.get(STRATEGY.SAMPLINGTIME))
            samplingindex = Sampler.get_sampling_times(starttime=samplingtime,
                                                       duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.loc[samplingindex,:].mean()

        else:
            print("strategy requires key 'kind'")
            return False

        sample_concs.rename("concentration", inplace=True)
        return sample_concs


def preparations():
    env = Environment()
    env.read_swmmoutfile(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp.out")
    graph = DirectedTree.from_swmm(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp.inp")
    node_data = pd.read_csv(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp/pop_node_data.csv")
    node_data = node_data.set_index("NAME").to_dict(orient="index")
    graph.add_nodevalues(node_data)
    env.add_graph(graph)
    router = Router()
    router.add_environment(env)
    packets = router.environment.get_packets()
    routetable = router.route(packets=packets)
    processed = router.postprocess(routetable, packets, DEFAULT.DEFAULT_CONSTITUENTS.get(CONSTITUENT.COV))
    return processed


def main():
    import matplotlib.pyplot as plt
    ncols = 100


    # creating test timeseries
    dtindex = pd.date_range("2000-01-01", periods=8640, freq="10S")

    dfdict = {i: np.cumsum(0.005*(np.random.normal(loc=-2+i/(ncols/4), scale=100, size=8640) - 0.5)) + 100\
              for i in range(ncols)}
    df_means = pd.DataFrame([(-2+i/(ncols/4)) for i in range(ncols)], columns=["means"])
    df_timeseries = pd.DataFrame(dfdict, index=dtindex)

    sampler = Sampler()
    flows = pd.Series(-0.5/(4320**2)*(np.arange(8640)-4320)**2 + 1.25, index=dtindex)
    sampler.add_flows(flows)

    fig, ax = plt.subplots(ncols=2, nrows=3, figsize=[12,12], facecolor="white")

    for k, (name, strategy) in enumerate(strategies.items()):
        if k < 6:
            i, j = k//2, k%2
            samples = sampler.sample(df_timeseries, strategy)
            temp = df_means.join(samples)

            temp.plot(x="means",y="concentration", kind="scatter", ax=ax[i, j])

            try:
                textstr = '\n'.join([f'freq={float(strategy.get("samplingfreq").strip("H")):.2f} H',
                                     f'dur={strategy.get("samplingtime")} s',
                                     f'corr={temp["means"].corr(temp["concentration"]):.2f}'])
            except:
                textstr = "..."

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