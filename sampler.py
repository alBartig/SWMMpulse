import pandas as pd
import datetime as dt
import numpy as np
from environment import Environment, DirectedTree, PACKET, CONSTITUENT, DirectedTree, HYDRAULICS, GROUP, DEFAULT
from prouter import Router

strategies = {"A":{"kind":"time",
                   "samplingfreq":f"{1}H",
                   "samplingduration":60,
                   "volume":250},
              "B":{"kind":"time",
                   "samplingfreq":f"{1/3}H",
                   "samplingduration":60,
                   "volume":250},
              "C":{"kind":"flow",
                   "samplingfreq":f"{1}H",
                   "samplingduration":60,
                   "volume":200},
              "D":{"kind":"flow",
                   "samplingfreq":f"{1/3}H",
                   "samplingduration":60,
                   "volume":200},
              "E":{"kind":"volume",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":50},
              "F":{"kind":"volume",
                   "samplecount":72,
                   "samplingduration":60,
                   "volume":250},
              "G":{"kind":"grab",
                   "samplingtime":dt.time(hour=9),
                   "samplingduration":120,
                   "volume":1000},
              "H":{"kind":"grab",
                   "samplingtime":dt.time(hour=12),
                   "samplingduration":120,
                   "volume":1000}}

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
        start = dt.date(year=2000, month=1, day=1)
        hours = float(freq.strip("H"))
        return pd.date_range(start, periods=24/hours, freq=freq)

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