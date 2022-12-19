import pandas as pd
import datetime as dt
import numpy as np
from environment import Environment, DirectedTree, PACKET, CONSTITUENT, DirectedTree, HYDRAULICS, GROUP, DEFAULT
from prouter import Router
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

STRATEGIES = {"A":{"kind":"time",
                   "samplecount":24,
                   "sampledtime":24,
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0),
                   "sampledtime:":dt.timedelta(minutes=24)},
              "B":{"kind":"time",
                   "samplecount":48,
                   "sampledtime":24,
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0)},
              "C":{"kind":"time",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":250,
                   "start-time":dt.time(hour=6, minute=0),
                   "end-time":dt.time(hour=13, minute=0)},
              "D":{"kind":"flow",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":200,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0)},
              "E":{"kind":"flow",
                   "samplecount":72,
                   "samplingduration":60,
                   "volume":200,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0)},
              "F":{"kind":"flow",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":200,
                   "start-time":dt.time(hour=6, minute=0),
                   "end-time":dt.time(hour=13, minute=0)},
              "G":{"kind":"volume",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":50,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0)},
              "H":{"kind":"volume",
                   "samplecount":72,
                   "samplingduration":60,
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0)},
              "I":{"kind":"volume",
                   "samplecount":24,
                   "samplingduration":60,
                   "volume":250,
                   "start-time":dt.time(hour=6, minute=0),
                   "end-time":dt.time(hour=13, minute=0)},
              "J":{"kind":"grab",
                   "start-time":dt.time(hour=9),
                   "sampledtime": 1,
                   "samplingduration":120,
                   "volume":1000},
              "K":{"kind":"grab",
                   "start-time":dt.time(hour=12),
                   "sampledtime": 1,
                   "samplingduration":120,
                   "volume":1000},
              "L":{"kind":"time",
                   "samplecount":72,
                   "sampledtime":24,
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0)}}


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
    SAMPLEDTIME = "sampledtime"

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


    def dt_from_time(self, time):
        try:
            return pd.to_datetime(f"2000-01-01 {time.hour:02d}:{time.minute:02d}:{time.second:02d}")
        except:
            return None


    @staticmethod
    def get_sampling_hours(start, end, freq):
        """
        Creates a date_range with the points in time that are being sampled
        Args:
            freq (str): Frequency-string

        Returns: pd.date_range

        """
        return pd.date_range(start, end, freq=freq, inclusive="left")


    @staticmethod
    def get_sampling_times(starttime, duration, timeindexfreq="10S"):
        """
        returns the sampling interval for a given starttime
        Args:
            starttime (timestamp): 
            duration (dt.timedelta): 
            interval (str): 

        Returns:

        """
        if type(duration) == int or type(duration) == float:
            duration = dt.timedelta(seconds=duration)
        nsteps = np.around(duration / pd.to_timedelta(timeindexfreq))
        return pd.date_range(starttime, periods=nsteps, freq=timeindexfreq)


    def sampling_index_time(self, start, end, samplingfreq, samplecount=None, duration=dt.timedelta(seconds=60), timeindexfreq="10S"):
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

        # check if samplecount matches starttimes
        if samplecount is None:
            print("No samplecount was handed to sampling index function, so no check could be completed")
        else:
            if len(starttimes) != samplecount:
                raise ValueError("The resulting sampling index does not match the samplecount")

        sampling_index = []

        if type(duration) is not dt.timedelta:
            duration = dt.timedelta(seconds=duration)

        for time in starttimes:
            sampling_index += Sampler.get_sampling_times(time, duration, timeindexfreq).to_list()

        return pd.DatetimeIndex(sampling_index)
            

    def sampling_index_flow(self, start, end, samplingfreq, samplecount=None, duration=dt.timedelta(seconds=60), timeindexfreq="10S"):
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

        # check if samplecount matches starttimes
        if samplecount is None:
            print("No samplecount was handed to sampling index function, so no check could be completed")
        else:
            if len(starttimes) != samplecount:
                raise ValueError("The resulting sampling index does not match the samplecount")

        weights = self.flows.loc[starttimes]
        weights = 1 + (weights - weights.mean()) / (weights.max() - weights.min())
        sampling_index = []

        if type(duration) is not dt.timedelta:
            duration = dt.timedelta(seconds=duration)

        durations = weights * duration
        for dur, time in zip(durations, starttimes):
            sampling_index += Sampler.get_sampling_times(time, dur, timeindexfreq).to_list()
        
        return pd.DatetimeIndex(sampling_index)
    

    def sampling_index_volume(self, start, end, samplecount, duration=dt.timedelta(seconds=120), freq="10S"):
        """
        Creates a volumne-weighted sampling index
        Args:
            samplecount (int): Number of samples to be taken over 24 hours
            duration (dt.timedelta): timedelta, over which sample is to be taken 
            freq (str): Frequency of the timeseries

        Returns:
            pd.DatetimeIndex
        """
        #start = dt.datetime.combine(dt.date(year=2000, month=1, day=1), start)
        #end = dt.datetime.combine(dt.date(year=2000, month=1, day=1), end)

        flows = self.flows[start:end]

        date = flows.index[0].date()
        starttimes = (flows.cumsum()//(flows.sum()/(samplecount-1))).drop_duplicates().index

        # check if samplecount matches starttimes
        if samplecount is None:
            print("No samplecount was handed to sampling index function, so no check could be completed")
        else:
            if len(starttimes) != samplecount:
                raise ValueError("The resulting sampling index does not match the samplecount")

        sampling_index = []
        for time in starttimes:
            sampling_index += Sampler.get_sampling_times(time, duration, freq).to_list()
        sampling_index = pd.DatetimeIndex(sampling_index)
        return sampling_index[sampling_index.date == date]


    def sampling_index(self, strategy):
        dt_index = pd.date_range("01-01-2000", periods=8640, freq="10S")
        # ------------------------------ read strategy information ---------------------------------------
        # read kind of sampling strategy 
        kind = strategy.get("kind", "")
        sampledtime = strategy.get(STRATEGY.SAMPLEDTIME, None)
        if sampledtime is None:
            print("no sampled time given for strategy. Value of 24 minutes used")
            sampledtime = dt.timedelta(minutes=24)
        if type(sampledtime) == int:
            sampledtime = dt.timedelta(minutes=sampledtime)
        samplecount = strategy.get(STRATEGY.SAMPLECOUNT, None)
        if samplecount is None:
            print("no sample count given for strategy. Value of 24 samples used")
            samplecount = 24
        # get start and end time of sampling strategy
        starttime = self.dt_from_time(strategy.get(STRATEGY.START, None))
        if starttime is None:
            print("no starttime given for strategy. Starttime of 0:00 used")
            starttime = pd.to_datetime("2000-01-01 00:00:00")
        endtime = self.dt_from_time(strategy.get(STRATEGY.END, None))
        if endtime is None:
            print("no endtime given for strategy. Endtime of 24:00 used")
            endtime = pd.to_datetime("2000-01-02 00:00:00")
        # add 1 day if endtime equals starttime
        if endtime == starttime:
            endtime = endtime + dt.timedelta(days=1)

        # ----------------------------- calculate strategy properties ------------------------------------
        # sampling time window in hours
        window_hours = np.around((endtime-starttime).total_seconds() / 3600, 2)
        # calculate sampling frequency as sampled time window / number of samples
        samplingfreq = f"{10 * np.around((window_hours * 3600) / (10 * samplecount))}S" # calculate sampling frequency in 10S steps
        # calculate sampled time of one sample by dividing sampled time / number of samples
        sampleduration = sampledtime.total_seconds() / samplecount # the case of grab sample uses sampled time

        if kind == "time":
            samplingindex = self.sampling_index_time(start=starttime,
                                                        end=endtime,
                                                        samplingfreq=samplingfreq,
                                                        duration=sampleduration)
        elif kind == "flow":
            samplingindex = self.sampling_index_flow(start=starttime,
                                                        end=endtime,
                                                        samplingfreq=samplingfreq,
                                                        duration=sampleduration)
        elif kind == "volume":
            samplingindex = self.sampling_index_volume(start=starttime,
                                                       end=endtime,
                                                       samplecount=strategy.get(STRATEGY.SAMPLECOUNT),
                                                       duration=sampleduration)
        elif kind == "grab":
            samplingindex = Sampler.get_sampling_times(starttime=starttime,
                                                       duration=sampledtime)
        else:
            raise ValueError("strategy requires key 'kind'")
        # prepare metadata
        smeta = {STRATEGY.SAMPLECOUNT: samplecount,
                 STRATEGY.SAMPLEDTIME: sampledtime,
                 "kind": strategy["kind"],
                 STRATEGY.SAMPLINGDURATION: sampleduration,
                 STRATEGY.START: starttime,
                 STRATEGY.END: endtime}

        # check if samplingindex matches sampled time
        if abs(dt.timedelta(seconds=len(samplingindex) * 10) - sampledtime) > dt.timedelta(minutes=1):
            raise ValueError("calculated samplingindex is off more than 1 minute from supposed sampled time")

        return samplingindex, smeta


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

        samplingindex = self.sampling_index(strategy)

        if kind == "time":
            samplingfreq = f"{24 / strategy.get(STRATEGY.SAMPLECOUNT)}H"
            # samplingindex = Sampler.sampling_index_time(start = strategy.get(STRATEGY.START),
            #                                             end = strategy.get(STRATEGY.END),
            #                                             samplingfreq=samplingfreq,
            #                                             duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.loc[samplingindex,:].mean()

        elif kind == "flow":
            samplingfreq = f"{24 / strategy.get(STRATEGY.SAMPLECOUNT)}H"
            # samplingindex = Sampler.sampling_index_time(start = strategy.get(STRATEGY.START),
            #                                             end = strategy.get(STRATEGY.END),
            #                                             samplingfreq=samplingfreq,
            #                                             duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.multiply(self.flows, axis="index")
            sample_concs = sample_concs.div(self.flows.mean(), axis="index").loc[samplingindex,:].mean()

        elif kind == "volume":
            # samplingindex = self.sampling_index_volume(start = strategy.get(STRATEGY.START),
            #                                             end = strategy.get(STRATEGY.END),
            #                                             samplecount=strategy.get(STRATEGY.SAMPLECOUNT),
            #                                             duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.loc[samplingindex,:].mean()

        elif kind == "grab":
            # samplingtime = dt.datetime.combine(df_timeseries.index[0].date(), strategy.get(STRATEGY.SAMPLINGTIME))
            # samplingindex = Sampler.get_sampling_times(starttime=samplingtime,
            #                                            duration=strategy.get(STRATEGY.SAMPLINGDURATION))
            sample_concs = df_timeseries.loc[samplingindex,:].mean()

        else:
            print("strategy requires key 'kind'")
            return False

        sample_concs.rename("concentration", inplace=True)
        return sample_concs
    
    def plot_strategy(self, strategy, ax=None):
        plotting = 0 # variable to store whether plot is returned or axes is assigned
        if ax is None:
            fig, ax = plt.subplots()
            plotting = 1
            
        sampling_index, smeta = self.sampling_index(strategy)
        for indextime in sampling_index.values:
            start = mdates.date2num(pd.to_datetime(indextime))
            end = mdates.date2num(pd.to_datetime(indextime) + dt.timedelta(seconds=10))
            width = end - start
            rect = Rectangle((start, 0), width, 1, color="mistyrose")
            ax.add_patch(rect)
        # assign date locator / formatter to the x-axis to get proper labels
        locator = mdates.AutoDateLocator(minticks=3)
        formatter = mdates.AutoDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
        ax.set(xlim=[mdates.date2num(pd.to_datetime("2000-01-01 00:00:00")),
                         mdates.date2num(pd.to_datetime("2000-01-02 00:00:00"))],
                   ylim=[0, 1])
        # add strategy metadata
        textstr = f"{smeta['kind']} weighted\n" \
                  f"{'sampled time:':<15}{smeta[STRATEGY.SAMPLEDTIME].total_seconds()/60:>8.0f} min\n" \
                  f"{'sample count:':<15}{smeta[STRATEGY.SAMPLECOUNT]:>9d}\n" \
                  f"{'sample duration:':<15}{smeta[STRATEGY.SAMPLINGDURATION]:>6.0f} sec\n" \
                  f"{'window:':<15}{smeta[STRATEGY.START].strftime('%H:%M')} - {smeta[STRATEGY.END].strftime('%H:%M')}"
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='square', facecolor='mistyrose', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
                    verticalalignment='top', bbox=props, zorder=15)
        if plotting == 0:
            return ax
        else:
            plt.show()
            return None


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
    ncols = 100
    strategies = list(STRATEGIES.keys())
    n_strategies = len(strategies)
    # creating test timeseries
    dtindex = pd.date_range("2000-01-01", periods=8640, freq="10S")

    dfdict = {i: np.cumsum(0.005*(np.random.normal(loc=-2+i/(ncols/4), scale=100, size=8640) - 0.5)) + 100\
              for i in range(ncols)}
    df_means = pd.DataFrame([(-2+i/(ncols/4)) for i in range(ncols)], columns=["means"])
    df_timeseries = pd.DataFrame(dfdict, index=dtindex)

    sampler = Sampler()
    flows = pd.Series(-0.5/(4320**2)*(np.arange(8640)-4320)**2 + 1.25, index=dtindex)
    sampler.add_flows(flows)

    fig, axs = plt.subplots(ncols=1, nrows=n_strategies, figsize=[8,n_strategies*1.2], constrained_layout=True)

    for i, sname in enumerate(strategies):
        strategy = STRATEGIES.get(sname)
        sampler.plot_strategy(strategy, axs[i])
        # sampling_index, smeta = sampler.sampling_index(strategy)
        # for indextime in sampling_index.values:
        #     start = mdates.date2num(pd.to_datetime(indextime))
        #     end = mdates.date2num(pd.to_datetime(indextime) + dt.timedelta(seconds=10))
        #     width = end-start
        #     rect = Rectangle((start, 0), width, 1, color = "mistyrose")
        #     axs[i].add_patch(rect)
        # print(sampling_index.values)
        # # assign date locator / formatter to the x-axis to get proper labels
        # locator = mdates.AutoDateLocator(minticks=3)
        # formatter = mdates.AutoDateFormatter(locator)
        # axs[i].xaxis.set_major_locator(locator)
        # axs[i].xaxis.set_major_formatter(formatter)
        # axs[i].set(xlim=[mdates.date2num(pd.to_datetime("2000-01-01 00:00:00")),
        #                  mdates.date2num(pd.to_datetime("2000-01-02 00:00:00"))],
        #            ylim=[0,1])
        # # add strategy metadata
        # textstr = f"{smeta['kind']} weighted\n" \
        #           f"{'sampled time:':<15}{smeta[STRATEGY.SAMPLEDTIME]:>8.0f} min\n" \
        #           f"{'sample count:':<15}{smeta[STRATEGY.SAMPLECOUNT]:>9d}\n" \
        #           f"{'sample duration:':<15}{smeta[STRATEGY.SAMPLINGDURATION]:>6.0f} sec"
        # # these are matplotlib.patch.Patch properties
        # props = dict(boxstyle='square', facecolor='mistyrose', alpha=0.5)
        # # place a text box in upper left in axes coords
        # axs[i].text(0.05, 0.95, textstr, transform=axs[i].transAxes, fontsize=9,
        #         verticalalignment='top', bbox=props, zorder=15)
        #print("stop")

    # fig, ax = plt.subplots(ncols=2, nrows=3, figsize=[12,12], facecolor="white")

    # for k, (name, strategy) in enumerate(STRATEGIES.items()):
    #     if k < 6:
    #         i, j = k//2, k%2
    #         samples = sampler.sample(df_timeseries, strategy)
    #         temp = df_means.join(samples)
    #
    #         temp.plot(x="means",y="concentration", kind="scatter", ax=ax[i, j])
    #
    #         try:
    #             textstr = '\n'.join([f'freq={float(strategy.get("samplingfreq").strip("H")):.2f} H',
    #                                  f'dur={strategy.get("samplingtime")} s',
    #                                  f'corr={temp["means"].corr(temp["concentration"]):.2f}'])
    #         except:
    #             textstr = "..."
    #
    #         # these are matplotlib.patch.Patch properties
    #         props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #
    #         # place a text box in upper left in axes coords
    #         ax[i, j].text(0.05, 0.95, textstr, transform=ax[i, j].transAxes, fontsize=14,
    #                       verticalalignment='top', bbox=props)
    #         ax[i, j].set(title=f"strategy {name}", xlabel="fraction of infected")
    plt.show()
    print("stop")


if __name__ == "__main__":
    main()