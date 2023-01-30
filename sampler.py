import pandas as pd
import datetime as dt
import numpy as np
from environment import Environment, DirectedTree, PACKET, CONSTITUENT, DirectedTree, HYDRAULICS, GROUP, DEFAULT
from prouter import Router
from matplotlib.patches import Rectangle, Circle
from matplotlib.ticker import AutoMinorLocator
from matplotlib.collections import PatchCollection
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

STRATEGIES = {"A":{"kind":"time",
                   "samplecount":24,
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0),
                   "sampledtime":dt.timedelta(minutes=24)},
              "B":{"kind":"flow",
                   "samplecount":24,
                   "sampledtime":dt.timedelta(minutes=24),
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0)},
              "C":{"kind":"volume",
                   "samplecount":24,
                   "sampledtime":dt.timedelta(minutes=24),
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0)},
              "D":{"kind":"time",
                   "samplecount":48,
                   "sampledtime":dt.timedelta(minutes=24),
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0)},
              "E":{"kind":"time",
                   "samplecount":72,
                   "sampledtime":dt.timedelta(minutes=24),
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0)},
              "F":{"kind":"time",
                   "samplecount":24,
                   "sampledtime":dt.timedelta(minutes=24),
                   "volume":250,
                   "start-time":dt.time(hour=4, minute=0),
                   "end-time":dt.time(hour=16, minute=0)},
              "G":{"kind":"time",
                   "samplecount":24,
                   "sampledtime":dt.timedelta(minutes=24),
                   "volume":250,
                   "start-time":dt.time(hour=5, minute=0),
                   "end-time":dt.time(hour=11, minute=0)},
              "H":{"kind":"time",
                   "samplecount":24,
                   "sampledtime":dt.timedelta(minutes=48),
                   "volume":250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time":dt.time(hour=0, minute=0)},
              "I":{"kind": "time",
                   "samplecount": 24,
                   "volume": 250,
                   "start-time":dt.time(hour=0, minute=0),
                   "end-time": dt.time(hour=0, minute=0),
                   "sampledtime": dt.timedelta(minutes=72)},
              "J":{"kind":"grab",
                   "samplecount":1,
                   "sampledtime":dt.timedelta(minutes=1),
                   "volume":1000,
                   "start-time":dt.time(hour=9, minute=0)},
              "K":{"kind":"grab",
                   "samplecount":1,
                   "start-time":dt.time(hour=12, minute=0),
                   "sampledtime":dt.timedelta(minutes=4)},
              "Y":{"kind": "time",
                   "samplecount": 72,
                   "volume": 250,
                   "start-time":dt.time(hour=5, minute=0),
                   "end-time": dt.time(hour=11, minute=0),
                   "sampledtime": dt.timedelta(minutes=24)}}


class STRATEGY:
    KIND = "kind"
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
    SAMPLETIMES = "sampletimes"
    SAMPLEWEIGHTS = "sampleweights"

class WEIGHTING:
    VOLUME = "volume"
    FLOW = "flow"
    TIME = "time"
    GRAB = "grab"

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


    def get_starttime_from_strategy(self, strategy):
        # get starttime of sampling strategy
        starttime = self.dt_from_time(strategy.get(STRATEGY.START, None))
        if starttime is None:
            print("no starttime given for strategy. Starttime of 0:00 used")
            starttime = pd.to_datetime("2000-01-01 00:00:00")
        return starttime


    def get_endtime_from_strategy(self, strategy):
        # get endtime of sampling strategy
        endtime = self.dt_from_time(strategy.get(STRATEGY.END, None))
        if endtime is None:
            print("no endtime given for strategy. Endtime of 24:00 used")
            endtime = pd.to_datetime("2000-01-02 00:00:00")
        # add 1 day if endtime equals starttime
        starttime = self.get_starttime_from_strategy(strategy)
        if endtime == starttime:
            endtime = endtime + dt.timedelta(days=1)
        return endtime


    def get_sampledtime_from_strategy(self, strategy):
        sampledtime = strategy.get(STRATEGY.SAMPLEDTIME, None)
        if sampledtime is None:
            print("no sampled time given for strategy. Value of 24 minutes used")
            sampledtime = dt.timedelta(minutes=24)
        elif type(sampledtime) == int or type(sampledtime) == float:
            sampledtime = dt.timedelta(minutes=sampledtime)
        return sampledtime


    def get_samplecount_from_strategy(self, strategy):
        samplecount = strategy.get(STRATEGY.SAMPLECOUNT, None)
        if samplecount is None:
            print("no sample count given for strategy. Value of 24 samples used")
            samplecount = 24
        return samplecount


    def get_samplingfreq_from_strategy(self, strategy):
        """
        calculates samplingfrequency from sampling window (end - start) and samplecount in 10S resolution
        Args:
            strategy (dict): must contain start, end, samplecount

        Returns:
            string
        """
        start = self.get_starttime_from_strategy(strategy)
        end = self.get_endtime_from_strategy(strategy)
        samplecount = strategy.get(STRATEGY.SAMPLECOUNT)
        # sampling time window in seconds
        window_hours = np.around((end - start).total_seconds(), 0)
        # calculate sampling frequency as sampled time window / number of samples
        # calculate sampling frequency in 10S steps
        samplingfreq = f"{10 * np.around((window_hours) / (10 * samplecount))}S"
        return samplingfreq
    
    
    def get_samplingduration_from_strategy(self, strategy):
        """
        checks strategy for sampling duration, if not calculated yet, it calculates it and returns it
        Args:
            strategy (dict): must contain: samplecount and sampledtime

        Returns:
            dt.timedelta
        """
        duration = strategy.get(STRATEGY.SAMPLINGDURATION)
        
        if type(duration) is not dt.timedelta:
            samplecount = strategy.get(STRATEGY.SAMPLECOUNT)
            sampledtime = strategy.get(STRATEGY.SAMPLEDTIME)

            if type(sampledtime) is not dt.timedelta:
                sampledtime = dt.timedelta(seconds=sampledtime)

            duration = sampledtime / samplecount
                
            strategy[STRATEGY.SAMPLINGDURATION] = duration
        
        return duration

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


    def calculate_sampleweights(self, strategy):
        """
        calculates sampleweights, by which to weight each sample in composite sample
        Args:
            strategy (dict): must contain kind, sampletimes

        Returns:
            pd.Series: weights of length of sampletimes with mean of 1
        """
        kind = strategy.get(STRATEGY.KIND)
        sampletimes = strategy.get(STRATEGY.SAMPLETIMES)

        if kind == WEIGHTING.FLOW:
            weights = self.flows.loc[sampletimes]
            weights = 1 + (weights - weights.mean()) / (weights.max() - weights.min())

        elif kind == WEIGHTING.TIME or kind == WEIGHTING.VOLUME:
            weights = pd.Series(np.ones(len(sampletimes)), name="weights")
        
        elif kind == WEIGHTING.GRAB:
            weights = pd.Series([1], name="weights")

        else:
            raise ValueError(f"strategy kind needs to be of different type than {type(kind)}")

        return weights


    def calculate_sampletimes(self, strategy):
        """
        calculates the times at which a new sample for the composite sample is taken
        Args:
            strategy (dict): must contain: kind, start, end, samplecount

        Returns:
            pd.DateRange: Daterange with starttimes of samples
        """
        kind = strategy.get(STRATEGY.KIND)
        start = self.get_starttime_from_strategy(strategy)
        end = self.get_endtime_from_strategy(strategy)
        samplecount = strategy.get(STRATEGY.SAMPLECOUNT)

        if kind == WEIGHTING.TIME or kind == WEIGHTING.FLOW:
            samplingfreq = self.get_samplingfreq_from_strategy(strategy)
            starttimes = pd.date_range(start, end, freq=samplingfreq, inclusive="left")

        elif kind == WEIGHTING.VOLUME:
            samplecount = self.get_samplecount_from_strategy(strategy)
            flows = self.flows[start:end]
            # date = flows.index[0].date()
            starttimes = (flows.cumsum()//(flows.sum()/(samplecount))).drop_duplicates().index[:-1]
            
        elif kind == WEIGHTING.GRAB:
            starttimes = pd.Series([start], name="sampletimes")

        else:
            raise ValueError(f"strategy kind needs to be of different type than {type(kind)}")

        # check if samplecount matches starttimes
        if samplecount is None:
            print("No samplecount was handed to sampling index function, so no check could be completed")
        else:
            if len(starttimes) != samplecount:
                raise ValueError("The resulting sampling index does not match the samplecount")

        return starttimes


    def calculate_samplingindex(self, strategy, timeindexfreq="10S"):
        """
        Calculates samplingindex from the sampletimes and sampleweights
        Args:
            strategy (dict): must contain sampletimes, sampleweights and duration
            timeindexfreq: timestep resolution of timeseries

        Returns:
            pd.DatetimeIndex: index which is sampled of timeseries
        """
        sampletimes = strategy.get(STRATEGY.SAMPLETIMES)
        sampleweights = strategy.get(STRATEGY.SAMPLEWEIGHTS)
        duration = self.get_samplingduration_from_strategy(strategy)

        sampling_index = []
        
        durations = sampleweights * duration
        for dur, time in zip(durations, sampletimes):
            sampling_index += Sampler.get_sampling_times(time, dur, timeindexfreq).to_list()

        return pd.DatetimeIndex(sampling_index)


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
        sampledtime = self.get_sampledtime_from_strategy(strategy)
        samplecount = self.get_samplecount_from_strategy(strategy)
        # get start and end time of sampling strategy
        starttime = self.get_starttime_from_strategy(strategy)
        endtime = self.get_endtime_from_strategy(strategy)

        # ----------------------------- calculate strategy properties ------------------------------------
        # sampling time window in hours
        window_hours = np.around((endtime-starttime).total_seconds() / 3600, 2)
        # calculate sampling frequency as sampled time window / number of samples
        samplingfreq = f"{10 * np.around((window_hours * 3600) / (10 * samplecount))}S" # calculate sampling frequency in 10S steps
        # calculate sampled time of one sample by dividing sampled time / number of samples
        strategy[STRATEGY.SAMPLINGDURATION] = sampledtime.total_seconds() / samplecount # the case of grab sample uses sampled time

        strategy[STRATEGY.SAMPLETIMES] = self.calculate_sampletimes(strategy)
        strategy[STRATEGY.SAMPLEWEIGHTS] = self.calculate_sampleweights(strategy)
        samplingindex = self.calculate_samplingindex(strategy)

        # if kind == "time":
        #     samplingindex = self.sampling_index_time(start=starttime,
        #                                                 end=endtime,
        #                                                 samplingfreq=samplingfreq,
        #                                                 duration=sampleduration)
        # elif kind == "flow":
        #     samplingindex = self.sampling_index_flow(start=starttime,
        #                                                 end=endtime,
        #                                                 samplingfreq=samplingfreq,
        #                                                 duration=sampleduration)
        # elif kind == "volume":
        #     samplingindex = self.sampling_index_volume(start=starttime,
        #                                                end=endtime,
        #                                                samplecount=strategy.get(STRATEGY.SAMPLECOUNT),
        #                                                duration=sampleduration)
        # elif kind == "grab":
        #     samplingindex = Sampler.get_sampling_times(starttime=starttime,
        #                                                duration=sampledtime)
        # else:
        #     raise ValueError("strategy requires key 'kind'")

        # check if samplingindex matches sampled time
        if abs(dt.timedelta(seconds=len(samplingindex) * 10) - sampledtime) > dt.timedelta(minutes=1):
            raise ValueError("calculated samplingindex is off more than 1 minute from supposed sampled time")

        return samplingindex, strategy


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

        samplingindex, smeta = self.sampling_index(strategy)
        sample_concs = df_timeseries.loc[samplingindex,:].mean()
        sample_concs.rename("concentration", inplace=True)
        return sample_concs, smeta
    
    def plot_strategy(self, strategy, ax=None, y0=0, hmax=1.0, h0=0.5, legend=True, scaling=1):
        plotting = 0 # variable to store whether plot is returned or axes is assigned
        if ax is None:
            fig, ax = plt.subplots()
            plotting = 1
            
        sampling_index, smeta = self.sampling_index(strategy)
        sampleweights = smeta.get(STRATEGY.SAMPLEWEIGHTS).values
        sampletimes = smeta.get(STRATEGY.SAMPLETIMES).values
        
        h0 = h0 * hmax * scaling

        for sampletime, sampleweight in zip(sampletimes, sampleweights):
            start = mdates.date2num(pd.to_datetime(sampletime))
            end = mdates.date2num(pd.to_datetime(sampletime) + dt.timedelta(seconds=120))
            width = end - start
            height = h0 * sampleweight
            rect = Rectangle((start, y0+hmax/2-height/2), width, height, color="lightcoral", zorder=10)
            ax.add_patch(rect)
            
        # assign date locator / formatter to the x-axis to get proper labels
        locator = mdates.AutoDateLocator(minticks=3)
        formatter = mdates.DateFormatter('%H:%M')#locator
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_minor_locator(AutoMinorLocator(n=3))
        ax.set(xlim=[mdates.date2num(pd.to_datetime("2000-01-01 00:00:00")),
                         mdates.date2num(pd.to_datetime("2000-01-02 00:00:00"))])
        # add strategy metadata
        if legend is True:
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
    
    fig = plt.figure(figsize=[8,5*4], constrained_layout=True)
    #fig.suptitle()
    subfigs = fig.subfigures(nrows=5, ncols=1)
    combs = [["A","B","C"],["A","D","E"],["A","F","G"],["A","H"],["A","I","J"]]
    for k, comb in enumerate(combs):
        axs = subfigs[k].subplot_mosaic([[0, 1, 2],
                                         [3, 3, 3]],
                                        gridspec_kw={'height_ratios': [3,1]})
        for i, sname in enumerate(comb):
            # get and visualize strategy
            strategy = STRATEGIES.get(sname)
            sampler.plot_strategy(strategy, axs[3], y0=i, legend=False)
            axs[3].set(ylim=[0,3], yticks=[0.5, 1.5, 2.5], yticklabels=["a)","b)","c)"])
            axs[3].yaxis.set_minor_locator(AutoMinorLocator(n=2))
            axs[3].grid(axis="y", which="minor")
            axs[3].grid(axis="x", which="major")

            # sample and plot samples
            samples = sampler.sample(df_timeseries, strategy)
            temp = df_means.join(samples)
            temp.plot(x="means",y="concentration", kind="scatter", ax=axs[i])

    plt.show()

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