import numpy as np
import pandas as pd
import datetime as dt
#from prouter import _Route_table, Postprocessing, preparation
#from pconstants import Loading

class Sampler:
    def __init__(self, duration=30, frequency="H", volume=0.25, resolution=10):

        #temporary start dates
        today = dt.datetime.today()
        start = dt.datetime.combine(today, dt.time(hour=0,minute=0,second=0))
        end = dt.datetime.combine(today, dt.time(hour=23,minute=59,second=59))

        self.settings = {'duration':duration, 'frequency':frequency, 'volume':volume, 'resolution':resolution}
        #temporary start times, date needs to be stripped later and replaced with date in given datetimeindex
        self.sample_times = pd.date_range(start, end, freq=frequency)

    @property
    def duration(self):
        return self.settings['duration']

    @property
    def frequency(self):
        return self.settings['frequency']

    @property
    def volume(self):
        return self.settings['volume']

    @property
    def resolution(self):
        return self.settings['resolution']

    def sample(self, df):
        if isinstance(df.columns, pd.MultiIndex):
            rv = self._sample_mi(df)
        else:
            rv = self.extract_series(self._sample_si(df))
        return rv

    def _sample_mi(self, midf):
        rd = {}
        i0 = midf.columns[0][:-1]
        df = self.extract_series(self._sample_si(midf.xs(i0, axis=1)), miprefix=i0)
        stepsize = len(midf.columns.levels[midf.columns.nlevels-1])
        for i in midf.columns[stepsize::stepsize]:
            miprefix = tuple(i[:-1])
            samples = self._sample_si(midf.xs(miprefix, axis=1))
            df = df.join(self.extract_series(samples, miprefix=miprefix))
        return df

    def _sample_si(self, sidf, as_df=True):
        date = sidf.index[0].date()
        n = np.floor(self.duration/self.resolution)
        v = self.volume/n
        samples = []
        for time in self.sample_times:
            t = dt.datetime.combine(date, time.time())
            slots = pd.date_range(t, periods=n, freq="10S")
            constituents = dict()
            for series in sidf:
                constituents[series] = np.mean([sidf[series][i] for i in slots])
            samples.append(Sample(time, self.volume, constituents))
        return samples

    @staticmethod
    def extract_series(samples, key="all", miprefix=()):
        if key == "all":
            constituents = set.union(*[set(sample.constituents.keys()) for sample in samples])
        else:
            constituents = [key]
        c = {const:[] for const in constituents}
        t = []
        for sample in samples:
            t.append(sample.time)
            [c[const].append(sample.constituents.get(const,0)) for const in constituents]
        mi = pd.MultiIndex.from_tuples([miprefix + (const,) for const in constituents])
        rdf = pd.DataFrame(c, index=t).columns=mi
        mi = pd.MultiIndex.from_tuples([miprefix + (const,) for const in constituents])
        rdf = pd.DataFrame(c, index=t)
        rdf.columns = mi
        return rdf


class Sample:
    def __init__(self, time, volume, constituents):
        self.time = time
        self.volume = volume
        self.constituents = constituents

    def __repr__(self):
        print(f"Sample volume: {self.volume}")
        print(f"Sample time:   {self.time}")
        print(f"Sample contents:")
        for key, value in self.constituents.items():
            print(f"{key}: {value}")

class Composite(Sample):
    def __init__(self, samples, weights=None):
        if weights is None:
            weights = [1] * len(samples)
        constituents = {c:[] for c in set.union(*[set(sample.constituents.keys()) for sample in samples])}
        time = min([sample.time for sample in samples])
        volume = samples[0].volume
        for sample in samples:
            for constituent, concentration in sample.constituents.items():
                constituents[constituent].append(concentration)
        for constituent, concentrations in constituents.items():
            constituents[constituent] = np.sum([c*w for c,w in list(zip(concentrations,weights))])/np.sum(weights)
        super().__init__(time, volume, constituents)


def load_pproc(qlut, graph, load=True):
    print('Test Postprocessing')
    evalnode = 'MH327-088-1'
    #routetable = _Route_table.from_file('C:/Users/alber/Documents/swmm/swmmpulse/route_tables/routetable.pickle')
    routetable = _Route_table.from_file('/mnt/c/Users/albert/Documents/SWMMpulse/route_tables/routetable.pickle')
    pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)

    #pproc.process_constituent(Loading.FECAL, entry_loc='C:/Users/alber/Documents/swmm/swmmpulse/entries/entries.pickle', load=load)
    pproc.process_constituent(Loading.FECAL, entry_loc='/mnt/c/Users/albert/Documents/SWMMpulse/entries/entries.pickle', load=load)
    print("Test postprocessing finished")
    return pproc

def test_sampler():
    qlut, graph = preparation()
    pproc = load_pproc(qlut, graph, True)
    sampler = Sampler(30,'H', 0.25)
    x,y = pproc.Fecal_Matter.timeseries()
    s = pd.Series(y,index=x,name='Fecal Matter')
    df = pd.DataFrame(s)
    samples = sampler.sample(df)
    print(sampler.settings)


def main():
    test_sampler()

if __name__ == "__main__":
    main()