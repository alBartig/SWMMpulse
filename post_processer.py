from pconstants import Discretization
from pconstants import Loading
import copy
import datetime

class Postprocessing:
    def __init__(self, routed_node, load_info, tlookup):
        self.constituents = load_info[Loading.CONSTITUENTS]
        self.dispersion = load_info[Loading.DISPERSION]
        self.packets = routed_node
        self.tlookup = tlookup

    def process_arrivals(self):
        for packet,ta in self.packets:
            packet.calculate_age(ta)

    def calculate_constituent_masses(self, constituent, graph, tlookup):
        link = graph.get_nodeinlet(node)
        ts = np.zeros((Discretization.SERIESLENGTH / Discretization.TIMESTEPLENGTH))
        for packet in self.packets:
            ts += packet.timeseries_masses(tlookup, link)

    def calculate_constituent_concentrations(self, constituent):
        #load load information
        spec_load = self.load_info[constituent][Loading.LOAD]
        spec_degr = self.load_info[constituent][Loading.DEGRADATION]
        spec_classes = self.load_info[constituent][Loading.CLASSIFICATIONS]
        #calculate loads
        packets = copy.copy(self.packets)
        [packet.calculate_load(spec_load,spec_degr) for packet in packets]

        tssize = int(Discretization.SERIESLENGTH / Discretization.TIMESTEPLENGTH)
        tstime = np.array
        #tsconcentrations =


        self.timeseries = timeseries

if __name__ == "__main__":
    start = datetime.datetime(year=2021,month=4,day=1,hour=10,minute=0,second=0)
    end = datetime.datetime(year=2021,month=4,day=1,hour=11)
    diff = end - start
    size = int(diff/Discretization.TIMESTEPLENGTH + 1)
    timesteps = [start+datetime.timedelta(seconds=step*Discretization.TIMESTEPLENGTH) for step in range(size)]
