import datetime as dt
import numpy as np
import random
import pandas as pd
from swmm_api import read_out_file


class UNITS:
    GRAM = 'g'
    MILLIGRAM = 'mg'
    MICROGRAM = 'ug'
    COUNT = '#'
    UNIT = "unit"
    DATE = "date"


class PACKET:
    CLASSIFICATION = "classification"
    PACKETID = "pid"
    CONSTITUENTS = "constituents"
    ORIGIN = "origin"
    T0 = "t0"
    AGE = "age"
    ARRIVAL_TIME = "arrival_time"
    LINK = "link"
    NODE = "node"


class CONSTITUENT:
    DECAY_RATE = "decay_rate"
    SPECIFIC_LOAD = "specific_load"
    DISPERSION_RATE = "dispersion_rate"
    NAME = "name"
    FECAL = 'Fecal-Matter'
    COV = 'Cov-RNA'
    PEP = 'Pepper-virus'


class GROUP:
    GROUPS = "groups"
    NAME = "name"
    DAILYPOOPS = "dailypoops"
    WEIGHT = "weight"
    CONSTITUENTS = "constituents"
    PATTERN = "pattern"


class DEFAULT:
    CONST_FECAL = {CONSTITUENT.NAME: CONSTITUENT.FECAL, CONSTITUENT.SPECIFIC_LOAD: 200,
                   UNITS.UNIT: UNITS.GRAM, CONSTITUENT.DECAY_RATE: 0.114}
    CONST_COV = {CONSTITUENT.NAME: CONSTITUENT.COV, CONSTITUENT.SPECIFIC_LOAD: 1000, UNITS.UNIT: UNITS.COUNT}
    CONST_PEP = {CONSTITUENT.NAME: CONSTITUENT.PEP, CONSTITUENT.SPECIFIC_LOAD: 1000, UNITS.UNIT: UNITS.COUNT}

    DEFAULT_CONSTITUENTS = [CONST_FECAL, CONST_COV]

    PATTERN_BRISTOL_MEN = [1.4, 0.3, 0.1, 0.0, 0.3, 1.7, 9.1, 21, 13, 9, 6.9, 4.9, 1.9, 3.6, 2.5, 2, 2.9, 2.3, 4.1, 4.0,
                           2.7,
                           2.1, 2.2, 2.0]
    PATTERN_BRISTOL_WOMEN = [2.7, 0.1, 0.1, 0.08, 0.02, 0.2, 3.3, 16.9, 20.1, 12.4, 8.8, 5.0, 2.5, 2.6, 3.1, 2.5, 3.0,
                             2.2, 4.3,
                             3.3, 2.1, 1.5, 2.0, 1.2]
    PATTERN_BRISTOL = [round((w * 0.5 + v * 0.5), 3) for w, v in zip(PATTERN_BRISTOL_MEN, PATTERN_BRISTOL_WOMEN)]

    GROUP_HEALTHY = {GROUP.NAME: "Healthy", GROUP.WEIGHT: 0.8, GROUP.DAILYPOOPS: 1,
                     GROUP.CONSTITUENTS: [CONST_FECAL], GROUP.PATTERN: PATTERN_BRISTOL}
    GROUP_INFECTED = {GROUP.NAME: "Infected", GROUP.WEIGHT: 0.2, GROUP.DAILYPOOPS: 1,
                      GROUP.CONSTITUENTS: [CONST_FECAL, CONST_COV], GROUP.PATTERN: PATTERN_BRISTOL}

    ENVIRONMENT = {
        UNITS.DATE: dt.date(day=1, month=1, year=2000),
        GROUP.GROUPS: [GROUP_HEALTHY, GROUP_INFECTED],
        CONSTITUENT.DISPERSION_RATE: 0.16
    }


class Environment:
    def __init__(self, information=None):
        if information is None:
            information = DEFAULT.ENVIRONMENT
        self.information = information
        self.time_range = pd.date_range(information.get(UNITS.DATE), periods=24 * 60 * 6, freq="10S")
        # standardize patterns in groups
        for group in self.information.get(GROUP.GROUPS):
            group[GROUP.PATTERN] = self._standardize_pattern(group.get(GROUP.PATTERN))
        
    def add_graph(self, graph):
        """
        Stores graph data within the environment
        Args:
            graph (DirectedTree: graph, swmm-model

        Returns:
            None
        """
        self.graph = graph
        return None

    def _standardize_pattern(self, pattern):
        """
        Takes the input weights and interpolates weights for a frequency of 10S
        Args:
            pattern (list): list of weights

        Returns:
            list
        """
        date = self.information.get(UNITS.DATE)
        pattern = list(pattern)  # convert pattern to list, may be ndarray instead
        pattern.append(pattern[0])  # end pattern as it started
        dest_range = pd.date_range(date, date + dt.timedelta(days=1),
                                   periods=len(self.time_range) + 1,
                                   inclusive="both")  # date_range in destination frequency
        orig_range = pd.date_range(date, date + dt.timedelta(days=1),
                                   periods=len(pattern),
                                   inclusive="both")  # date_range with period of input pattern
        s_weights = pd.Series(pattern, index=orig_range)
        s_weights = s_weights.reindex(dest_range)
        s_weights.interpolate(limit_direction="both", inplace=True)
        s_weights.drop(s_weights.index[-1], inplace=True)
        return s_weights.values

    def _packet_gen(self, i, group, node=None):
        """
        creates a packet dictionary with an id and an origin time. The time is assigned randomly
        according to the time pattern in the environment
        Args:
            i (int): index number for packet, should be unique in the list

        Returns:
            dict
        """
        packet = {PACKET.PACKETID: f"P{int(i):d}",
                  PACKET.T0: random.choices(self.time_range.values, weights=group.get(GROUP.PATTERN), k=1),
                  PACKET.CLASSIFICATION: group.get(GROUP.NAME),
                  PACKET.CONSTITUENTS: group.get(GROUP.CONSTITUENTS),
                  PACKET.ORIGIN: node}
        return packet

    def get_packets(self, population=None):
        """
        Returns a dictionary with number of packets of the given environment.
        A pool of individuals is created by creating a list with the node-strings in it, as often as there are people at
        that node. This pool is shuffled. For each group, the weighted slice is taken from the pool and then
        extended/shrunken according to the dailypoops of that group. Then a dictionary with all packets and their
        metadata is created for each poop.
        Args:
            population (list): list of tuples (node, population)

        Returns:
            dict
        """
        if population is None:
            population = [(node, valdict.get("POP",0)) for node, valdict in self.graph.adjls.items()]
        # create a list of strings (node-names) were each string represens an individual in the environment.
        # From the list, individuals are drawn at random later
        people = []
        ls_nodes = [node[0] for node in population] #list of nodes in graph
        ls_pop = [node[1] for node in population] #list of populations in graph
        pop_tot = sum(ls_pop)

        # calculate sum of weights of groups
        grp_weights_tot = sum(group.get(GROUP.WEIGHT) for group in self.information.get(GROUP.GROUPS))

        #prepare lists for packet attributes
        arr_t0s = [] #prepare list arrival times
        arr_class = [] #prepare list of classifications
        arr_contents = [] #prepare list of contents
        arr_nodes = [] #prepare list of origin nodes

        #assign attributes for each group
        groups = self.information.get(GROUP.GROUPS)
        for group in groups:
            # calculate number of people to pick
            dailypoops = group.get(GROUP.DAILYPOOPS)
            number = int(round(dailypoops * group.get(GROUP.WEIGHT) * pop_tot / grp_weights_tot, 0))
            arr_t0s += random.choices(self.time_range.values, weights=group.get(GROUP.PATTERN), k=number)
            arr_class += [group.get(GROUP.NAME)] * number
            arr_contents = [group.get(GROUP.CONSTITUENTS)] * number
            arr_nodes += random.choices(ls_nodes, weights=ls_pop, k=number)
        arr_pids = [f"P{i}" for i in range(len(arr_t0s))]  #counter for packets

        packets = pd.DataFrame(data=[arr_pids, arr_class, arr_nodes, arr_t0s, arr_contents],
                               index=[PACKET.PACKETID, PACKET.CLASSIFICATION, PACKET.ORIGIN, PACKET.T0, PACKET.CONSTITUENTS]).T
        return packets

    def read_swmmoutfile(self, outfile_path):
        """
        Reads swmm-outfile, interpolates to 10s frequency and reindexes to chosen environment date
        Args:
            outfile_path (str): path to swmm-outfile

        Returns:
            None
        """
        with read_out_file(outfile_path) as out:
            # import outfile as dataframe
            df = out.to_frame()
            df = df.loc[:, ("link", slice(None), ["flow", "velocity"])].droplevel(0, axis=1)
            # reindex and interpolate dataframe
            _date = df.index[0].date()
            _dtindex = pd.date_range(_date, periods=8641, freq="10S")
            df = df.reindex(_dtindex)
            df.loc[df.index[-1], :] = df.loc[df.index[0], :]
            df.interpolate(direction="both", inplace=True)
            df.drop(df.index[-1], inplace=True)
            date = self.information.get(UNITS.DATE)
            df.index = df.index.map(lambda d: d.replace(year=date.year, month=date.month, day=date.day))
            self.flow_rates = df.xs("flow", axis=1, level=1).to_dict(orient="list")
            self.flow_velocities = df.xs("velocity", axis=1, level=1).to_dict(orient="list")
            self.TIMESTEPS = len(df.index)
        return None

    def _calc_index(self, time):
        hr, mn, sk = time.hour, time.minute, time.second  # calculate index of tm in timeseries
        i0 = (hr * 360 + mn * 6 + sk) % self.TIMESTEPS
        return i0

    def get_velocity(self, time, link):
        return self.flow_velocities.get(link)[self._calc_index(time)]

    def get_flow_rate(self, time, link):
        return self.flow_rates.get(link)[self._calc_index(time)]

    def read_swmminpfile(self, inpfile_path):
        self.graph = DirectedTree.from_swmm(inpfile_path)
        return None


class DirectedTree:
    def __init__(self, links, nodes=None):
        self.nodeindex = 0
        self.linkindex = 0
        self.links = {}
        self.adjls = {}

        for link in links:
            self.add_link(link)

        if nodes is not None:
            for n, p in nodes:
                self.add_nodevalues(n, [('name', n), ('population', p)])

    @classmethod
    def from_swmm(cls, path):
        from swmm_api.input_file.section_labels import JUNCTIONS, OUTFALLS, STORAGE, ORIFICES, DIVIDERS,\
            CONDUITS, WEIRS, PUMPS, OUTLETS
        from swmm_api import read_inp_file

        linktypes = [CONDUITS, PUMPS, ORIFICES, WEIRS]
        nodetypes = [JUNCTIONS, OUTFALLS, STORAGE, DIVIDERS]
        inp = read_inp_file(path)
        links = []
        for linktype in linktypes:
            if linktype in inp:
                try:
                    links += list(
                        inp[linktype].get_dataframe(set_index=False)[["Name", "FromNode", "ToNode", "Length"]].values)
                except:
                    links += list(inp[linktype].get_dataframe(set_index=False)[["Name", "FromNode", "ToNode"]].values)
        for i in range(len(links)):
            try:
                links[i][3] = ("length", links[i][3])
            except:
                pass
        return cls(links)

    def add_link(self, link):
        name = link[0]
        inlet = link[1]
        outlet = link[2]

        if self.links.__contains__(name) == False:
            self.links[name] = {'name': name, 'inletnode': inlet,
                                'outletnode': outlet, 'linkindex': self.linkindex}
            self.linkindex += 1

            if len(link) > 3:
                for attribute in link[3:]:
                    field = attribute[0]
                    value = attribute[1]
                    self.links[name][field] = value

            if self.adjls.__contains__(inlet):
                self.adjls[inlet]['outlets'].append([outlet, name])
            else:
                self.adjls[inlet] = {'inlets': [], 'outlets': [[outlet, name]],
                                     'nodeindex': self.nodeindex}
                self.nodeindex += 1

            if self.adjls.__contains__(outlet):
                self.adjls[outlet]['inlets'].append([inlet, name])
            else:
                self.adjls[outlet] = {'inlets': [[inlet, name]], 'outlets': [],
                                      'nodeindex': self.nodeindex}
                self.nodeindex += 1

    def add_nodevalue(self, node, field, value, overwrite=False):
        if self.adjls.__contains__(node):
            if overwrite == False and self.adjls[node].__contains__(field):
                self.adjls[node][field] += value
            else:
                self.adjls[node][field] = value
            return True
        else:
            return False

    def add_nodevalues(self, dictionary):
        """
        Takes a dictionary {node: {key: value},...} to add values to nodes
        Args:
            dictionary (dict):

        Returns:
            None
        """
        for node, valdict in dictionary.items():
            try:
                self.adjls[node].update(valdict)
            except:
                print(f"node {node} not found in graph")
        return None

    def add_linkvalue(self, link, field, value):
        if self.links.__contains__(link):
            self.links[link][field] = value
            return True
        else:
            return False

    def check_node(self,node):
        if self.adjls.__contains__(node):
            return True
        else:
            return False

    def get_inletnodes(self,node):
        if self.check_node(node):
            inlets = self.adjls[node]['inlets']
            inletnodes = [inlet[0] for inlet in inlets]
            return inletnodes
        else:
            return False

    def get_nodevalue(self, node, value):
        if self.adjls.__contains__(node):
            if self.adjls[node].get(value) != None:
                return self.adjls[node].get(value)
            else:
                return False
        else:
            return False

    def get_nodeindex(self, node):
        return self.adjls[node]['nodeindex']
    
    def order_shreve(self):
        """
        assigns Shreve order values to graph-nodes
        Returns:
            None
        """
        import numpy as np

        def dfs(junction,value,target):
            #import list of visited nodes
            nonlocal visited
            #prepare array for upstream acc_values
            upstream_values = []

            for node in self.get_inletnodes(junction):
                at = self.get_nodeindex(node)
                if visited[at] != True:
                    visited[at] = True
                    dfs(node,value,target)
                    upstream_values.append(self.get_nodevalue(node,target))

            #Check if node has value assigned
            try:
                local_acc = sum(upstream_values)+self.get_nodevalue(junction,value)
            except:
                local_acc = sum(upstream_values)

            self.add_nodevalue(junction,target,local_acc)

        if end is None:
            end = self.root

        visited = np.zeros(self.nodeindex+1,dtype=bool)
        #name target field
        target = "shreve"
        dfs(end,value,target)
        return True


def test_flows():
    env = Environment()
    env.read_swmmoutfile(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp.out")
    print("finished reading output file")
    graph = DirectedTree.from_swmm(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp.inp")
    node_data = pd.read_csv(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp/pop_node_data.csv")
    node_data = node_data.set_index("NAME").to_dict(orient="index")
    graph.add_nodevalues(node_data)
    env.add_graph(graph)
    print("finished reading input file")
    env.get_packets()
    print("finished testing sequence")


def main():
    test_flows()


if __name__ == '__main__':
    main()
