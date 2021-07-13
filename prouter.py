import geopandas
from network_lib import DirectedTree as ntwk, GraphConstants as gc
from timeseries import QSeries, TSeries, TDict
import create_loadings
from pconstants import Discretization, Loading
import datetime
from environment import Environment
from tqdm import tqdm
import pickle
from Exceptions import RoutingError, PlausibilityError
import pandas as pd
import numpy as np

class PRouter:
    def __init__(self,graph,qlookup):
        self.qlookup = qlookup
        self.graph = graph
        self.env = Environment()

    def route_packet(self, packet, path):
        ct = packet.t0
        stops = [(packet.origin,ct)]
        for stop in path[1:]:
            node = stop[0]
            link = stop[1]
            #lookup velocity, if v == 0 break routing
            v,delay = self.qlookup.lookup_v(link, ct)
            if v is None:
                break
            try:
                td = int(round(self.graph.get_linkvalue(link,gc.LENGTH) / v))
            except:
                print(f"Link: {link}, Conduit length: {self.graph.get_linkvalue(link,'LENGTH')}, Velocity: {v}")
                raise ValueError
            if td < 0:
                print('flowtime negative')
                raise PlausibilityError
            ct = ct + datetime.timedelta(seconds=td) + delay
            stops.append((node,ct))
            #check for plausibility:
            for stop in stops:
                if stop[1] < packet.t0:
                    raise RoutingError(packet, stops)
        return (packet,stops)

    def route_plist(self,plist, path):
        return [self.route_packet(p, path) for p in plist]

    def node_plist(self, nname, npop):
        """
        Creates packetlist (plist) for node
        Args:
            node (dict): {'name':nodename,'population':overall population at node}

        Returns:
            list: list with pobjects for each group in self.env
        """
        plist = []
        [plist.extend(group.plist(nname,npop,self.env.date)) for group in self.env.groups]
        #print('Packetlist created successfully')
        return plist

    def node_path(self, node):
        '''
        function takes node, generates seedlist via stored community info and returns dictionary with a list of stops
        and the a list with the routed packets
        :param node: str, node that is to be simulated
        :return: dict, {path:list of nodes, packets:list of tuples (packet, list with arrival times)}
        '''
        ppath = self.graph.trace_path(node)
        #print('Path created successfully')
        return ppath

    def route(self):
        '''
        routes packets from each node and appends them to DataObject Route_Table
        :return: returns True if successful
        '''
        print('Extracting nodes')
        nodes = self.graph.nodes
        route_table = _Route_table(nodes)
        print('Creating packets and routing Nodes')
        for node in tqdm(nodes):
            pop = self.graph.get_nodevalue(node,'population')
            plist = self.node_plist(node,pop)
            ppath = self.node_path(node)
            rplist = self.route_plist(plist,ppath)
            route_table.append_rplist([stop[0] for stop in ppath],rplist)
        print(f'Route_Table created. {len(route_table.content)} Packets appended')
        return route_table

    def to_file(self, fpath):
        with open(fpath,'wb') as fobj:
            pickle.dump(self,fobj)
        print('PRouter Object saved')
        return True

    @staticmethod
    def from_file(fpath):
        with open(fpath,'rb') as fobj:
            rp = pickle.load(fobj)
        print('PRouter Object loaded')
        return rp

class _Route_table:
    def __init__(self, nodes):
        self.nodes = nodes
        self.columns = ['packets', *self.nodes]
        self.size = len(self.columns)
        self.content = []

    def append_rplist(self,path_nodes,rplist):
        columns = self._find_columns(path_nodes)
        try:
            for routed_packet in rplist:
                self._append_rpacket(routed_packet, columns=columns)
        except:
            print('Could not append routed packetlist')
            raise BaseException

    def _append_rpacket(self, packet, columns=False):
        packet, stops = packet[0],packet[1]
        #checks whether columns to sort values into were already given, if not, columns are generated
        if columns == False:
            path = [stop[0] for stop in stops]
            columns = self._find_columns(path)
        #array of zeros of correct size is generated
        line = [0] * self.size
        #add packet to first field
        line[0] = packet
        #fill columns that are being passed
        try:
            for i, stop in enumerate(stops):
                line[columns[i]] = stop[1]
        except:
            raise AssertionError
        self.content.append(line)
        return True

    def _find_columns(self, path):
        '''
        :param path: list, keys to search for in self.columns
        :return: list, list of indexes of keys given in path
        '''
        try:
            return [self.columns.index(stop) for stop in path]
        except ValueError:
            print('Stop not in columns')
            print(f'Path: {path}')
            print(f'Columns: {self.columns}')
            quit

    def _extract_node(self, location, disp):
        '''
        :param node: name of the node that is to be passed to postprocessing
        :return: list, list containing all arriving packets with information in list: (packet, ta) of
        the packets at this node
        '''
        from copy import copy
        node,link = location[0],location[1]
        index = self._find_columns((node,))[0]
        nlist = [copy(row[0]).set_arrival(row[index]) for row in self.content\
                 if row[index] != 0]
        [p.set_dispersion(disp) for p in nlist]
        return {'node':node, 'link':link, 'packets':nlist}

    def to_file(self, fpath):
        with open(fpath,'wb') as fobj:
            pickle.dump(self,fobj)
        print('RTable Object saved')
        return True

    @staticmethod
    def from_file(fpath):
        with open(fpath,'rb') as fobj:
            rp = pickle.load(fobj)
        print('RTable Object loaded')
        return rp

class Postprocessing(TDict):
    def __init__(self, extracted_node, qlookup):
        self.node = extracted_node['node']
        self.link = extracted_node['link']
        self.packets = extracted_node['packets']
        self.qlut = qlookup
        start, end = qlookup.timestamps[0], qlookup.timestamps[0]+datetime.timedelta(days=1)
        timestamps = [start+i*Discretization.TIMESTEPLENGTH for i in range(int((end-start)/Discretization.TIMESTEPLENGTH))]
        super().__init__(timestamps)

    @classmethod
    def from_prouter(cls, prouter, node):
        location = [node,prouter.graph.get_outletlinks(node)[0]]
        rt_slice = prouter.route_table._extract_node(location,prouter.env.dispersion)
        return cls(rt_slice, prouter.qlookup)

    @classmethod
    def from_rtable(cls, rtable, node, qlookup, graph, env=Environment()):
        location = [node,graph.get_outletlinks(node)[0]]
        rt_slice = rtable._extract_node(location,env.dispersion)
        return cls(rt_slice, qlookup)

    def _create_tseries(self, constituent, entry_loc, load=False):
        if load is False:
            print('processing packets\n')
            entries = [p.get_loads(constituent, self.link, self.qlut, self) for p in tqdm(self.packets) if p.__contains__(constituent)]
            with open(entry_loc, 'wb') as fobj:
                pickle.dump(entries, fobj)
                print('entries saved')
        else:
            with open(entry_loc, 'rb') as fobj:
                entries = pickle.load(fobj)
                print('entries loaded')
        return TSeries(self.timestamps,entries)

    def process_constituent(self, constituent, entry_loc , load=False):
        ts = self._create_tseries(constituent, entry_loc, load)
        self.__setattr__(constituent,ts)
        return ts

    def as_conc(self, values):
        qseries = self.qlut._explode_eid(self.link)
        try:
            c = [w/(q*Discretization.TIMESTEPLENGTH.seconds) for q,w in list(zip(qseries,values))]
        except:
            pass
        return c

    def sample_times(self, duration=120, frequency="H"):
        n = np.floor(duration / 10)
        starts = pd.date_range(self.timestamps[0], self.timestamps[-1], freq=frequency)
        slots = []
        for t in starts:
            slots += pd.date_range(t, periods=n, freq="10S")
        return pd.DatetimeIndex(slots)

def test_routing(qlut, graph, save=False):
    print("Test Routing PRouter")
    router = PRouter(graph=graph, qlookup=qlut)
    routetable = router.route()
    if save:
        #routetable.to_file('C:/Users/alber/Documents/swmm/swmmpulse/route_tables/routetable.pickle')
        routetable.to_file('/mnt/c/Users/albert/Documents/SWMMpulse/route_tables/routetable.pickle')
    print('Finished Test Routing PRouter')

def test_postprocessing(qlut, graph, load=True):
    print('Test Postprocessing')
    evalnode = 'MH327-088-1'
    #routetable = _Route_table.from_file('C:/Users/alber/Documents/swmm/swmmpulse/route_tables/routetable.pickle')
    routetable = _Route_table.from_file('/mnt/c/Users/albert/Documents/SWMMpulse/route_tables/routetable.pickle')
    pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)

    #pproc.process_constituent(Loading.FECAL, entry_loc='C:/Users/alber/Documents/swmm/swmmpulse/entries/entries.pickle', load=load)
    pproc.process_constituent(Loading.FECAL, entry_loc='/mnt/c/Users/albert/Documents/SWMMpulse/entries/entries.pickle', load=load)
    print("Test postprocessing finished")
    return pproc

def preparation():
    lpath = '/mnt/c/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.out'
    gpath = '/mnt/c/Users/albert/documents/SWMMpulse/HS_calib_120_simp/'
    #gpath = 'C:/Users/alber/documents/swmm/swmmpulse/HS_calib_120_simp/'
    #lpath = 'C:/Users/alber/Documents/swmm/swmmpulse/HS_calib_120_simp.out'
    qlut = QSeries(lpath)
    graph = ntwk.from_directory(gpath)
    return qlut, graph

def main():
    qlut, graph = preparation()
    test_routing(qlut,graph,True)
    test_postprocessing(qlut, graph, load=True)

if __name__ == '__main__':
    main()