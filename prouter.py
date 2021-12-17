#import geopandas
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
import os
from copy import copy, deepcopy

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

    def node_plist(self, nname, npop, startid):
        """
        Creates packetlist (plist) for node
        Args:
            node (dict): {'name':nodename,'population':overall population at node}

        Returns:
            list: list with pobjects for each group in self.env
        """
        plist = []
        for group in self.env.groups:
            new_plist = group.plist(nname,npop,self.env.date,startid)
            startid += len(new_plist)
            plist += new_plist
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

    def route(self, rtname):
        '''
        routes packets from each node and appends them to DataObject Route_Table
        :return: returns True if successful
        '''
        #print('Extracting nodes')
        nodes = self.graph.nodes
        route_table = _Route_table(nodes, rtname)
        startid = 0
        #print('Creating packets and routing Nodes')
        for node in nodes:
        #for node in tqdm(nodes, desc="Routing packets..."):
            try:
                pop = int(self.graph.get_nodevalue(node,'POP'))
            except:
                pop = 0
            plist = self.node_plist(node,pop,startid)
            startid += len(plist)
            ppath = self.node_path(node)
            rplist = self.route_plist(plist,ppath)
            route_table.append_rplist([stop[0] for stop in ppath],rplist)
        #print(f'Route_Table created. {len(route_table.content)} Packets appended')
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
    def __init__(self, nodes, name):
        self.nodes = nodes
        self.columns = ['packets', *self.nodes]
        self.size = len(self.columns)
        self.content = []
        self.name = name

    def append_rplist(self,path_nodes,rplist):
        columns = self._find_columns(path_nodes)
        for routed_packet in rplist:
            self._append_rpacket(routed_packet, columns=columns)

    def _append_rpacket(self, packet, columns=False):
        packet, stops = packet[0],packet[1]
        #checks whether columns to sort values into were already given, if not, columns are generated
        if columns == False:
            path = [stop[0] for stop in stops]
            columns = self._find_columns(path)
        #array of nan of correct size is generated
        line = self.size*[np.NaN]
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
        node,link = location[0],location[1]
        index = self._find_columns((node,))[0]
        nlist = [copy(row[0]).set_arrival(row[index]) for row in self.content if type(row[index]) != float]
        [p.set_dispersion(disp) for p in nlist]
        return {'node':node, 'link':link, 'packets':nlist}

    def to_file(self, fpath):
        with open(fpath,'wb') as fobj:
            pickle.dump(self,fobj)
        print('RTable Object saved')
        return True

    def to_parquet(self, fpath, compression="brotli"):
        #import time
        #extract packets and write to table
        packets = [r[0].to_list() for r in self.content]
        #packets = [r[0].to_list() for r in tqdm(self.content, desc="Saving packets...")]
        dfp = pd.DataFrame(packets, columns=self.content[0][0].tags)
        dfp.to_parquet(os.path.join(fpath,self.name+"_packets.parquet"), compression=compression)

        #write routetable to table
        #start = time.time()
        _content = []
        for row in self.content:
            _content.append([row[0],*row[1:]])
        #print(f"Copied rt-contents: {time.time()-start:.3f} s")
        for row in _content:
        #for row in tqdm(_content, desc="Saving routetable..."):
            row[0] = row[0].pid
        dfrt = pd.DataFrame(_content, columns=self.columns).set_index("packets").astype("str")
        dfrt.to_parquet(os.path.join(fpath,self.name+"_routetable.parquet"), compression="brotli")
        #print('RTable Object saved')
        return True

    def from_parquet(self, fpath):
        dfrt = pd.read_parquet(fpath)
        #dfp =

    @staticmethod
    def from_file(fpath):
        with open(fpath,'rb') as fobj:
            rp = pickle.load(fobj)
        print('RTable Object loaded')
        return rp

    def postprocess(self, location, disp):
        location = [node,graph.get_outletlinks(node)[0]]
        rt_slice = rtable._extract_node(location,env.dispersion)

class Postprocessing:
    def __init__(self, extracted_node, qlookup, name):
        self.node = extracted_node['node']
        self.link = extracted_node['link']
        self.packets = extracted_node['packets']
        self.qlut = qlookup
        self.name = name

    @classmethod
    def from_prouter(cls, prouter, node):
        location = [node,prouter.graph.get_outletlinks(node)[0]]
        rt_slice = prouter.route_table._extract_node(location,prouter.env.dispersion)
        return cls(rt_slice, prouter.qlookup)

    @classmethod
    def from_rtable(cls, rtable, node, qlookup, graph, env=Environment()):
        location = [node,graph.get_outletlinks(node)[0]]
        rt_slice = rtable._extract_node(location,env.dispersion)
        return cls(rt_slice, qlookup, "pproc"+rtable.name.strip("rt"))

    @classmethod
    def from_file(cls, fpath):
        with open(fpath,'rb') as fobj:
            rp = pickle.load(fobj)
        print(f'PProc Object loaded from: {fpath}')
        return rp

    def _create_tseries(self, constituent, entry_loc, load=False, name=None):
        if load is False:
            #print('processing packets\n')
            entries = []
            for packet in self.packets:
            #for packet in tqdm(self.packets, f"Calculating loads from packets for {constituent}..."):
                if packet.__contains__(constituent):
                    entries.append(packet.get_entry(constituent, self.link, self.qlut))
            # with open(entry_loc, 'wb') as fobj:
            #     pickle.dump(entries, fobj)
            #     print('entries saved')
        else:
            with open(entry_loc, 'rb') as fobj:
                entries = pickle.load(fobj)
                print('entries loaded')
        tagnames = list(set(entries[0].keys()).difference(["pid","values","timestamps"]))
        return TSeries(self.qlut.timestamps, entries, tagnames=tagnames, name="_".join([name,constituent]))

    def process_constituent(self, constituent, entry_loc , load=False):
        ts = self._create_tseries(constituent, entry_loc, load, "ts" + self.name.strip("pproc"))
        self.__setattr__(constituent, ts)
        return ts

    def as_conc(self, values):
        qseries = self.qlut._explode_eid(self.link)
        try:
            c = [w/(q*Discretization.TIMESTEPLENGTH.seconds) for q,w in list(zip(qseries,values))]
        except:
            pass
        return c

    def save(self, opath):
        with open(opath,'wb') as fobj:
            pickle.dump(self,fobj)
        print(f'PProc Object saved: {opath}')
        return True

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