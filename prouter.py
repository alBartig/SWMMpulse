import geopandas
from network_lib import DirectedTree as ntwk, GraphConstants as gc
from timeseries import QSeries, TSeries, TDict
import create_loadings
from pconstants import Discretization, Loading
import datetime
from environment import Environment
from tqdm import tqdm
import pickle

class PRouter:
    def __init__(self,graph_location,qlookup):
        self.qlookup = qlookup
        self.graph = self.load_graph(graph_location)
        self.env = Environment()

    def load_graph(self, fpath):
        links = []
        linkfiles = ['Conduits','Pumps','Orifices']
        nodefiles = ['Junctions']
        for file in linkfiles[:1]:
            gdf = geopandas.read_file(fpath+file +'.shp')
            [links.append(list(link)) for link in gdf[['NAME','INLETNODE','OUTLETNODE','LENGTH']].to_numpy()]
        for i in range(len(links)):
            links[i][3]=('length',links[i][3])
        graph = ntwk(links)
        nodes = []
        for file in nodefiles:
            gdf = geopandas.read_file(fpath+file +'.shp')
            [nodes.append(list(node)) for node in gdf[['NAME', 'POP']].to_numpy()]
        for n,p in nodes:
            graph.add_nodevalues(n,[('name',n),('population',p)])
        return graph

    def set_community_info(self, groups:iter, weights:iter, **kwargs):
        '''
        Sets standard information on the modeled communities
        :param groups: classification of mutually exclusive groups, eg: ('healthy','infected')
        :param weights: weights by which a given population is divided into the groups, eg: (0,9,0.1)
        :param kwargs: keyword arguments that are optional:
            distribution:iter, tuple of applicable pattern distribution, if none are given, 'Bristol' is applied
            ppp:iter, poops per person per day, number of defecation events per person per day, default is 1
        :return: stores community information as a dictionary within pattern_router object
        '''
        com_info = {}
        weights_adjusted = [weight/sum(weights) for weight in weights]
        for i in range(len(groups)):
            com_info[groups[i]] = {'weight':weights_adjusted[i],'distribution':kwargs.get('distribution')[i],\
                                   'ppp':kwargs.get('poops_per_person')[i]}
        self.community_info = com_info

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
            ct = ct + datetime.timedelta(seconds=td) + delay
            stops.append((node,ct))
        return (packet,stops)

    def route_plist(self,plist, path):
        return [self.route_packet(p, path) for p in plist]

    def create_community(self,node):
        com = create_loadings.community(node)
        pop = self.graph.get_nodevalue(node,'population')
        for group in self.community_info:
            com['population'] = pop * group['weight']
            com.group_from_dict(group)

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
        self.route_table = route_table
        print(f'Route_Table created. {len(route_table.content)} Packets appended')
        return True

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
        nlist = [copy(row[0]).set_arrival(row[index], disp) for row in self.content\
                 if row[index] != 0]
        return {'node':node, 'link':link, 'packets':nlist}

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

    def _create_tseries(self, constituent):
        entries = [p.get_loads(constituent, self.link, self.qlut) for p in self.packets]
        return TSeries(self.timestamps,entries)

    def process_constituent(self,constituent):
        ts = self._create_tseries(constituent)
        self.__setattr__(constituent,ts)
        return ts

if __name__ == '__main__':
    lpath = '/mnt/c/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.out'
    qlut = QSeries(lpath)
    skip = True
    if skip is True:
        print('Skipping Test PRouter')
    else:
        #gpath = 'C:/Users/alber/documents/swmm/swmmpulse/HS_calib_120_simp/'
        #lpath = 'C:/Users/alber/Documents/swmm/swmmpulse/HS_calib_120_simp.out'
        gpath = '/mnt/c/Users/albert/documents/SWMMpulse/HS_calib_120_simp/'
        sim = PRouter(graph_location=gpath, qlookup=qlut)
        sim.route()
        sim.to_file('/mnt/c/Users/albert/Documents/SWMMpulse/prouter.pickle')
    #print('Finished Test PRouter')
    print('Test Postprocessing')
    evalnode = 'MH327-088-1'
    sim = PRouter.from_file('/mnt/c/Users/albert/Documents/SWMMpulse/prouter.pickle')
    pproc = Postprocessing.from_prouter(sim, evalnode)
    pproc.process_constituent(Loading.FECAL)
    print('Test Postprocessing finished')