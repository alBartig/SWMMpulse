import geopandas
from network_lib import DirectedTree as ntwk
from timeseries import qlookup
import create_loadings
import numpy as np
import datetime
from environment import Environment
from post_processer import Postprocessing

class Packet_router:
    def __init__(self,graph_location,lookup_location):
        self.lookup = qlookup(lookup_location)
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
        stops = [(packet.origin,packet.t0)]
        ct = packet.t0
        for stop in path[1:]:
            node = stop[0]
            link = stop[1]
            v = self.lookup.lookup_v(link,ct)
            td = self.graph.get_linkvalue(link,'LENGTH') / v
            ct = ct + datetime.timedelta(seconds=td)
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
        return plist

    def node_path(self, node):
        '''
        function takes node, generates seedlist via stored community info and returns dictionary with a list of stops
        and the a list with the routed packets
        :param node: str, node that is to be simulated
        :return: dict, {path:list of nodes, packets:list of tuples (packet, list with arrival times)}
        '''
        return self.graph.trace_path(node)

    def route(self):
        '''
        routes packets from each node and appends them to DataObject Route_Table
        :return: returns True if successful
        '''
        try:
            nodes = self.graph.nodes
            route_table = Route_table(nodes)
            for node in nodes:
                pop = self.graph.get_nodevalue(node,'population')
                plist = self.node_plist(node,pop)
                ppath = self.node_path(node)
                rplist = self.route_plist(plist,ppath)
                route_table.append_rplist(ppath)
            return self.route_table
        except Exception:
            return False

class Route_table:
    def __init__(self, nodes):
        self.nodes = nodes
        self.columns = ['packets', *self.nodes]
        self.size = len(self.columns)
        self.content = []

    def append_rplist(self,ppath,rplist):
        columns = self._find_columns(ppath)
        for routed_packet in rplist:
            self._append_rpacket(ppath, routed_packet, columns=columns)

    def _append_rpacket(self, path, packet, columns=False):
        try:
            #checks whether columns to sort values into were already given, if not, columns are generated
            if columns == False:
                cols = self._find_columns(path)
            else:
                cols = columns
            arrivaltimes = packet[1]
            packet = packet[0]
            #array of zeros of correct size is generated
            line = np.zeros(self.size)
            #add packet to first field
            line[0] = packet
            #fill columns that are being passed
            for i, value in enumerate(arrivaltimes):
                line[cols[i]] = value
            self.content.append(line)
            return True
        except Exception:
            return False

    def _find_columns(self, path):
        '''
        :param path: list, keys to search for in self.columns
        :return: list, list of indexes of keys given in path
        '''
        try:
            return [self.columns.index(stop) for stop in path]
        except Exception:
            return False

    def extract_node(self, node):
        '''
        :param node: name of the node that is to be passed to postprocessing
        :return: list, list containing all arriving packets with information in list: (packet, ta) of
        the packets at this node
        '''
        index = self._find_columns((node,))[0]
        nlist = [(row[0],row[index]) for row in self.content]
        return {'node':node, 'packets':nlist}

if __name__ == '__main__':
    gpath = '/mnt/c/Users/albert/documents/SWMMpulse/HS_calib_120_simp/'
    lpath = '/mnt/c/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.out'
    evalnode = 'MH327-088-1'
    sim = Packet_router(graph_location=gpath,lookup_location=lpath)
    rtable = sim.route()
    pproc = Postprocessing(rtable.extract_node(evalnode),sim.env)