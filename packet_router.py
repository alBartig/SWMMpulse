import geopandas
import pandas as pd
import network_lib as ntwk
import timeseries
import create_loadings
import numpy as np
from swmm_api import read_out_file
import datetime
from route_table import Route_table
from pconstants import Discretization

class Packet_router:
    def __init__(self,graph_location,lookup_location):
        self.lookup = timeseries.qlookup(lookup_location)
        self.graph = self.load_graph(graph_location)

    def load_graph(self, fpath):
        links = []
        linkfiles = ['Conduits','Pumps','Orifices']
        for file in linkfiles:
            gdf = geopandas.read_file(fpath+linkfiles[0] +'.shp')
            [links.append(list(link)) for link in gdf[['NAME','INLETNODE','OUTLETNODE','LENGTH']].to_numpy()]
        for i in range(len(links)):
            links[i][3]=('length',links[i][3])
        graph = ntwk.Network(links)
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

    def init_route_table(self):
        self.route_table = Route_table(self.graph.get_nodes())

    def route_packet(self, packet, path):
        routed_packet = [(path[0][0],self.t0)]
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

    def create_community(self,node):
        com = create_loadings.community(node)
        pop = self.graph.get_nodevalue(node,'population')
        for group in self.community_info:
            com['population'] = pop * group['weight']
            com.group_from_dict(group)

    def generate_seeds(self, node):
        com = self.create_community(node)
        return com.community_timeseries()

    def simulate_node(self, node):
        '''
        function takes node, generates seedlist via stored community info and returns dictionary with a list of stops
        and the a list with the routed packets
        :param node: str, node that is to be simulated
        :return: dict, {path:list of nodes, packets:list of tuples (packet, list with arrival times)}
        '''
        routed_list = []
        path = self.graph.trace_path(node)
        seedlist = self.generate_seeds(node)
        for packet in seedlist:
            routed_packet = self.route_packet(packet,path=path)
            routed_list.append(routed_packet)
        return {path:[stop[0] for stop in path],'packets':routed_list}

    def simulate(self):
        '''
        routes packets from each node and appends them to DataObject Route_Table
        :return: returns True if successful
        '''
        try:
            for node in self.graph.get_nodes():
                routed_node = self.simulate_node(node)
                self.route_table.append_routed_node(routed_node)
            return True
        except Exception:
            return False

if __name__ == '__main__':
    gpath = 'C:/Users/alber/Documents/SWMM/Training6_GA/Training6_GA/'
    lpath = 'C:/Users/alber/Documents/SWMM/Training6_GA/'
    sim = Packet_router(graph_location=gpath,lookup_location=lpath)
    date = datetime.time(hour=0,minute=15,second=0)