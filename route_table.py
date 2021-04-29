import numpy as np

class Route_table:
    def __init__(self, nodes):
        self.nodes = nodes
        self.columns = ['packets', *self.nodes]
        self.size = len(columns)
        self.content = []

    def append_routed_node(self, routed_node):
        '''
        takes routed_node dictionary from DataType Route_table and appends it to itself
        :param routed_node: dictionary, API for Route_table
        :return: False if run unsuccessful
        '''
        path = routed_node['path']
        packets = routed_node['packets']
        columns = self.find_columns(path)
        try:
            for packet in packets:
                self.append_routed_packet(path,packet,columns=columns)
        except Exception:
            return False

    def append_routed_packet(self,path, packet, columns=False):
        try:
            #checks whether columns to sort values into were already given, if not, columns are generated
            if columns == False:
                cols = self.find_columns(path)
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

    def find_columns(self, path):
        '''
        :param path: list, keys to search for in self.columns
        :return: list, list of indexes of keys given in path
        '''
        try:
            return [self.columns.index(stop) for stop in path]
        except Exception:
            return False

    def pass_routed_node(self, node):
        '''
        :param node: name of the node that is to be passed to postprocessing
        :return: list, list containing all arriving packets with information in list: (packet, ta) of
        the packets at this node
        '''
        node_column = self.find_columns((node,))[0]
        arrivals = [(row[0],row[node_columns]) for row in self.content]
        return arrivals

