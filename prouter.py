import datetime
from environment import Environment, PACKET, CONSTITUENT, DirectedTree, HYDRAULICS, GROUP, DEFAULT
from Exceptions import RoutingError, PlausibilityError
import pandas as pd
import numpy as np

class Router:
    def __init__(self):
        pass

    def add_environment(self, env):
        """
        Adds environment to the router
        Args:
            env (Environment): environment

        Returns:
            None
        """
        self.environment = env
        return None

    def _route_packet(self, packet):
        """
        routs each packet by tracing the path in the system and then adding up
        traveltimes between each node
        Args:
            packet (dict): dictionary of the packet to be routed

        Returns:
            list
        """
        path = self.graph.trace_path(packet.get(PACKET.ORIGIN)) #trace path
        ct = packet.get(PACKET.T0) #get origin time
        stops = [(packet.get(PACKET.ORIGIN), ct)] #prepare list for stops
        for stop in path[1:]:
            node, link = stop[0], stop[1]
            #lookup velocity, if v == 0 break routing
            v = self.environment.get_velocity(ct, link)
            #v,delay = self.flows.lookup_v(link, ct) #lookup flowvelocity (and delay until v > 0)
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
            ct = ct + datetime.timedelta(seconds=td)# + delay
            stops.append((node, ct)) #append stop
            #check for plausibility:
            t0 = packet.get("t0")
            for stop in stops: #this could maybe be eliminated!!
                if stop[0] < t0:
                    raise RoutingError(packet, stops)
        return stops

    def route(self, packets=None, as_dataframe=False):
        """
        Create routetable.
        Prepare input data ahead with: add_flow(), add_graph(), add_environment()
        Returns:
            dict
        """
        if packets is None:
            #generate packets from environment
            packets = self.environment.get_packets()
        #order graph
        self.environment.graph.order_shreve()
        nodes = [(name, values["shreve"], values["outlets"]) for name, values in self.environment.graph.adjls.items()]
        nodes.sort(key=lambda x: x[1]) #sort list of nodes by shreve order
        columns = {node[0]:{} for node in nodes} #prepare container for routed packets
        for node in nodes[:-1]: #iterate through columns until last column
            #load packets that originate from node into column
            node, onode, olink = node[0], node[2][0][0], node[2][0][1] #node=current node, onode=next node, olink=next link
            columns[node].update(packets.loc[packets[PACKET.ORIGIN] == node].
                                 set_index(PACKET.PACKETID)[PACKET.T0].to_dict()) #load packets originating from current node
            distance = self.environment.graph.get_linkvalue(olink, HYDRAULICS.LENGTH) #get distance for link
            starttimes = list(columns[node].values()) #get starttimes from current node
            packetids = list(columns[node].keys()) #get packetids from current node
            #get velocities for outlet link by using starttimes as indexes. Round starttimes and take modulo 8640
            velocities = self.environment.flow_velocities.get(olink)[np.rint(starttimes).astype(int) % 8640]
            #calculate arrival times at next node by adding distance/velocity/10s
            endtimes = starttimes + distance / velocities / 10
            #endtimes = (starttimes + np.rint(distance / velocities).astype("int")) % 8640
            columns[onode].update({packetid: endtime for packetid, endtime in zip(packetids, endtimes)})

        if as_dataframe:
            routetable = pd.DataFrame.from_dict(columns) #generate DataFrame from columns-dictionary
            return routetable
        else:
            return columns

    def postprocess(self, routetable, packets, constituent, node=None, as_df=False):
        """
        Creates a timeseries-table pd.DataFrame from a given routetable for a given node and constituent.
        The DataFrame contains a  column for each packet routed with numerical values for timesteps
        that contain the constituent and np.NaN values for timesteps that don't.
        Args:
            constituent (dict): dict of the constituent for which to create the timeseries
            routetable (pd.DataFrame): routetable DataFrame, returned from self.route()
            node (str): name of the node for which to create the timeseries

        Returns:
            pd.DataFrame
        """
        if node is None:
            node = self.environment.graph.root
            link = self.environment.graph.get_inletlinks(node)[0]
        else:
            link = self.environment.graph.get_outlets(node)[0][1]

        # unpack constituent parameters ------------------------------------------------------
        fractions = constituent.get(CONSTITUENT.FRACTIONS)
        skewedness = constituent.get(CONSTITUENT.SKEWEDNESS)
        skewedness_left = [s[0] for s in skewedness]
        skewedness_right = [s[1] for s in skewedness]
        decay_rate = constituent.get(CONSTITUENT.DECAY_RATE) #decay-rate in 1/d
        specific_load = constituent.get(CONSTITUENT.SPECIFIC_LOAD)
        groups = constituent.get(GROUP.GROUPS)
        dispersion_rate = self.environment.information.get(HYDRAULICS.DISPERSION_RATE)

        # prepare dataframe for postprocessing
        df_post = packets.loc[packets["classification"].isin(groups)].set_index("pid")
        df_post["tm"] = pd.Series(routetable[self.environment.graph.root])
        arr_tm = np.rint(df_post["tm"].values).astype(int) % 8640 + 8640

        # calculate additional values for packets
        arr_age = (df_post["tm"].values - df_post["t0"].values) * 10# age of packets in seconds
        arr_load = specific_load * np.e ** (-decay_rate * arr_age / 86400)# reduced load
        arr_2sd_m = 2 * (2 * dispersion_rate * arr_age) ** 0.5
        arr_flow_velocity = self.environment.flow_velocities.get(link)[np.rint(arr_tm).astype(int)%8640]
        # two standard devitions in 10s steps
        arr_2sd_s = np.repeat([np.ceil(arr_2sd_m / arr_flow_velocity / 10)], len(fractions), axis=0).T

        # prepare numpy arrays for calculation
        arr_dxleft = arr_2sd_s * skewedness_left
        arr_dxright = arr_2sd_s * skewedness_right
        arr_borderleft = np.floor(np.repeat([arr_tm], len(fractions), axis=0) - arr_dxleft.T)
        arr_borderright = np.ceil(np.repeat([arr_tm], len(fractions), axis=0) + arr_dxright.T)

        #calculate fractions
        arr_load_frctns = np.repeat([arr_load], len(fractions), axis=0).T * fractions
        arr_peak_frctns = ((2*arr_load_frctns).T / arr_2sd_s.T).T /\
                          (np.sum(skewedness, axis=0) * arr_2sd_s)
        arr_dcleft = arr_peak_frctns / arr_dxleft
        arr_dcright = arr_peak_frctns / arr_dxright
        arr_tm_frctns = np.repeat([arr_tm], len(fractions), axis=0)

        arr_timeseries = np.zeros([len(df_post), len(fractions), 8640 * 3])# prepare
        for i, packet in enumerate(zip(arr_timeseries, arr_tm_frctns.T, arr_dxleft, arr_dxright, arr_borderleft.T, arr_borderright.T,
                                       arr_peak_frctns, arr_dcleft, arr_dcright)):
            self._process_packet(packet)

        arr_timeseries = np.sum(arr_timeseries, axis=1)
        ts1, ts2, ts3 = np.split(arr_timeseries.T, 3)
        arr_timeseries = ts1 + ts2 + ts3
        dic_timeseries = {pid: timeseries for pid, timeseries in zip(df_post.index, arr_timeseries.T)}

        if not as_df:
            return dic_timeseries
        else:
            df = pd.DataFrame(dic_timeseries, index = pd.date_range("2000-01-01", freq="10S", periods=8640))
            return df

    def _process_packet(self, packet):
        t, tm, dxleft, dxright, borderleft, borderright, peaks_fracts, dcleft, dcright = packet
        for fraction in zip(t, tm, dxleft, dxright, borderleft, borderright, peaks_fracts, dcleft, dcright):
            self._process_fraction(fraction)
        return None

    def _process_fraction(self, fraction):
        x, tm, dxleft, dxright, borderleft, borderright, peak, dcleft, dcright = fraction
        try:
            x[borderleft:tm] = -dcleft * abs(np.arange(borderleft - tm, 0)) + peak
        except:
            print(f"{borderleft} : {tm}\n",
                  f"{tm - borderleft} : 0")
            x[borderleft:tm] = -dcleft * abs(np.arange(borderleft - tm, 0)) + peak
        try:
            x[tm:borderright] = -dcright * abs(np.arange(0, borderright - tm)) + peak
        except:
            print(f"{tm} : {borderright}\n",
                  f"{borderright - tm} : 0")
            x[tm:borderright] = -dcright * abs(np.arange(0, borderright - tm)) + peak
        return x


def preparation_env():
    print(f"Preparing environment")
    env = Environment()
    print(f"Reading swmm-outfile")
    env.read_swmmoutfile(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp.out")
    print(f"Preparing graph")
    graph = DirectedTree.from_swmm(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp.inp")
    node_data = pd.read_csv(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp/pop_node_data.csv")
    node_data = node_data.set_index("NAME").to_dict(orient="index")
    graph.add_nodevalues(node_data)
    env.add_graph(graph)
    print(f"finished preparing environment")
    return env

def testing():
    import time
    env = preparation_env()
    print(f"Preparing router")
    router = Router()
    router.add_environment(env)
    print(f"finished preparing router")
    print(f"beginning router testing")
    start = time.time()
    packets = router.environment.get_packets()
    routetable = router.route(packets=packets)
    print(f"time for routing: {time.time()-start} seconds")
    print(f"testing postprocessing")
    processed = router.postprocess(routetable, packets, DEFAULT.DEFAULT_CONSTITUENTS.get(CONSTITUENT.COV))
    print(f"finished router testing")

def main():
    testing()

if __name__ == '__main__':
    main()