from timeseries import QSeries, QFrame
from pobject import PObject
from prouter import PRouter, Postprocessing, _Route_table
from network_lib import DirectedTree as ntwk, GraphConstants as gc
from environment import Environment
import os
from pconstants import Loading
import matplotlib.pyplot as plt
from multiprocessing import Pool
from tqdm import tqdm
from helpers import series_to_dfmi
import logging

def prepare_environment(lpath='C:/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.out',\
                        gpath = 'C:/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.inp',
                        pop_data = 'C:/Users/albert/Documents/SWMMpulse/HS_calib_120_simp/Junctions.dbf',
                        logger=logging.getLogger("root")):
    #lpath = '/mnt/c/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.out'
    #lpath = 'C:/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.out'
    logger.info("...loading hydraulic lookup")
    qlut = QFrame(lpath)
    #gpath = '/mnt/c/Users/albert/documents/SWMMpulse/HS_calib_120_simp.inp'
    #pop_data = '/mnt/c/Users/albert/documents/SWMMpulse/HS_calib_120_simp/Junctions.dbf'
    #gpath = 'C:/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.inp'
    #pop_data = 'C:/Users/albert/Documents/SWMMpulse/HS_calib_120_simp/Junctions.dbf'
    logger.info("...loading graph model")
    graph = ntwk.from_swmm(gpath)
    logger.info("...loading population data")
    graph.nodedata_from_dbf(pop_data, cols="POP")
    return graph,qlut

def generate_rts(graph, qlut):
    fractions = [0.00070,0.00050,0.00030,0.00010]
    for fraction in fractions:
        env = Environment()
        env.groups[0].weight = (1 - fraction)
        env.groups[1].weight = fraction

        for i in range(25):
            fname = f"rt_f{f'{fraction:.5f}'[-5:]}_{str(i).zfill(2)}.pickle"
            router = PRouter(graph=graph, qlookup=qlut)
            router.env = env
            routetable = router.route()
            #routetable.to_file(f'C:/Users/alber/documents/swmm/swmmpulse/route_tables/{fname}')
            routetable.to_file(f'/mnt/c/Users/albert/Documents/SWMMpulse/route_tables/{fname}')
            #pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)

def evaluate_rts(graph, qlut, evalnode):
    #path = 'C:/Users/alber/documents/swmm/swmmpulse/'
    path = '/mnt/c/Users/albert/documents/SWMMpulse/'
    files = [p for p in os.listdir(os.path.join(path, "route_tables")) if p[:3] == 'rt_']
    for file in files[1:]:
        routetable = _Route_table.from_file(os.path.join(path,"route_tables",file))
        pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)
        eval_constituents = [Loading.FECAL, Loading.COV]
        for const in eval_constituents:
            entry_loc = os.path.join(path, 'entries', f"{file.rstrip('.pickle')}_{const}.pickle")
            pproc.process_constituent(const, entry_loc, load=False)
            pproc.Fecal_Matter.timeseries(wpath=os.path.join(path,'timeseries',f"{file.rstrip('.pickle')}_{const}.csv"))

            print('postprocesser loaded')
    print('finish')

def create_cseries(graph, qlut, evalnode):
    #path = 'C:/Users/alber/documents/swmm/swmmpulse/'
    path = '/mnt/c/Users/albert/documents/SWMMpulse/'
    files = [p for p in os.listdir(os.path.join(path, "route_tables")) if p[:3] == 'rt_']
    ls = []
    sid = 0
    for file in files[:]:
        routetable = _Route_table.from_file(os.path.join(path,"route_tables",file))
        ratio = round(len([p for p in routetable.content if p[0].classification != "Healthy"])/len(routetable.content),6)
        sid += 1
        temp = {"node":evalnode,"ratio":ratio,"id":sid}
        pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)
        eval_constituents = [Loading.FECAL, Loading.COV]
        for const in eval_constituents:
            entry_loc = os.path.join(path, 'entries', f"{file.rstrip('.pickle')}_{const}.pickle")
            pproc.process_constituent(const, entry_loc, load=True)
            t,m = getattr(pproc,const).timeseries()
            c = pproc.as_conc(m)
            temp[const] = c
        ls.append(temp)
    return ls, t

def load_pproc(rtable):
    #path = 'C:/Users/alber/documents/swmm/swmmpulse/'
    path = '/mnt/c/Users/albert/documents/SWMMpulse/'
    pfile = Postprocessing.from_rtable(rtable)

def plot_data(graph, qlut, evalnode):
    path = 'C:/Users/alber/documents/swmm/swmmpulse'
    file = 'rt_f.0001_00.pickle'
    routetable = _Route_table.from_file(os.path.join(path,file))
    pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)
    pproc.process_constituent(Loading.COV, entry_loc=os.path.join(path,'entries/rt_f.0001_00_Cov_RNA.pickle'), load=True)
    pproc.process_constituent(Loading.FECAL, entry_loc=os.path.join(path,'entries/rt_f.0001_00_Fecal_Matter.pickle'), load=True)
    t,c = pproc.Fecal_Matter.timeseries()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    plt.plot(t,c)

    print('Series loaded')

    print('shit')

def testdrive(graph, qlut, evalnode):
    fpath = f'/mnt/c/Users/albert/documents/SWMMpulse/'
    #fpath = f'C:/Users/alber/documents/swmm/swmmpulse/'
    router = PRouter(graph=graph, qlookup=qlut)
    router.env = Environment()
    routetable = router.route()
    routetable.to_file(os.path.join(fpath,'route_tables','rtb_testdrive.pickle'))
    pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)
    eval_constituent = Loading.FECAL
    entry_loc = os.path.join(fpath, 'entries', "entries_testdrive.pickle")
    pproc.process_constituent(eval_constituent, entry_loc, load=False)
    pproc.Fecal_Matter.timeseries(wpath=os.path.join(fpath, 'timeseries', 'timeseries_testdrive.pickle'))
    print("Function finished")

def test_cov(graph, qlut, evalnode):
    path = 'C:/Users/alber/documents/swmm/swmmpulse'
    file = 'route_tables/rt_f00003_00.pickle'
    routetable = _Route_table.from_file(os.path.join(path,file))
    pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)
    pproc.process_constituent(Loading.COV, entry_loc=os.path.join(path,'entries/rt_f00003_00_Cov_RNA.pickle'), load=False)
    print('test_cov finished')

def simulation(sim_dict, debug=False):
    i = sim_dict["i"]
    graph = sim_dict["graph"]
    qlut =sim_dict["qlut"]
    env = Environment()
    path = sim_dict["path"]
    evalnode = sim_dict["evalnode"]
    constituents = sim_dict["constituents"]
    prefix = sim_dict["prefix"]
    logger = sim_dict["logger"]

    #print(f"VARIATION {i}\n")
    fraction = 0.0001 + 0.0006 * i / 1000
    env.groups[0].weight = (1 - fraction)
    env.groups[1].weight = fraction
    fname = f"{prefix}_variation_{str(i).zfill(4)}"

    # Create routetable
    logger.debug("...creating Packetrouter")

    router = PRouter(graph=graph, qlookup=qlut)
    router.env = env

    logger.debug("...creating Routetable")

    routetable = router.route(fname)
    routetable.to_parquet(os.path.join(path, "route_tables"))

    logger.debug("...creating Postprocessing from Routetable")

    # Slice routetable
    pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)

    logger.debug("...calculate Timeseries for constituents")

    # Create timeseries for each constituent
    for constituent in constituents:
        pproc.process_constituent(constituent, os.path.join(path, "entries"))
        pproc.__getattribute__(constituent).to_parquet(os.path.join(path, "entries"))

    logger.debug("check timeseries!")
    logger.info("simulation finished")

def mp_simulation():
    logging.basicConfig(level=50, force=True)
    logger = logging.getLogger("Testsimulation")
    # Settings
    #path = f'/mnt/c/Users/albert/documents/SWMMpulse/'
    path = f'C:/Users/albert/documents/SWMMpulse/'
    n = 1000
    evalnode = 'MH327-088-1'
    env = Environment()
    graph, qlut = prepare_environment()
    constituents = [Loading.COV]
    prefix = "s3"

    sims = [{"i":i, "graph":graph, "qlut":qlut, "constituents":constituents, "evalnode":evalnode, "path":path,\
             "env":env, "prefix":prefix, "logger":logger} for i in tqdm(range(n), desc="Preparing inputs for variations...")]
    pool = Pool(processes=6)
    for _ in tqdm(pool.imap(simulation, sims, chunksize=6), desc="Processing variations...", total=len(sims)):
        pass
    pool.close()
    pool.join()
    print("Multiprocessing finished")

def test_sim():
    logging.basicConfig(level=10)
    logger = logging.getLogger("Testsimulation")

    path = f'C:/Users/albert/documents/SWMMpulse/'
    evalnode = 'MH327-088-1'
    prefix = "testsim"
    logger.info("...building environment")
    env = Environment()
    graph, qlut = prepare_environment()
    constituents = [Loading.COV, Loading.FECAL]
    logger.info("...simulate")
    simulation({"i":1000, "graph":graph, "qlut":qlut, "constituents":constituents, "evalnode":evalnode, "path":path,\
             "env":env, "prefix":prefix, "logger":logger}, debug=True)
    logger.info("process ended")

def post_calculate_fecal_ts():
    path = f'C:/Users/albert/documents/SWMMpulse/'
    evalnode = 'MH327-088-1'
    env = Environment()
    graph, qlut = prepare_environment()
    constituents = [Loading.FECAL]

    def calc_ts(sim_dict):
        rtfile = sim_dict["rtfile"]
        graph = sim_dict["graph"]
        qlut =sim_dict["qlut"]
        path = sim_dict["path"]
        evalnode = sim_dict["evalnode"]
        constituents = sim_dict["constituents"]
        # Slice routetable
        routetable = _Route_table.from_file(rtfile)
        pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)

        # Create timeseries for each constituent
        for constituent in constituents:
            pproc.process_constituent(constituent, os.path.join(path, "entries"))
            pproc.__getattribute__(constituent).to_parquet(os.path.join(path, "entries"))

    rtfiles = [file for file in os.listdir(os.path.join(path, "route_tables")) if len(file)==36]
    sims = [{"rtfile": rtfile, "graph": graph, "qlut": qlut, "constituents": constituents, "evalnode": evalnode, "path": path, \
             "env": env} for rtfile in tqdm(rtfiles, desc="Loading rtfiles...")]

    pool = Pool(processes=6)
    for _ in tqdm(pool.imap(calc_ts, sims, chunksize=6), desc="Processing variations...", total=len(sims)):
        pass

def bulk_simulation():
    # Settings
    #path = f'/mnt/c/Users/albert/documents/SWMMpulse/'
    path = f'C:/Users/albert/documents/SWMMpulse/'
    n = 10
    evalnode = 'MH327-088-1'
    env = Environment()
    graph, qlut = prepare_environment()
    constituents = [Loading.FECAL, Loading.COV]

    for i in range(n):
        print(f"SCENARIO {i}\n")
        fraction = 0.1
        env.groups[0].weight = (1 - fraction)
        env.groups[1].weight = fraction
        fname = f"rt_scenario_{str(i).zfill(3)}"

        # Create routetable
        print("...initializing route table")
        router = PRouter(graph=graph, qlookup=qlut)
        router.env = env
        print("...creating and routing packets")
        routetable = router.route(fname)
        print("...saving route table")
        routetable.to_parquet(os.path.join(path, "route_tables"))

        # Slice routetable
        print("...creating postprocessing object")
        pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)

        # Create timeseries for each constituent
        for constituent in constituents:
            print(f"...calculating timeseries for {constituent}")
            pproc.process_constituent(constituent, os.path.join(path, "entries"))
            print(f"...saving timeseries for {constituent}")
            pproc.__getattribute__(constituent).to_parquet(os.path.join(path, "entries"))

def aggregate_timeseries():
    import os
    import pandas as pd
    from tqdm import tqdm

    path = r"C://users/albert/documents/swmmpulse/entries"
    files = [file for file in os.listdir(path) if (file[:4] == "tss3") and (file.strip(".parquet")[-4:] == "vals")]

    fid = int(files[0].split("_")[2])
    tempseries = pd.read_parquet(os.path.join(path, files[0])).fillna(0.0).sum(axis=1).rename(str(fid))
    df_series = pd.DataFrame(tempseries)

    for file in tqdm(files[1:]):
        fid = int(file.split("_")[2])
        tempseries = pd.read_parquet(os.path.join(path, file)).fillna(0.0).sum(axis=1).rename(str(fid))
        df_series = df_series.join(tempseries)

    df_series.to_parquet(os.path.join(path, "aggregated_timeseries_ts3.parquet"), compression="LZ4")
    print("finished")

def extract_infected():
    import os
    import pandas as pd
    from tqdm import tqdm

    def create_names(fname):
        fid = int(fname.split("_")[2])
        pk = fname.replace("routetable", "packets")
        return fid, fname, pk

    path = r"C:\Users\albert\Documents\SWMMpulse\route_tables"
    files = [file for file in os.listdir(path) if (file[:2] == "s3") and (file.split(".")[0][-10:] == "routetable")]
    print(files)
    print("test")
    ls = []
    for file in tqdm(files):
        fid, rt, pk = create_names(file)
        df_packets = pd.read_parquet(os.path.join(path, pk))
        n_inf = df_packets.loc[df_packets["contents"].map(lambda x: x.__contains__("Cov-RNA"))]["pid"].count()
        n_tot = df_packets["pid"].count()
        r_inf = n_inf/n_tot
        ls.append({"sim_id": fid, "n_tot": n_tot, "n_inf": n_inf, "r_inf": r_inf})
    df = pd.DataFrame(ls)
    df.to_parquet(os.path.join(path, "infected_rates.parquet"), compression="LZ4")

if __name__ == "__main__":
    #test_sim()
    #mp_simulation()
    #aggregate_timeseries()
    extract_infected()

    print('finished')