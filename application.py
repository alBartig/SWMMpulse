from timeseries import QSeries
from pobject import PObject
from prouter import PRouter, Postprocessing, _Route_table
from network_lib import DirectedTree as ntwk, GraphConstants as gc
from environment import Environment
import os
from pconstants import Loading
import matplotlib.pyplot as plt

def prepare_environment():
    lpath = '/mnt/c/Users/albert/Documents/SWMMpulse/HS_calib_120_simp.out'
    #lpath = 'C:/Users/alber/Documents/swmm/swmmpulse/HS_calib_120_simp.out'
    qlut = QSeries(lpath)
    gpath = '/mnt/c/Users/albert/documents/SWMMpulse/HS_calib_120_simp/'
    #gpath = 'C:/Users/alber/documents/swmm/swmmpulse/HS_calib_120_simp/'
    graph = ntwk.from_directory(gpath)
    return graph,qlut

def generate_rts(graph, qlut):
    fractions = [0.00030,0.00010,0.00005,0.00003]
    for fraction in fractions:
        env = Environment()
        env.groups[0].weight = (1 - fraction)
        env.groups[1].weight = fraction

        for i in range(25):
            fname = f"rt_f{f'{fraction:.5f}'[-5:]}_{str(i).zfill(2)}.pickle"
            router = PRouter(graph=graph, qlookup=qlut)
            router.env = env
            routetable = router.route()
            routetable.to_file(f'C:/Users/alber/documents/swmm/swmmpulse/route_tables/{fname}')
            #routetable.to_file(f'/mnt/c/Users/albert/Documents/SWMMpulse/{fname}')
            #pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)

def evaluate_rts(graph, qlut, evalnode):
    path = 'C:/Users/alber/documents/swmm/swmmpulse/'
    files = [p for p in os.listdir(os.path.join(path, "route_tables")) if p[:2] == 'rt']
    for file in files:
        routetable = _Route_table.from_file(os.path.join(path,"route_tables",file))
        pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)
        eval_constituents = [Loading.FECAL, Loading.COV]
        for const in eval_constituents:
            entry_loc = os.path.join(path, 'entries', f"{file.rstrip('.pickle')}_{const}.pickle")
            pproc.process_constituent(const, entry_loc)
    print('finish')

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

def test_cov(graph, qlut, evalnode):
    path = 'C:/Users/alber/documents/swmm/swmmpulse'
    file = 'route_tables/rt_f00003_00.pickle'
    routetable = _Route_table.from_file(os.path.join(path,file))
    pproc = Postprocessing.from_rtable(routetable, evalnode, qlut, graph)
    pproc.process_constituent(Loading.COV, entry_loc=os.path.join(path,'entries/rt_f00003_00_Cov_RNA.pickle'), load=False)
    print('test_cov finished')


if __name__ == "__main__":
    evalnode = 'MH327-088-1'
    graph, qlut = prepare_environment()
    test_cov(graph, qlut, evalnode)

    print('finished')