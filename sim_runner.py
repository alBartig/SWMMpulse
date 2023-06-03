import pandas as pd
import numpy as np
from environment import Environment, DirectedTree, DEFAULT, CONSTITUENT
from prouter import Router
import json
from tqdm import tqdm
from pathlib import Path


def preparation_env(p_swmm_out=None, p_swmm_inp=None, p_pop_data=None):
    """
    Prepares the environment with given pathes
    Args:
        p_swmm_out (str): Path of the swmm-out-file used for the hydraulic information
        p_swmm_inp (str): Path of the swmm-inp-file used to construct the graph
        p_pop_data (str): Path of the csv-file containing the population data

    Returns:
        Environment
    """
    env = Environment()
    # read hydraulics from swmm-inp and add to env
    if not p_swmm_out:
        p_swmm_out = r"C:\Users\albert\Documents\SWMMpulse\05_swmm_model\02_sensitivity_analysis\HS_calib_120_simp_005lps.out"
    # env.read_swmmoutfile(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp.out")
    env.read_swmmoutfile(p_swmm_out)
    # read nodes and edges from swmm-out and construct graph
    if not p_swmm_inp:
        p_swmm_inp = r"C:\Users\albert\Documents\SWMMpulse\05_swmm_model\02_sensitivity_analysis\HS_calib_120_simp_005lps.inp"
    # graph = DirectedTree.from_swmm(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp.inp")
    graph = DirectedTree.from_swmm(p_swmm_inp)
    # read pop-data from csv-file and add to graph
    if not p_pop_data:
        p_pop_data = r"pop_node_data.csv"
    # node_data = pd.read_csv(r"C:\Users\albert/Documents/SWMMpulse/HS_calib_120_simp/pop_node_data.csv")
    node_data = pd.read_csv(p_pop_data)
    node_data = node_data.set_index("NAME").to_dict(orient="index")
    graph.add_nodevalues(node_data)
    # add graph to environment
    env.add_graph(graph)
    return env


def set_infectivity(r_inf, env):
    """
    sets fraction of shedding population in an environment env to r_inf
    Args:
        r_inf (float): fraction of RNA-shedding population
        env (Environment): environment to change

    Returns:
        Environment
    """
    env.information["groups"][0]["weight"] = 1 - r_inf
    env.information["groups"][1]["weight"] = r_inf
    return env


def monte_carlo_variable_shedding_fraction(n=None, sim_dir=None):
    """
    n Simulations are run for infection rate scenarios between a start rate end an end rate. in between values are
    interpolated linearly by simulation id.
    Args:
        n: number of simulations run

    Returns:

    """
    if not n:
        n = 10000
    # create subdirectories
    (sim_dir/"00_packets").mkdir(parents=True, exist_ok=True)
    (sim_dir/"01_routetables").mkdir(parents=True, exist_ok=True)
    (sim_dir/"02_processed").mkdir(parents=True, exist_ok=True)

    # runnr = 0
    # run = f"00_run{runnr - 1:d}"
    # create router-object
    router = Router()
    # prepare simulation-environment
    env = preparation_env()
    env.information["dispersion_rate"] = 0.16
    groot = env.graph.root

    cov = env.information["constituents"]["Cov-RNA"]
    cov["fractions"] = [1]
    cov["skewedness"] = [[1, 1]]
    env.information["constituents"]["Cov-RNA"] = cov

    # bins = pd.interval_range(start=0, end=7200, freq=1200)
    # bins = bins.append(pd.interval_range(start=7200, end=14400, freq=3600))

    # prepare simulation counter
    simid = 0
    # run n-simulations
    for i in tqdm(range(n), desc="Running scenarios..."):
        # calculate fraction of shedding population of current simulation
        r_inf = np.around(0.0002 + i * (0.01 - 0.0002) / n, 4)

        # update environment with new infectivity
        env = set_infectivity(r_inf, env)
        router.add_environment(env)

        # generate packets
        df_packets = router.environment.get_packets()
        fname = f"HS120_packets_{simid:06d}.parquet"
        df_packets.to_parquet(sim_dir / "00_packets" / fname)

        # generate routetable
        dict_routetable = router.route(packets=df_packets, as_dataframe=False)

        # write routetable at network root to file
        fname = f"HS120_routetable_{simid:06d}.json"
        with open(sim_dir / "01_routetables" / fname, "w") as jobj:
            jobj.write(json.dumps(dict_routetable[groot]))

        # calculate timeseries data at network root
        dict_processed = router.postprocess(dict_routetable, df_packets, env.information["constituents"].get(CONSTITUENT.COV))
        fname = f"HS120_processed_{simid:06d}.json"

        # create dataframe
        df_timeseries = pd.DataFrame(dict_processed, index=pd.date_range("2000-01-01", freq="10S", periods=8640))
        # calculate timeseries for entire catchment
        df_timeseries.sum(axis=1).to_json(sim_dir / "02_processed" / fname)

        simid += 1
    print("----------------------- MonteCarlo runs finished ----------------------------------")

def main():
    sim_dir = Path(r"C:\Users\albert\Documents\SWMMpulse\03_sensitivity_analysis\03_avgflow_005")
    monte_carlo_variable_shedding_fraction(sim_dir=sim_dir, n=10000)
    pass


if __name__ == '__main__':
    main()
