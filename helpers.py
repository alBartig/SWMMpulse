import random
import numpy as np
import pandas as pd

def get_weighted_population(pop, fraction):
    n = 0
    for i in range(pop):
        p = random.random()
        if p <= fraction:
            n += 1
    return n

def series_to_dfmi(ls, timeindex):
    """
    turns list of dicts into multiindex dataframe
    dicts need to contain timeseries and other descriptors
    Args:
        ls (list): list of dicts
        timeindex: datetime index

    Returns:
        pd.DataFrame
    """
    l = len(timeindex)
    # check iterables and single values
    consts = []
    lists = []
    seltypes = [list, np.ndarray]
    for key, value in ls[0].items():
        if type(value) in seltypes and len(value) == l:
            lists.append(key)
        elif type(value) != list:
            consts.append(key)

    # create temporary list of rows
    temp = []
    for entry in ls:
        for i, time in enumerate(timeindex):
            row = [time]
            [row.append(entry[ls][i]) for ls in lists]
            [row.append(entry[key]) for key in consts]
            temp.append(row)

    df = pd.DataFrame(temp, columns=["time"] + lists + consts).set_index(consts + ["time"])
    return df