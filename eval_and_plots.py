import pandas as pd
import os
import datetime as dt
from sampler import Sampler, STRATEGY, STRATEGIES
import matplotlib.pyplot as plt
import statsmodels.api as sm
import numpy as np
import json
from pathlib import Path
from tqdm import tqdm
from win10toast import ToastNotifier
import seaborn as sns
toaster = ToastNotifier()

if __name__ == "__main__":
    fdir = Path(r"C:\Users\albert\Documents\SWMMpulse")
    path = Path(r"C:\Users\albert\Documents\SWMMpulse\02_processed\00_run1")

    files = os.listdir(path)
    flist = []
    for file in files:
        flist.append(file.replace(".", "_").split("_") + [file])
    dffiles = pd.DataFrame.from_records(flist,
                                        columns=["catchment", "infection rate", "step", "simID", "file-extension",
                                                 "filename"])
    dffiles["infection rate"] = dffiles["infection rate"].astype(int) * 10 ** (-3)
    dffiles.head()