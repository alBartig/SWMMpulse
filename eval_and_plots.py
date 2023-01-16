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
from matplotlib.ticker import AutoMinorLocator
toaster = ToastNotifier()

FDIR = Path(r"C:\Users\albert\Documents\SWMMpulse")
RUNNR = 2
DTINDEX = pd.date_range("2000-01-01", periods=8640, freq="10S")

PLOTS = [{"strategies": ["A", "I", "J"], # GRAB SAMPLES
          "plot_title": "Evaluation of grab sampling",
          "axtitles": ["Reference Regime", "Grab sampling: 9:00 am", "Grab sampling: 12:00 am"],
          "save_location": FDIR / "plots" / "grfk_grab_samples",
          "plot_name": "impact_grab_sampling",
          "ylim": [0,20]},
         {"strategies": ["A", "B", "C"], # WEIGHTING METHODS
          "plot_title": "sample concentrations over infection rates by weighting methods",
          "axtitles": ["time-weighted", "flow-weighted", "volume-weighted"],
          "save_location": FDIR/"plots"/"grfk_weighting_method",
          "plot_name":"imcact_weighting_method"},
         {"strategies": ["A", "D", "E"], # SAMPLE COUNT
          "plot_title": "sample concentrations over infection rates by sample count",
          "axtitles": ["24 samples / day", "48 samples / day", "72 samples / day"],
          "save_location": FDIR / "plots" / "grfk_samplecount",
          "plot_name": "impact_samplecount"},
         {"strategies": ["A", "F", "G"], # SAMPLING WINDOW
          "plot_title": "sample concentrations over infection rates by chosen time window",
          "axtitles": ["24 hours: 0:00 - 24:00", "12 hours: 4:00 - 16:00", "6 hours: 5:00 - 11:00"],
          "save_location": FDIR / "plots" / "grfk_sampling_period",
          "plot_name": "impact_sampling_period",
          "ylim": [0, 16]},
         {"strategies": ["A", "H", "K"], # SAMPLED TIME
          "plot_title": "sample concentrations over sampled time",
          "axtitles": ["24 min sampled time", "48 min sampled time", "72 min sampled time"],
          "save_location": FDIR / "plots" / "grfk_sampledtime",
          "plot_name": "impact_sampledtime"}]


def linreg_samples(df):
    # Prepare Data
    X = df["infection rate"].values
    X = sm.add_constant(X)
    y = df["concentration"].values
    # Fit model and extract data
    model = sm.OLS(y, X)
    results = model.fit()
    const, beta = results.params
    # calculate additional series prediction and errors
    df["pred_concentration"] = df["infection rate"] * beta + const
    df["errors"] = df["concentration"] - df["pred_concentration"]
    standard_error = np.std(df["errors"])
    pcc = df["infection rate"].corr(df["concentration"])
    df["pred_concentration_upper"] = df["pred_concentration"] + standard_error
    df["pred_concentration_lower"] = df["pred_concentration"] - standard_error
    return df, standard_error, pcc


def arange_eval_strategy(title=None):
    fig = plt.figure(figsize=[8.5, 4], facecolor="white", constrained_layout=True)
    axs = fig.subplot_mosaic([[0, 1, 2],
                              [3, 3, 3]],
                             gridspec_kw={'height_ratios': [3, 1]})
    fig.suptitle(title)
    return fig, axs


def prep_axs(ax, axttl=None, strategies=["A", "B", "C"], prefix="VC-"):
    ax[0].set(title=f"{prefix+strategies[0]:<4}: {axttl[0]:<20}", ylabel="sample concentrations [cp/l]")
    ax[1].set(title=f"{prefix+strategies[1]:<4}: {axttl[1]:<20}", xlabel="Fraction of shedding people in the catchment [-]")
    ax[2].set(title=f"{prefix+strategies[2]:<4}: {axttl[2]:<20}")

    ax[3].set(ylim=[0, 3], yticks=[0.5, 1.5, 2.5], yticklabels=["VC-"+s+":" for s in strategies])
    ax[3].yaxis.set_minor_locator(AutoMinorLocator(n=2))
    ax[3].grid(axis="y", which="minor", zorder=0)
    ax[3].grid(axis="x", which="major", zorder=0)
    return ax


def plot_samples(df, ax, **kwargs):
    # plot data
    # plot linear regression
    ax.plot("infection rate", "pred_concentration", data=df, linewidth=2, zorder=10, linestyle="solid")
    # plot samples
    #ax.scatter(x="infection rate", y="concentration", data=df, alpha=0.3, marker="+", s=30, zorder=5, linewidths=1)
    ax.scatter(x="infection rate", y="concentration", data=df, alpha=0.15, marker=".", s=5, zorder=5, linewidths=1)
    # plot area between standard error
    ax.fill_between(df["infection rate"], df["pred_concentration_upper"], df["pred_concentration_lower"],
                    alpha=0.2, zorder=0, color="lightcoral")
    # prepare legend
    # prepare text
    standard_error = np.std(df["errors"]) / df["concentration"].mean()
    pcc = df["infection rate"].corr(df["concentration"])
    textstr = f"{'RMSE:':<5}{standard_error:>6.2f}\n{'PCC:':<5}{pcc:>8.2f}"
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='square', facecolor='mistyrose', alpha=0.5)
    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=props, zorder=15)
    # additional ax settings
    ylim = kwargs.get("ylim", [0, 8])
    ax.set(xlim=[0.002,0.01], ylim=ylim, xticks=np.arange(0.002, 0.012, 0.002))


def prepare_sampler():
    sampler = Sampler()
    flows = pd.read_csv(r"C:\Users\albert\Documents\SWMMpulse\Flow\com2a_hydraulic.csv")
    flows["time"] = DTINDEX
    flows = flows.set_index("time")["flow-rate [l/s]"]

    sampler.add_flows(flows)
    sampler.flows.head()
    return sampler


def read_processed_file_directory():
    files = os.listdir(FDIR / "02_processed" / f"00_run{RUNNR}")
    flist = []
    for file in files:
        flist.append(file.replace(".", "_").split("_") + [file])
    dffiles = pd.DataFrame.from_records(flist,
                                        columns=["catchment", "step", "simID", "file-extension",
                                                 "filename"])
    #dffiles["infection rate"] = dffiles["infection rate"].astype(int) * 10 ** (-3)
    return dffiles


def calculate_infected_rates(dffiles):
    packet_dir = FDIR / "00_packets" / f"00_run{RUNNR}"
    infected_rates = []
    # iterate over filenames
    for fname in tqdm(dffiles.loc[:,"filename"].values, desc="Calculating true infection rate..."):
            # modify filename so it matches packet-files
            fname = fname.replace("processed", "packets")
            # read each packet-file
            df_p = pd.read_json(packet_dir/fname)
            # get number of packets in packet-file
            rowcount, colcount = df_p.shape
            # count packets classified as "infected"
            infectedcount = df_p["classification"].value_counts()["infected"]
            # append calculated infected rate
            infected_rates.append(infectedcount/rowcount)
    # append new column to dffiles
    dffiles["infection rate"] = infected_rates
    toaster.show_toast("finished calculating infection rates")
    return dffiles


def read_timeseries(dffiles):
    timeseries_dir = FDIR / "02_processed" / f"00_run{RUNNR}"
    # dtindex = pd.date_range("2000-01-01", periods=8640, freq="10S")
    df_timeseries = pd.DataFrame(index=DTINDEX)

    for simID, tsfile in tqdm(dffiles.loc[:, ["simID", "filename"]].values, desc="reading timeseries files..."):
        df_timeseries = df_timeseries.join(pd.read_json(timeseries_dir / tsfile, typ="series").rename(simID))

    print("----------------reading timeseries finished---------------------")
    toaster.show_toast("Timeseries read", "All timeseries have been read into memory")
    return df_timeseries


def eval_plot(dffiles, df_timeseries, sampler, plot_data=None):
    title = plot_data.get("plot_title")
    axtitles = plot_data.get("axtitles")
    strategies = plot_data.get("strategies")
    save_loc = plot_data.get("save_location")
    plot_name = plot_data.get("plot_name")

    fig, axs = arange_eval_strategy(title)

    for k, name in enumerate(strategies):
        # get strategy
        strategy = STRATEGIES.get(name)
        # sample timeseries
        samples, smeta = sampler.sample(df_timeseries, strategy)
        df_samples = pd.DataFrame(samples)
        # link sample concentrations to corresponding infection rates
        df_samples = df_samples.join(dffiles.set_index("simID")["infection rate"])
        # store samples, sampletimes and weights in parquet files
        df_samples.to_parquet(save_loc / f"{plot_name}_{name}_samples.parquet")
        df_s = pd.DataFrame(smeta.get(STRATEGY.SAMPLETIMES)).join(smeta.get(STRATEGY.SAMPLEWEIGHTS))
        df_s.columns = df_s.columns.astype("string")
        df_s.to_parquet(save_loc / f"{plot_name}_{name}_strategy.parquet")
        # regression analysis
        df_samples, standard_error, pcc = linreg_samples(df_samples)
        # plot data
        plot_samples(df_samples, axs[k])
        # plot strategy
        sampler.plot_strategy(strategy, axs[3], y0=k, legend=False)
        axs[3].set(ylim=[0, 3], yticks=[0.5, 1.5, 2.5], yticklabels=["a)", "b)", "c)"])
        axs[3].yaxis.set_minor_locator(AutoMinorLocator(n=2))
        axs[3].grid(axis="y", which="minor")
        axs[3].grid(axis="x", which="major")

    prep_axs(axs, axtitles, strategies)

    plt.savefig(save_loc / f"{plot_name}.png", dpi=300)
    plt.savefig(save_loc / f"{plot_name}.svg")
    return None


def replot_data(sampler, plot_data=None):
    """
    Uses the saved data to plot evaluation plots again
    Args:
        sampler (Sampler):
        plot_data (dict):

    Returns: None

    """
    title = plot_data.get("plot_title")
    axtitles = plot_data.get("axtitles")
    strategies = plot_data.get("strategies")
    save_loc = plot_data.get("save_location")
    plot_name = plot_data.get("plot_name")

    fig, axs = arange_eval_strategy(title)

    for k, name in enumerate(strategies):
        # get strategy
        strategy = STRATEGIES.get(name)
        df_s = pd.read_parquet(save_loc / f"{plot_name}_{name}_strategy.parquet")
        df_samples = pd.read_parquet(save_loc / f"{plot_name}_{name}_samples.parquet")
        df_samples, standard_error, pcc = linreg_samples(df_samples)
        # plot data
        plot_samples(df_samples, axs[k], **plot_data)
        # plot strategy
        sampler.plot_strategy(strategy, axs[3], y0=k, legend=False)

    prep_axs(axs, axtitles, strategies)

    plt.savefig(save_loc / f"replotted_{plot_name}.png", dpi=300)
    plt.savefig(save_loc / f"replotted_{plot_name}.svg")
    return None


def eval_weighting(dffiles, df_timeseries, sampler):
    selected_strategies = ["A", "B", "C"]
    axtitles = ["time-weighted", "flow-weighted", "volume-weighted"]

    plot_data = {"strategies": selected_strategies,
                 "plot_title": "sample concentrations over infection rates by weighting methods",
                 "axtitles": axtitles,
                 "save_location": FDIR/"plots"/"grfk_weighting_method",
                 "plot_name":"imcact_weighting_method"}

    eval_plot(dffiles, df_timeseries, sampler, plot_data=plot_data)  # Marlene
    return None


def eval_samplecount(dffiles, df_timeseries, sampler):
    selected_strategies = ["A", "D", "E"]
    axtitles = ["24 samples / day", "48 samples / day", "72 samples / day"]
    title = "sample concentrations over infection rates by sample count"

    plot_data = {"strategies": selected_strategies,
                 "plot_title": title,
                 "axtitles": axtitles,
                 "save_location": FDIR / "plots" / "grfk_samplecount",
                 "plot_name": "impact_samplecount"}

    eval_plot(dffiles, df_timeseries, sampler, plot_data=plot_data)
    return None


def eval_samplingwindow(dffiles, df_timeseries, sampler):
    title = "sample concentrations over infection rates by chosen time window"
    axtitles = ["24 hours: 0:00 - 24:00", "12 hours: 4:00 - 16:00", "6 hours: 5:00 - 11:00"]
    selected_strategies = ["A", "F", "G"]

    plot_data = {"strategies": selected_strategies,
                 "plot_title": title,
                 "axtitles": axtitles,
                 "save_location": FDIR / "plots" / "grfk_sampling_period",
                 "plot_name": "impact_sampling_period"}

    eval_plot(dffiles, df_timeseries, sampler, plot_data=plot_data)
    return None


def eval_grab_sample(dffiles, df_timeseries, sampler):
    title = "Evaluation of grab sampling"
    selected_strategies = ["A", "I", "J"]
    axtitles = ["Reference Regime", "Grab sampling: 9:00 am", "Grab sampling: 12:00 am"]

    plot_data = {"strategies": selected_strategies,
                 "plot_title": title,
                 "axtitles": axtitles,
                 "save_location": FDIR / "plots" / "grfk_grab_samples",
                 "plot_name": "impact_grab_sampling"}

    eval_plot(dffiles, df_timeseries, sampler, plot_data=plot_data)
    return None


def eval_sampledtime(dffiles, df_timeseries, sampler):
    axtitles = ["24 min sampled time", "48 min sampled time", "72 min sampled time"]
    title = "sample concentrations over sampled time"
    selected_strategies = ["A", "H", "K"]

    plot_data = {"strategies": selected_strategies,
                 "plot_title": title,
                 "axtitles": axtitles,
                 "save_location": FDIR / "plots" / "grfk_sampledtime",
                 "plot_name": "impact_sampledtime"}

    eval_plot(dffiles, df_timeseries, sampler, plot_data=plot_data)
    return None


def evaluate_sim_data():
    # read filelist with processed timeseries
    dffiles = read_processed_file_directory()
    # calculate actual fraction of shedding population
    dffiles = calculate_infected_rates(dffiles)
    # read timeseries files
    df_timeseries = read_timeseries(dffiles)

    # prepare sampler
    sampler = prepare_sampler()
    # iterate through plotting dicts and evaluate each
    for plot_dict in PLOTS:
        eval_plot(dffiles, df_timeseries, plot_dict)

    return None


def replot_evaluations():
    # prepare sampler
    sampler = prepare_sampler()
    # iterate through plotting dicts and replot evaluations
    for plot_dict in PLOTS:
        replot_data(sampler, plot_data=plot_dict)

    return None

def main():
    replot_evaluations()


if __name__ == "__main__":
    main()
