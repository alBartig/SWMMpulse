import pandas as pd
import os
import datetime as dt
from sampler import Sampler, STRATEGY, STRATEGIES
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import statsmodels.api as sm
import numpy as np
import json
from pathlib import Path
from tqdm import tqdm
from win10toast import ToastNotifier
import seaborn as sns
from matplotlib.ticker import AutoMinorLocator
toaster = ToastNotifier()

RUNNR = 2
DTINDEX = pd.date_range("2000-01-01", periods=8640, freq="10S")

PLOTS = [{"strategies": ["A", "J", "K"], # GRAB SAMPLES
          "plot_title": "Evaluation of grab sampling",
          "axtitles": ["Reference Regime", "Grab sampling: 9:00 am", "Grab sampling: 12:00 am"],
          "plot_name": "impact_grab_sampling",
          "ylim": [0,20]},
         {"strategies": ["A", "B", "C"], # WEIGHTING METHODS
          "plot_title": "sample concentrations over infection rates by weighting methods",
          "axtitles": ["time-weighted", "flow-weighted", "volume-weighted"],
          "plot_name":"imcact_weighting_method"},
         {"strategies": ["A", "D", "E"], # SAMPLE COUNT
          "plot_title": "sample concentrations over infection rates by sample count",
          "axtitles": ["24 samples / day", "48 samples / day", "72 samples / day"],
          "plot_name": "impact_samplecount",
          "scaling":[1, 0.5, 1/3]},
         {"strategies": ["A", "F", "G"], # SAMPLING WINDOW
          "plot_title": "sample concentrations over infection rates by chosen time window",
          "axtitles": ["24 hours: 0:00 - 24:00", "12 hours: 4:00 - 16:00", "6 hours: 5:00 - 11:00"],
          "plot_name": "impact_sampling_period",
          "ylim": [0, 16]},
         {"strategies": ["A", "H", "I"], # SAMPLED TIME
          "plot_title": "sample concentrations over sampled time",
          "axtitles": ["24 min sampled time", "48 min sampled time", "72 min sampled time"],
          "plot_name": "impact_sampledtime",
          "scaling": [1/3, 2/3, 1]},
         {"strategies": ["A", "X", "Y"],  # SAMPLED TIME
          "plot_title": "Best results witout investing additional resources",
          "axtitles": ["Reference scenario", "Intermediate scenario", "Most elaborate scenario"],
          "plot_name": "impact_scenarios",
          "scaling": [1, 1/2, 1/3]}
         ]




def arange_eval_strategy(**kwargs):
    title = kwargs.get("title")
    figsize = kwargs.get("figsize", [8.5, 4])
    fig = plt.figure(figsize=figsize, facecolor="white", constrained_layout=True)
    axs = fig.subplot_mosaic([[0, 1, 2],
                              [3, 3, 3]],
                             gridspec_kw={'height_ratios': [3, 1]})
    fig.suptitle(title)
    return fig, axs


def prep_axs(ax, axttl=None, strategies=["A", "B", "C"], prefix="VC-"):
    ax[0].set(title=f"{prefix+strategies[0]:<4}: {axttl[0]:<20}", ylabel="sample concentrations [gc/L]")
    ax[1].set(title=f"{prefix+strategies[1]:<4}: {axttl[1]:<20}", xlabel="Fraction of people shedding the virus in the catchment [-]")
    ax[2].set(title=f"{prefix+strategies[2]:<4}: {axttl[2]:<20}")
    #ax[0].grid(zorder=0, color="lightgrey")

    ax[3].set(ylim=[0, 3], yticks=[0.5, 1.5, 2.5], yticklabels=["VC-"+s+":" for s in strategies],
              xlabel="Time [HH:MM]", title="Sampling strategy visualization")
    ax[3].yaxis.set_minor_locator(AutoMinorLocator(n=2))
    ax[3].grid(axis="y", which="minor", color="lightgrey", zorder=0)
    ax[3].grid(axis="x", which="major", color="lightgrey", zorder=0)
    return ax


def plot_samples(df, ax, **kwargs):
    # plot data
    # plot linear regression
    ax.plot("infection rate", "pred_concentration", data=df, linewidth=2, zorder=10, linestyle="solid")
    # plot samples
    #ax.scatter(x="infection rate", y="concentration", data=df, alpha=0.3, marker="+", s=30, zorder=5, linewidths=1)
    ax.scatter(x="infection rate", y="concentration", data=df, alpha=0.13,
               marker=".", s=5, zorder=0, linewidths=1, rasterized=True)
    # plot area between standard error
    ax.fill_between(df["infection rate"], df["pred_concentration_upper"], df["pred_concentration_lower"],
                    alpha=0.3, zorder=5, color="lightcoral")
    # prepare legend
    # prepare text
    standard_error = np.std(df["errors"]) / df["concentration"].mean()
    pcc = df["infection rate"].corr(df["concentration"])
    textstr = f"{'NRMSE:':<6}{standard_error:>6.2f}\n{'PCC:':<7}{pcc:>8.2f}"
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='square', facecolor='mistyrose', alpha=0.5)
    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=props, zorder=15)
    # additional ax settings
    ylim = kwargs.get("ylim", [0, 8])
    ax.set(xlim=[0.0002,0.01], ylim=ylim, xticks=[0.0002]+np.arange(0.002, 0.012, 0.002).tolist())

    def small_formatter(x, pos):
        return f"{x*100:3.2f}"

    ax.xaxis.set_major_formatter(small_formatter)


def eval_plot(dffiles, df_timeseries, sampler, plot_data=None):
    #title = plot_data.get("plot_title")
    axtitles = plot_data.get("axtitles")
    strategies = plot_data.get("strategies")
    save_loc = plot_data.get("save_location")
    plot_name = plot_data.get("plot_name")
    scaling = plot_data.get("scaling", [1, 1, 1])

    plot_dir = save_loc / plot_name
    plot_dir.mkdir(exist_ok=True, parents=True)

    #fig, axs = arange_eval_strategy(title)
    fig, axs = arange_eval_strategy()

    for k, (name, scalar) in enumerate(zip(strategies, scaling)):
        # get strategy
        strategy = STRATEGIES.get(name)
        # sample timeseries
        samples, smeta = sampler.sample(df_timeseries, strategy)
        df_samples = pd.DataFrame(samples)
        # link sample concentrations to corresponding infection rates
        df_samples = df_samples.join(dffiles.set_index("simID")["infection rate"])
        # store samples, sampletimes and weights in parquet files
        df_samples.to_parquet(plot_dir / f"{plot_name}_{name}_samples.parquet")
        df_s = pd.DataFrame(smeta.get(STRATEGY.SAMPLETIMES)).join(smeta.get(STRATEGY.SAMPLEWEIGHTS))
        df_s.columns = df_s.columns.astype("string")
        df_s.to_parquet(plot_dir / f"{plot_name}_{name}_strategy.parquet")
        # regression analysis
        df_samples, standard_error, pcc = _linreg_samples(df_samples)
        # plot data
        plot_samples(df_samples, axs[k])
        # plot strategy
        sampler.plot_strategy(strategy, axs[3], y0=k, legend=False, scaling=scalar)
        axs[3].set(ylim=[0, 3], yticks=[0.5, 1.5, 2.5], yticklabels=["a)", "b)", "c)"])
        axs[3].yaxis.set_minor_locator(AutoMinorLocator(n=2))
        axs[3].grid(axis="y", which="minor")
        axs[3].grid(axis="x", which="major")

    prep_axs(axs, axtitles, strategies)
    plt.savefig(plot_dir / f"{plot_name}.png", dpi=700)
    plt.savefig(plot_dir / f"{plot_name}.svg")
    return None


def replot_data(sampler, plot_data=None, **grafics_dict):
    """
    Uses the saved data to plot evaluation plots again
    Args:
        sampler (Sampler):
        plot_data (dict):

    Returns: None

    """
    #title = plot_data.get("plot_title")
    axtitles = plot_data.get("axtitles")
    strategies = plot_data.get("strategies")
    save_loc = plot_data.get("save_location")
    plot_name = plot_data.get("plot_name")
    scaling = plot_data.get("scaling", [1, 1, 1])

    #fig, axs = arange_eval_strategy(title=title, **grafics_dict)
    fig, axs = arange_eval_strategy(**grafics_dict)

    for k, (name, scalar) in enumerate(zip(strategies, scaling)):
        # get strategy
        strategy = STRATEGIES.get(name)
        df_s = pd.read_parquet(save_loc / f"{plot_name}_{name}_strategy.parquet")
        df_samples = pd.read_parquet(save_loc / f"{plot_name}_{name}_samples.parquet")
        df_samples, standard_error, pcc = _linreg_samples(df_samples)
        # plot data
        plot_samples(df_samples, axs[k], **plot_data)
        # plot strategy
        sampler.plot_strategy(strategy, axs[3], y0=k, scaling=scalar, legend=False)

    prep_axs(axs, axtitles, strategies)

    plt.savefig(save_loc / f"replotted_{plot_name}.png", dpi=grafics_dict.get("dpi", 300))
    plt.savefig(save_loc / f"replotted_{plot_name}.svg")
    return None


def _read_timeseries(dffiles, data_dir):
    # timeseries_dir = data_dir / "02_processed" / f"00_run{RUNNR}"
    timeseries_dir = data_dir / "02_processed"
    # dtindex = pd.date_range("2000-01-01", periods=8640, freq="10S")
    df_timeseries = pd.DataFrame(index=DTINDEX)

    for simID, tsfile in tqdm(dffiles.loc[:, ["simID", "filename"]].values, desc="reading timeseries files..."):
        df_timeseries = df_timeseries.join(pd.read_json(timeseries_dir / tsfile, typ="series").rename(simID))

    print("----------------reading timeseries finished---------------------")
    # toaster.show_toast("Timeseries read", "All timeseries have been read into memory")
    return df_timeseries


def _calculate_infected_rates(dffiles, data_dir):
    # packet_dir = data_dir / "00_packets" / f"00_run{RUNNR}"
    packet_dir = data_dir / "00_packets"
    infected_rates = []
    # iterate over filenames
    for fname in tqdm(dffiles.loc[:,"filename"].values, desc="Calculating true infection rate..."):
            # modify filename so it matches packet-files
            fname = fname.replace("processed", "packets")
            fname = fname.replace("json", "parquet")
            # read each packet-file
            # df_p = pd.read_json(packet_dir/fname)
            df_p = pd.read_parquet(packet_dir/fname)
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


def _read_processed_file_directory(data_dir):
    # files = os.listdir(data_dir / "02_processed" / f"00_run{RUNNR}")
    files = os.listdir(data_dir / "02_processed")
    flist = []
    for file in files:
        flist.append(file.replace(".", "_").split("_") + [file])
    dffiles = pd.DataFrame.from_records(flist,
                                        columns=["catchment", "step", "simID", "file-extension",
                                                 "filename"])
    #dffiles["infection rate"] = dffiles["infection rate"].astype(int) * 10 ** (-3)
    return dffiles


def prepare_filelist(data_dir):
    """
    Scans all file in data-directory and extracts metadata and fractions of RNA-shedding-population and stores in df
    Args:
        data_dir (Path): Path where simulation data is to be found

    Returns:
        None
    """
    plot_dir = data_dir / "04_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    # read filelist with processed timeseries
    df_files = _read_processed_file_directory(data_dir=data_dir)
    # calculate actual fraction of shedding population
    df_files = _calculate_infected_rates(df_files, data_dir=data_dir)
    df_files.to_parquet(plot_dir / "df_files.parquet")
    return


def prepare_sampler(flow_dir=None):
    """
    Ceates sampler object with hydraulic information
    Args:
        flow_dir (Path): Path where hydraulic information lies

    Returns:
        Sampler
    """
    if not flow_dir:
        # flow_dir = Path(r"C:\Users\albert\Documents\SWMMpulse\05_swmm_model\02_sensitivity_analysis\com2a_hydraulic_005lps.csv")
        flow_dir = Path(r"C:\Users\albert\Documents\SWMMpulse\05_swmm_model\00_flow_analysis\com2a_hydraulic_smooth.csv")
    sampler = Sampler()
    try:
        flows = pd.read_csv(flow_dir)
        flows["time"] = DTINDEX
        flows = flows.set_index("time")["flow-rate [l/s]"]
    except:
        flows = pd.read_csv(flow_dir, delimiter=";")
        flows["time"] = DTINDEX
        flows = flows.set_index("time")["flow-rate [l/s]"]

    sampler.add_flows(flows)
    sampler.flows.head()
    return sampler


def _linreg_samples(df):
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


def evaluate_strategy(dffiles, df_timeseries, sampler, strategy, save_location):
    """
    samples timeseries and returns samples and strategies
    Args:
        dffiles (pd.DataFrame): Dataframe of files
        df_timeseries (pd.DataFrame):
        sampler (Sampler): sampler object
        strategy (str): name of strategy ("A", "B", ...)

    Returns:
        df_samples (pd.DataFrame), df_strategy (pd.DataFrame)
    """
    # get strategy-dict
    strategy_dict = STRATEGIES.get(strategy)
    # sample timeseries
    samples, smeta = sampler.sample(df_timeseries, strategy_dict)
    df_samples = pd.DataFrame(samples)
    # link sample concentrations to corresponding infection rates
    df_samples = df_samples.join(dffiles.set_index("simID")["infection rate"])
    df_samples, standard_error, pcc = _linreg_samples(df_samples)
    df_samples.to_parquet(save_location / f"Strategy_{strategy}_samples.parquet")
    # store samples, sampletimes and weights in parquet files

    df_strategy = pd.DataFrame(smeta.get(STRATEGY.SAMPLETIMES)).join(smeta.get(STRATEGY.SAMPLEWEIGHTS))
    df_strategy.columns = df_strategy.columns.astype("string")
    df_strategy.to_parquet(save_location / f"Strategy_{strategy}_strategy.parquet")

    return df_samples, df_strategy


def evaluate_scenarios(data_dir, plots):
    """
    Calculates strategies and sample concentrations for the strategies given in plots
    Args:
        data_dir (Path): Path where the data to be evaluated is to be found
        plots (list): List of plots with specifications for each plot (mainly the strategies to be plotted)

    Returns:
        None
    """
    plot_dir = data_dir / "04_plots"
    # read filelist with processed timeseries
    df_files = pd.read_parquet(plot_dir / "df_files.parquet")
    # prepare sampler
    sampler = prepare_sampler()
    # read timeseries files
    df_timeseries = _read_timeseries(df_files, data_dir=data_dir)

    # iterate through plotting dicts and evaluate each
    for plot_dict in plots:
        strategies = plot_dict.get("strategies")

        for strategy_name in strategies:
            df_samples, df_strategy = evaluate_strategy(df_files, df_timeseries, sampler, strategy_name,
                                                        plot_dir)
            df_samples.to_parquet(plot_dir / f"strategy_{strategy_name}_samples.parquet")
            df_strategy.to_parquet(plot_dir / f"strategy_{strategy_name}_strategy.parquet")
    return


def plot_single_comparison(data_dict, sampler):
    axtitles = data_dict.get("axtitles")
    strategies = data_dict.get("strategies")
    plot_dir = data_dict.get("save_location")
    plot_name = data_dict.get("plot_name")
    scaling = data_dict.get("scaling", [1, 1, 1])

    # ---------------------------------------- 1. PART - PREPARE AXES -------------------------
    def prepare_comparison_plot(**kwargs):
        title = kwargs.get("title")
        figsize = kwargs.get("figsize", [8.5, 4])
        fig = plt.figure(figsize=figsize, facecolor="white", constrained_layout=True)
        axs = fig.subplot_mosaic([[0, 1, 2],
                                  [3, 3, 3]],
                                 gridspec_kw={'height_ratios': [3, 1]})
        fig.suptitle(title)
        return fig, axs

    fig, axs = prepare_comparison_plot()

    # ----------------------------------------- 2. PART - PLOT DATA ----------------------------
    def plot_model_runs(df, ax, **kwargs):
        # plot data
        # plot linear regression
        ax.plot("infection rate", "pred_concentration", data=df, linewidth=2, zorder=10, linestyle="solid")
        # plot samples
        # ax.scatter(x="infection rate", y="concentration", data=df, alpha=0.3, marker="+", s=30, zorder=5, linewidths=1)
        ax.scatter(x="infection rate", y="concentration", data=df, alpha=0.13,
                   marker=".", s=5, zorder=0, linewidths=1, rasterized=True)
        # plot area between standard error
        ax.fill_between(df["infection rate"], df["pred_concentration_upper"], df["pred_concentration_lower"],
                        alpha=0.3, zorder=5, color="lightcoral")
        # prepare legend
        # prepare text
        standard_error = np.std(df["errors"]) / df["concentration"].mean()
        pcc = df["infection rate"].corr(df["concentration"])
        textstr = f"{'NRMSE:':<6}{standard_error:>6.2f}\n{'PCC:':<7}{pcc:>8.2f}"
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='square', facecolor='mistyrose', alpha=0.5)
        # place a text box in upper left in axes coords
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=props, zorder=15)
        # additional ax settings
        # ylim = kwargs.get("ylim", [0, 8])
        ax.set(xlim=[0.0002, 0.01], xticks=[0.0002] + np.arange(0.002, 0.012, 0.002).tolist())
        # ax.set(ylim=ylim)

        def small_formatter(x, pos):
            return f"{x * 100:3.2f}"

        ax.xaxis.set_major_formatter(small_formatter)
        return

    for k, (strategy_name, scalar) in enumerate(zip(strategies, scaling)):
        df_samples = data_dict.get(strategy_name).get("df_samples")
        # plot data
        plot_model_runs(df_samples, axs[k])
        # plot strategy
        sampler.plot_strategy(STRATEGIES.get(strategy_name), axs[3], y0=k, legend=False, scaling=scalar)

    ylim_max = max([axs[k].get_ylim()[1] for k in range(len(strategies))])
    for k in range(len(strategies)):
        axs[k].set(ylim=[0, ylim_max])


    # --------------------------------------------- 3. PART - FINISH AXES --------------------------
    def adjust_comparison_plot_axis(ax, axttl=None, strategies=None, prefix="VC-"):
        ax[0].set(title=f"{prefix + strategies[0]:<4}: {axttl[0]:<20}", ylabel="sample concentrations [gc/L]")
        ax[1].set(title=f"{prefix + strategies[1]:<4}: {axttl[1]:<20}",
                  xlabel="Percentage of people shedding the virus in the catchment [%]")
        ax[2].set(title=f"{prefix + strategies[2]:<4}: {axttl[2]:<20}")
        # ax[0].grid(zorder=0, color="lightgrey")

        ax[3].set(ylim=[0, 3], yticks=[0.5, 1.5, 2.5], yticklabels=["VC-" + s + ":" for s in strategies],
                  xlabel="Time [HH:MM]", title="Sampling strategy visualization", ylabel="Campaigns")
        ax[3].yaxis.set_minor_locator(AutoMinorLocator(n=2))
        ax[3].grid(axis="y", which="minor", color="lightgrey", zorder=0)
        ax[3].grid(axis="x", which="major", color="lightgrey", zorder=0)
        return ax

    adjust_comparison_plot_axis(axs, axtitles, strategies)

    # --------------------------------------------- 4. PART - SAVE PLOTS --------------------------
    plt.savefig(plot_dir / f"{plot_name}.png", dpi=700)
    plt.savefig(plot_dir / f"{plot_name}.svg")
    return


def plot_scenarios(data_dir, plots):
    """
    Plots the data from evaluate_scenarios
    Args:
        data_dir (Path): Path, where the subdirectory "04_plots" should be created
        plots (list): list of dictionaries, constant in script

    Returns:
        None
    """
    plot_dir = data_dir / "04_plots"

    # prepare sampler
    sampler = prepare_sampler()

    for plot_dict in plots:
        strategies = plot_dict.get("strategies")
        plot_dict["save_location"] = plot_dir

        # load information for plot/comparison
        for strategy_name in strategies:
            plot_dict[strategy_name] = {}
            plot_dict[strategy_name]["df_samples"] = pd.read_parquet(plot_dir / f"strategy_{strategy_name}_samples.parquet")
            plot_dict[strategy_name]["df_strategy"] = pd.read_parquet(plot_dir / f"strategy_{strategy_name}_strategy.parquet")

        plot_single_comparison(plot_dict, sampler)

    return None


def replot_evaluations():
    # prepare sampler
    sampler = prepare_sampler()
    # iterate through plotting dicts and replot evaluations
    for plot_dict in PLOTS:
        replot_data(sampler, plot_data=plot_dict)

    return None


def plot_cumulative_pollutograph(data_dir):
    from environment import DEFAULT
    import matplotlib.dates as mdates

    plot_dir = data_dir / "04_plots"
    # read inflows
    sampler = prepare_sampler()
    # read filelist with processed timeseries
    df_files = pd.read_parquet(plot_dir / "df_files.parquet")
    # read timeseries files
    df_timeseries = _read_timeseries(df_files.iloc[::10,:], data_dir=data_dir)
    # prepare plot
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=[8.5, 4], constrained_layout=True)
    axs2 = axs[0].twinx()

    # plot DWF
    axs[0].plot(sampler.flows, color="tab:purple")
    axs[0].set(title="a) Input patterns", ylabel="Dry-Weather-Flow [L/s]")
    # assign date locator / formatter to the x-axis to get proper labels
    # locator = mdates.HourLocator(interval=6)
    locator = mdates.AutoDateLocator(minticks=4, maxticks=5)
    locator.intervald[4] = [6]
    formatter = mdates.DateFormatter('%H:%M')  # locator
    axs[0].xaxis.set_major_locator(locator)
    axs[0].xaxis.set_major_formatter(formatter)
    axs[0].xaxis.set_minor_locator(AutoMinorLocator(n=3))
    axs[0].set(xlim=[mdates.date2num(pd.to_datetime("2000-01-01 00:00:00")),
                     mdates.date2num(pd.to_datetime("2000-01-02 00:00:00"))],
               xlabel="Time [HH:MM]")
    axs[0].set(xlim=[pd.to_datetime("2000-01-01 00:00:00"), pd.to_datetime("2000-01-02 00:00:00")],
               ylim=[0, 350], yticks=np.arange(0, 350, 50))
    axs[0].grid()

    # plot defecation pattern
    t = pd.date_range("2000-01-01 00:30:00", "2000-01-01 23:30:00", periods=24)
    axs2.bar(x=t, height=DEFAULT.PATTERN_BRISTOL, width=0.04, color="tab:brown")
    axs2.set(ylabel="Temporal defecation distribution [%]",
             ylim=[0, 21], yticks=np.arange(0, 21, 3))

    # plot cumulative concentration-distribution
    df_cumsum = df_timeseries.div(df_timeseries.sum(axis=0), axis=1).cumsum(axis=0)
    df_cumsum.index = mdates.date2num(df_cumsum.index)
    df_cumsum.plot(linewidth=0.1, color="black", legend=False, ax=axs[1], alpha=0.3)
    locator = mdates.AutoDateLocator(minticks=4, maxticks=5)
    locator.intervald[4] = [6]
    formatter = mdates.DateFormatter('%H:%M')  # locator
    axs[1].set(ylim=[0, 1], ylabel="Cumulative scenario load [-]", xlabel="Time [HH:MM]")
    axs[1].set(xlim=[mdates.date2num(pd.to_datetime("2000-01-01 00:00:00")),
                     mdates.date2num(pd.to_datetime("2000-01-02 00:00:00"))],
               title="b) Model outputs at catchment outlet")
    axs[1].xaxis.set_major_locator(locator)
    axs[1].xaxis.set_major_formatter(formatter)
    axs[1].xaxis.set_minor_locator(AutoMinorLocator(n=3))
    axs[1].grid()

    custom_lines = [Line2D([0], [0], color="tab:purple", lw=2),
                    Line2D([0], [0], color="tab:brown", lw=2),
                    Line2D([0], [0], color="black", lw=0.5)]
    axs[1].legend(custom_lines, ["Dry-weather-flow", "Defecation pattern", "Model outputs"], loc="lower right")
    plt.savefig(plot_dir / f"cumulative_model_output.png", dpi=700)
    fig.show()
    return


def main():
    data_dir = Path(r"C:\Users\albert\Documents\SWMMpulse")
    # prepare_filelist(data_dir=data_dir)
    # evaluate_scenarios(data_dir=data_dir, plots=PLOTS[:])
    plot_scenarios(data_dir=data_dir, plots=PLOTS[:])
    # plot_cumulative_pollutograph(data_dir=data_dir)

    return


if __name__ == "__main__":
    main()
