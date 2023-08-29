"""
Tools to processing-process Raven output files.
"""
import glob
import re
import sys
from pathlib import Path

import matplotlib.animation as animation
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import raven_tools.config.variables as var

# model_path = Path("/media/mainman/Work/RAVEN/results/incoming")
model_path = ("/media/mainman/Work/RAVEN/results/incoming")
model_types = ["GR4J", "MOHYSE", "HYMOD", "HMETS", "HBV"]

csv_dict = {}
perf_metrics = ['Run', 'KGE_NP', 'PBIAS', 'RMSE', 'VE']

catchments_by_id = [key for key in var.catchments]
catchments_by_id = ["CH-0053"]
# model_name = "HYMOD"
file_list = []


# other = pd.read_csv(Path(model_path, "CH-0053", model_name, f"OstNonDomSolutions0.txt"), sep="\s+", skiprows=1)


def pareto(model_type: str):
    """Analyse Non-Dominated pareto solutions from a multi-threaded Ostrich run.

    Args:
        model_type: str
            Model type, e.g. 'GR4J'
    """
    cols = ["KGE_NP", "PBIAS", "RMSE", "VE"]
    df = pd.DataFrame(columns=cols)
    for sol in range(0, 8):
        ot = pd.read_csv(Path(model_path, "CH-0053", model_type, f"OstNonDomSolutions{sol}.txt"), sep="\s+",
                         skiprows=1)
        df = pd.concat([df, ot])
    df = df[cols]
    axs = pd.plotting.scatter_matrix(df, alpha=0.2, figsize=(6, 6), diagonal="kde")
    axs[0, 2].xaxis.set_label_text("test")
    axs[0, 1].xaxis.get_label()

    plt.savefig(f"/home/mainman/Documents/Studium/UniBe/Master's Thesis/data/figures/pareto.png", dpi=300)


def pareto2(model_name: str):
    """Analyse performance metrics from a multi-threaded Ostrich run.

    Args:
        model_name:
    """
    cols = ["KGE_NP", "PBIAS", "RMSE", "VE"]
    df = pd.DataFrame(columns=cols)
    for sol in range(0, 8):
        ot = pd.read_csv(Path(model_path, "CH-0053", model_name, f"OstModel{sol}.txt"), sep="\s+")
        df = pd.concat([df, ot])
    df = df[cols]
    df.index[df.duplicated(subset=cols)]
    pd.plotting.scatter_matrix(df, alpha=0.2, figsize=(6, 6), diagonal="kde")
    plt.savefig(f"/home/mainman/Documents/Studium/UniBe/Master's Thesis/data/figures/pareto2.png", dpi=300)


def ani_gr4J():
    """Shows a plot with performance metrics for a running, multi-threaded Ostrich run, which auto-updates

    Start the Ostrich run, then use this function the get auto-updated plot of the performance metrics.
    """
    model_name = "GR4J"
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ctm = "CH-0053"

    def animate(i):
        cols = ["Run", "KGE_NP", "PBIAS", "RMSE", "VE"]
        df = pd.DataFrame(columns=cols)
        for sol in range(0, 8):
            ot = pd.read_csv(Path(model_path, ctm, model_name, f"OstModel{sol}.txt"), sep="\s+", skiprows=0)
            df = pd.concat([df, ot])
        df = df[cols].sort_values(by="KGE_NP")
        df.index = np.arange(1, len(df) + 1)
        run = df.index.values
        kge_np = df['KGE_NP'].values
        # fig, axes = plt.plot()
        ax1.clear()
        ax1.scatter(run, kge_np)
        fig.supxlabel('Ostrich Iteration Steps')
        fig.suptitle(f"{model_name} - {ctm}")
        # axes = df.plot(kind='scatter', x='Ctm', y='KGE_NP', color='orange')

    ani = animation.FuncAnimation(fig, animate, interval=1000)
    plt.show()


# pareto2("HYMOD")


def perf_diag_across_ctms(model_type: str):
    """Compare performance metrics for multiple catchments

    Args:
        model_type: str
            Model type, e.g. 'GR4J'.
    """
    df_names = {'cali_df': 'cali', 'vali_df': 'vali'}
    df_list = []
    df_out = pd.DataFrame()
    for df, df_str in df_names.items():
        df_gen = df_diag_generator(df, model_type, cali_or_vali=df_str)
        df_out = pd.concat([df_out, df_gen])
        df_list.append(df_out)
    new_cols_clean = [clean_text(f"CH-", tx) for tx in df_out.columns.tolist()]
    df_out = df_out.rename(columns=dict(zip(df_out.columns.tolist(), new_cols_clean)))
    metrics_df = pd.DataFrame(columns=['Run', 'KGE_NP', 'PBIAS', 'RMSE', 'VE'])
    metrics_list = perf_metrics[1:5]
    # fig = plt.figure(figsize=(10, 10), dpi=150, layout='constrained')
    # ax = plt.axes()
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(13, 13), dpi=150, layout='constrained')
    marker_list = ['o', 'x']
    # for nn, ax in enumerate(axes.flat):
    # for df in df_list:
    # pd.plotting.parallel_coordinates(df_list, 'type', color=plt.cm.tab20(np.arange(len(df))), ls='',
    #                                  marker='x', axvlines=False, ax=ax)
    xlabels = ["KGE_NP (-)", "PBIAS (%)", "RMSE ()", "VE ()"]
    axtitles = ["Non-Parametric KGE (ranked)", "Percent Bias", "Root Mean Squared Error",
                "Volumetric Efficiency"]
    for nn, ax in enumerate(axes.flat):
        # if nn <= 6:
        # ax.set_xlim(left=xlim_left, right=xlim_right)
        # if nn == 0:
        #     ax.set_ylim(-1, 1)
        # ax.invert_yaxis()
        # else:
        # ax.hlines(bounds[nn], 0, 500, colors='red', linestyles='dotted')
        # pass
        # df_sol[f"{df_sol.columns[nn + 1]}"].plot(ax=ax, ylabel=xlabels[nn])
        pd.plotting.parallel_coordinates(df_out[df_out.index == metrics_list[nn]], 'type',
                                         color=plt.cm.Set1(np.arange(2)),
                                         ls='',
                                         marker='o',
                                         axvlines=False, ax=ax)
        ax.set_title(axtitles[nn])
        ax.set_ylabel(ylabel=xlabels[nn])
        ax.margins(x=0.2, y=0.6)
    # else:
    #     pass
    fig.supxlabel('Ostrich Iteration Steps')
    fig.suptitle(f"Model: {model_type}\nCalibration with 5000 Runs")
    plt.savefig(
        f"{model_path}/figures/perf_comp_{model_type}_cali.png")
    # fig.supxlabel('Ostrich Iteration Steps')
    # fig.suptitle(f"{m} - {c}")
    # plt.show()
    # for pm in perf_metrics:
    #     cali_list.plot.box(column=[pm])
    # fig, axes = plt.subplots(ncols=8, dpi=300)
    # fig.suptitle(model_name)
    # # plt.show()

    # plt.savefig(f"/home/mainman/Documents/Studium/UniBe/Master's Thesis/data/figures/{model_name}.png")


def df_diag_generator(df_name, model_type, cali_or_vali: str):
    """Reads performance metrics from CSV file into DataFrame

    Args:
        df_name:
            Name for DataFrame to use
        model_type: str
            Model type, e.g. 'GR4J'.
        cali_or_vali: str
            Calibration or Validation - 'cali' or 'vali'

    Returns:
        df: pd.DataFrame
            Dataframe with performance metrics.

    """
    df = pd.DataFrame(columns=perf_metrics)
    df.name = df_name
    if cali_or_vali == "cali":
        cali_vali = 'HYDROGRAPH_CALIBRATION'
        cali_vali_col = 'cali'
    if cali_or_vali == "vali":
        cali_vali = 'HYDROGRAPH_VALIDATION'
        cali_vali_col = 'vali'
    for c in catchments_by_id:
        file_path = Path(model_path, c, model_type, "processor_0", "model_best",
                         f"{c}_{model_type}_Diagnostics.csv")
        csv_pd = pd.read_csv(open(file_path), sep=",")
        d = csv_pd.loc[csv_pd['Run'] == cali_vali, perf_metrics]
        # vali = csv_pd.loc[csv_pd['Run'] == 'HYDROGRAPH_VALIDATION', perf_metrics]
        d.reset_index(inplace=True)
        # vali.reset_index(inplace=True)
        d.loc[0, perf_metrics[0]] = c
        # vali.loc[0, perf_metrics[0]] = c
        df = pd.concat([df, d])
    df = df[perf_metrics]
    df.set_index(perf_metrics[0], inplace=True)
    df = df.transpose()
    df['type'] = cali_vali_col
    # df.reset_index(inplace=True)
    return df


def df_sol_generator(ctm_ch_id: str, model_type: str):
    """Loads solutions from OstModel files and returns them in a single DataFrame.

    Args:
        ctm_ch_id: str
            Catchment id
        model_type: str
            Model type, e.g. 'GR4J'

    Returns:
        df: pd.DataFrame
            Combined solutions.

    """

    cols = ["Run", "KGE_NP_CALI", "PBIAS_CALI", "RMSE_CALI", "VE_CALI"]
    #    cols = ["Run", "KGE_NP", "PBIAS", "RMSE", "VE"]
    df = pd.DataFrame(columns=cols)
    #    for sol in range(0, 100):
    #        ot = pd.read_csv(Path(model_path, ctm, model_type, f"OstModel{sol}.txt"), sep="\s+", skiprows=0)
    #        df = pd.concat([df, ot])
    df = ostmodel_loader(model_path, ctm_ch_id, model_type)
    df = df[cols].sort_values(by=cols[1])
    df.index = np.arange(1, len(df) + 1)
    run = df.index.values
    # kge_np = df['KGE_NP'].values
    return df


def ostmodel_loader(model_path: str, ctm_ch_id: str, model_type: str):
    """Loads solutions from multi-threaded Ostrich OstModel files into single DataFrame.

    Args:
        model_path: str
            Path to model files
        ctm_ch_id: str
            Catchment id
        model_type: str
            Model type, e.g. 'GR4'

    Returns:
        df: pd.DataFrame
            DataFrame with OstModel outputs

    """
    df = pd.concat([pd.read_csv(f, sep="\s+", skiprows=0) for f in
                    glob.glob(model_path + "/" + ctm_ch_id + "/" + model_type + "/OstModel*.txt")])
    return df


def df_factors_generator(ctm_ch_id: str, model_type: str):
    """Loads performance metrics from OstModel files and returns them in a single DataFrame.

    Args:
        ctm_ch_id: str
            Catchment id.
        model_type: str
            Model type, e.g. 'GR4J'.

    Returns:
        df: pd.DataFrame
            Performance metrics combined in a single DataFrame.

    """
    cols = ["Run", "obj.function", "KGE_NP_CALI", "PBIAS_CALI", "RMSE_CALI", "VE_CALI", "KGE_NP_VALI", "PBIAS_VALI",
            "RMSE_VALI", "VE_VALI"]
    df = pd.DataFrame(columns=cols)
    new_cols = []
    df = ostmodel_loader(model_path, ctm_ch_id, model_type)
    for col in df.columns.tolist():
        if col in cols:
            pass
        else:
            new_cols.append(col)
    new_cols_clean = [clean_text(f"{model_type}_", tx) for tx in new_cols]
    df = df.sort_values(by="KGE_NP_CALI")
    df = df[new_cols].rename(columns=dict(zip(new_cols, new_cols_clean)))
    df.index = np.arange(1, len(df) + 1)
    return df


def clean_text(rgx_patt, text):
    """Cleans text

    Args:
        rgx_patt:
        text:

    Returns:

    """
    new_text = text
    new_text = re.sub(rgx_patt, '', new_text)
    return new_text


def perf_metrics_chart(model_type: str):
    """Creates plots with performance metrics evolving over the Ostrich iteration runs and saves it to file.

    Args:
        model_type: str
            Model type, e.g. 'GR4J'.
    """
    with plt.ioff():
        xlim_left = 0
        xlim_right = 5001
        # for m in model_types:
        # for c in catchments_by_id:
        for c in catchments_by_id:
            df_names = {'cali_df': 'cali', 'vali_df': 'vali'}
            ylabels = ["KGE_NP (-)", "PBIAS (%)", "RMSE ()", "VE ()"]
            axtitles = ["Non-Parametric KGE (ranked)", "Percent Bias", "Root Mean Squared Error",
                        "Volumetric Efficiency"]
            df_list = []

            # for df, df_str in df_names.items():
            #     df_out = pd.DataFrame()
            #     df_gen = df_diag_generator(df, model_type, type=df_str)
            #     df_out = pd.concat([df_out, df_gen])
            #     df_list.append(df_out)
            df_sol = df_sol_generator(c, model_type)
            fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharex=True, dpi=150, layout='constrained')

            for nn, ax in enumerate(axes.flat):
                if nn <= 6:
                    ax.set_xlim(left=xlim_left, right=xlim_right)
                    if nn == 0:
                        ax.set_ylim(-1, 1)
                        # ax.invert_yaxis()
                    else:
                        # ax.hlines(bounds[nn], 0, 500, colors='red', linestyles='dotted')
                        pass
                    try:
                        df_sol[f"{df_sol.columns[nn + 1]}"].plot(ax=ax, ylabel=ylabels[nn])
                    except TypeError:
                        pass
                    ax.set_title(axtitles[nn])
                else:
                    pass
            fig.supxlabel('Ostrich Iteration Steps')
            fig.suptitle(f"Model: {model_type}\nCatchment: {c}\nCalibration with {xlim_right - 1} Runs")
            plt.savefig(
                f"{model_path}/figures/perf_{model_type}_{c}_{xlim_right}_cali.png")


def bounds_reader(model_type: str):
    """Reads upper and lower bounds of default parameters values from config file.

    Args:
        model_type: str
            Model type, e.g. 'GR4J'

    Returns:
        upper: list
            List with upper limits.
        lower: list
            List with lower limits.

    """
    upper = list(var.default_params[f"{model_type}"]["upper"].values())
    lower = list(var.default_params[f"{model_type}"]["lower"].values())
    return upper, lower


def model_factors_chart(model_type: str):
    """Creates plots for each model type that show performance metrics and save to files.

    Args:
        model_type: str
            Model type, e.g. 'GR4J'
    """
    # with plt.ioff():
    xlim_left = 0
    xlim_right = 400
    bounds = [
        [],
        [0.01, 2.5],
        [-15, 10],
        [10, 700],
        [0.5, 7],
        [0.5, 30],
        [0, 1]
    ]
    nn_max = 0
    # for m in model_types:
    for c in catchments_by_id:
        #    c = "CH-0139"
        file_path = Path(model_path, c, model_type, "dds_status.out")
        factors = pd.read_csv(file_path, sep='\t')
        factors.set_index("STEP", inplace=True)
        df_names = {'cali_df': 'cali', 'vali_df': 'vali'}
        xlabels = ["KGE_NP (-)", "PBIAS (%)", "RMSE (?)", "VE (?)"]
        axtitles = ["Non-Parametric KGE", "Percent Bias", "Root Mean Squared Error", "Volumetric Efficiency"]
        df_list = []

        # for df, df_str in df_names.items():
        #     df_out = pd.DataFrame()
        #     df_gen = df_diag_generator(df, model_type, type=df_str)
        #     df_out = pd.concat([df_out, df_gen])
        #     df_list.append(df_out)
        df_sol = df_factors_generator(c, model_type)
        bnd_upper, bnd_lower = bounds_reader(model_type=model_type)
        if model_type == "GR4J":
            nrows = 2
            ncols = 3
            nn_max = 6
        elif model_type == "HBV":
            nrows = 7
            ncols = 3
            nn_max = 21
        elif model_type == "HMETS":
            nrows = 7
            ncols = 3
            nn_max = 21
        elif model_type == "HYMOD":
            nrows = 2
            ncols = 3
            nn_max = 9
        elif model_type == "MOHYSE":
            nrows = 4
            ncols = 3
            nn_max = 9
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(13, 13), sharex=True, dpi=150, layout='constrained')
        for nn, ax in enumerate(axes.flat):
            if nn <= nn_max:
                # else:

                #     pass
                # factors[f"{factors.columns[nn]}"].plot(ax=ax)
                # ax.set_title(factors.columns[nn])
                try:
                    df_sol[f"{df_sol.columns[nn]}"].plot(ax=ax)
                except IndexError:
                    pass
                except TypeError:
                    pass
                ax.set_xlim(left=xlim_left, right=xlim_right)
                ax.set_ylim(
                    [float(bnd_lower[nn]) - 0.2 * float(ax.get_ylim()[1] - ax.get_ylim()[0]),
                     float(bnd_upper[nn]) + 0.2 * float(ax.get_ylim()[1] - ax.get_ylim()[0])])
                ax.hlines(y=[float(bnd_lower[nn]), float(bnd_upper[nn])], xmin=0, xmax=xlim_right, colors='red',
                          linestyles='dotted')
                # ax.hlines(y=bnd_upper[nn], xmin=0, xmax=xlim_right, colors='red', linestyles='dotted')
                try:
                    ax.set_title(df_sol.columns[nn])
                except IndexError:
                    pass
            else:
                pass
        fig.supxlabel('Ostrich Iteration Steps')
        fig.suptitle(f"Model: {model_type}\nCatchment: {c}\nCalibration with {xlim_right - 1} Runs")
        #   plt.show()
        plt.savefig(
            f"{model_path}/figures/factors_{model_type}_{c}_{xlim_right}.png")
        plt.close()


def hydrograph():
    """Plot the hydrograph generated in a RAVEN run.

    """
    hydro = pd.read_csv(Path(model_path, "CH-0053/HBV/processor_0/model/output/CH-0053_HBV_Hydrographs.csv", sep=","))
    hydro['date'] = pd.to_datetime(hydro['date'], format='%Y-%m-%d')
    hydro.set_index('date', inplace=True)
    hydro.drop(columns=['time', 'hour'], inplace=True)
    hydro.plot()

    watershedStorage = pd.read_csv(
        Path(model_path, "CH-0053/HBV/processor_0/model/output/CH-0053_HBV_WatershedStorage.csv", sep=','))
    watershedStorage['date'] = pd.to_datetime(watershedStorage['date'], format='%Y-%m-%d')
    watershedStorage.set_index('date', inplace=True)
    watershedStorage.drop(columns=['time [d]', 'hour'], inplace=True)
    watershedStorage.plot()

    end_cali = watershedStorage.index.searchsorted("2001-01-01")
    start_vali = watershedStorage.index.searchsorted("2001-01-01")
    end_vali = watershedStorage.index.searchsorted("2021-01-01")
    simulations_cali = watershedStorage.iloc[:end_cali]
    observations_cali = watershedStorage.iloc[:end_cali]
    simulations_vali = watershedStorage.iloc[start_vali:end_vali]
    observations_vali = watershedStorage.iloc[start_vali:end_vali]

    end_test = watershedStorage.index.searchsorted('1981-12-31')
    df_test = watershedStorage.iloc[:end_test]
    df_test.plot()


# model_factors_chart("MOHYSE")
model_type = sys.argv[1]
perf_diag_across_ctms(model_type=model_type)
model_factors_chart(model_type=model_type)
perf_metrics_chart(model_type=model_type)
if __name__ == '__main__':
    pass
