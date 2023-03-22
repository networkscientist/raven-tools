"""
Tools to processing-process Raven output files.
"""
import os
import sys
from pathlib import Path

import HydroErr as he
import pandas as pd
from matplotlib import pyplot as plt

if __name__ == '__main__':
    df_cali = pd.read_csv(
        Path(
            "/media/mainman/Work/RAVEN/models/CH-0058/HYMOD/processor_0/model_best/CH-0058_HYMOD_Hydrographs.csv"),
        usecols=[1, 3, 4, 5], index_col='date')
    df_cali = df_cali.iloc[1:, ]
    df_cali.plot()
    plt.show()
    df_cali.set_index('date', inplace=True)
    df_vali = pd.read_csv(Path(os.path.dirname(__file__), sys.argv[1]), index_col=1, parse_dates=True)
    df_cali = pd.read_csv(
        Path("/media/mainman/Work/RAVEN/models/CH-0058/GR4J/processor_0/model/output/CH-0058_GR4J_Hydrographs.csv"),
        index_col=1)
    df_vali = pd.read_csv(Path(""), index_col=1)

    end_cali = df_cali.index.searchsorted("2000-12-31")
    start_vali = df_vali.index.searchsorted("2001-01-01")
    end_vali = df_vali.index.searchsorted("2020-12-31")
    simulations_cali = df_cali.iloc[:end_cali, 3]
    observations_cali = df_cali.iloc[:end_cali, 4]
    simulations_vali = df_vali.iloc[start_vali:end_vali, 3]
    observations_vali = df_vali.iloc[start_vali:end_vali, 4]

    he.ve(simulated_array=simulations_cali, observed_array=observations_cali)
