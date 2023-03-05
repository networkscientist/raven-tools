import csv
import os
import sys
from pathlib import Path

import hydroeval as he
import pandas as pd

if __name__ == '__main__':
    # df_cali = pd.read_csv(
    #     "/media/mainman/Work/RAVEN/models/Dischmabach/GR4J/model/output/Dischmabach_GR4J_Hydrographs.csv",
    #     index_col=1)
    # df_vali = pd.read_csv(
    #     "/media/mainman/Work/RAVEN/models/Dischmabach/GR4J/model/output/Dischmabach_GR4J_Hydrographs.csv",
    #     index_col=1)

    df_cali = pd.read_csv(Path(os.path.dirname(__file__), sys.argv[1]), index_col=1)
    df_vali = pd.read_csv(Path(os.path.dirname(__file__), sys.argv[1]), index_col=1)
    # simulations = df.iloc[:"2000-01-01",:4]
    end_cali = df_cali.index.searchsorted("2000-12-31")
    start_vali = df_vali.index.searchsorted("2001-01-01")
    end_vali = df_vali.index.searchsorted("2020-12-31")
    simulations_cali = df_cali.iloc[:end_cali, 3]
    observations_cali = df_cali.iloc[:end_cali, 4]
    simulations_vali = df_vali.iloc[start_vali:end_vali, 3]
    observations_vali = df_vali.iloc[start_vali:end_vali, 4]

    nse_cali = he.evaluator(he.nse, simulations_cali, observations_cali)
    nse_vali = he.evaluator(he.nse, simulations_vali, observations_vali)
    # kge_orig, r_orig, alpha_orig, beta_orig = he.evaluator(he.kge, simulations, observations)
    # kge_prime, r_prime, alpha_prime, beta_prime = he.evaluator(he.kgeprime, simulations, observations)
    kge_np_cali, r_np, alpha_np, beta_np = he.evaluator(he.kgenp, simulations_cali, observations_cali)
    kge_np_vali, r_np, alpha_np, beta_np = he.evaluator(he.kgenp, simulations_vali, observations_vali)

    fields = ['Run', 'NSE', 'KGE_NP']
    row_cali = ["HYDROGRAPH_CALIBRATION", nse_cali[0], kge_np_cali[0]]
    row_vali = ["HYDROGRAPH_VALIDATION", nse_vali[0], kge_np_vali[0]]

    with open(Path(os.path.dirname(__file__), sys.argv[2]),
              'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerow(row_cali)
        csvwriter.writerow(row_vali)
