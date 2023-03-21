import csv
import os
import sys
from pathlib import Path

import HydroErr as hr
import hydroeval as he
import pandas as pd

if __name__ == '__main__':
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

    # nse_cali = he.evaluator(he.nse, simulations_cali, observations_cali)
    # nse_vali = he.evaluator(he.nse, simulations_vali, observations_vali)

    kge_np_cali = he.evaluator(obj_fn=he.kgenp, simulations=simulations_cali, evaluation=observations_cali)[0][0]
    kge_np_vali = he.evaluator(obj_fn=he.kgenp, simulations=simulations_vali, evaluation=observations_vali)[0][0]
    pbias_cali = he.evaluator(obj_fn=he.pbias, simulations=simulations_cali, evaluation=observations_cali)[0]
    pbias_vali = he.evaluator(obj_fn=he.pbias, simulations=simulations_vali, evaluation=observations_vali)[0]
    rmse_cali = he.evaluator(obj_fn=he.rmse, simulations=simulations_cali, evaluation=observations_cali)[0]
    rmse_vali = he.evaluator(obj_fn=he.rmse, simulations=simulations_vali, evaluation=observations_vali)[0]
    ve_cali = hr.ve(simulated_array=simulations_cali, observed_array=observations_cali)
    ve_vali = hr.ve(simulated_array=simulations_vali, observed_array=observations_vali)

    fields = ['Run', 'KGE_NP', 'PBIAS', 'RMSE', 'VE']
    row_cali = ["HYDROGRAPH_CALIBRATION", kge_np_cali, pbias_cali, rmse_cali, ve_cali]
    row_vali = ["HYDROGRAPH_VALIDATION", kge_np_vali, pbias_vali, rmse_vali, ve_vali]

    with open(Path(os.path.dirname(__file__), sys.argv[2]),
              'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerow(row_cali)
        csvwriter.writerow(row_vali)
