import csv
import math
import os
import sys
from pathlib import Path

import HydroErr as hr
import hydroeval as he
import pandas as pd


# import cProfile
# from line_profiler_pycharm import profile


# @profile
def csv_import():
    df_cali = pd.read_csv(Path(os.path.dirname(__file__), sys.argv[1]), index_col='date',
                          usecols=['date', 'CH-0139 [m3/s]', 'CH-0139 (observed) [m3/s]'])
    df_vali = df_cali
    # simulations = df.iloc[:"2000-01-01",:4]
    end_cali = df_cali.index.searchsorted("2001-01-01")
    start_vali = df_vali.index.searchsorted("2001-01-01")
    end_vali = df_vali.index.searchsorted("2021-01-01")
    simulations_cali = df_cali.iloc[:end_cali, 0]
    observations_cali = df_cali.iloc[:end_cali, 1]
    simulations_vali = df_vali.iloc[start_vali:end_vali, 0]
    observations_vali = df_vali.iloc[start_vali:end_vali, 1]
    return simulations_cali, observations_cali, simulations_vali, observations_vali


# @profile
def calc_metrics(simulations_cali, observations_cali, simulations_vali, observations_vali):
    kge_np_cali, rs_cali, alpha_cali, beta_cali = \
        he.evaluator(obj_fn=he.kgenp, simulations=simulations_cali, evaluation=observations_cali)
    kge_np_vali, rs_vali, alpha_vali, beta_vali = \
        he.evaluator(obj_fn=he.kgenp, simulations=simulations_vali, evaluation=observations_vali)
    pbias_cali = he.evaluator(obj_fn=he.pbias, simulations=simulations_cali, evaluation=observations_cali)[0]
    pbias_vali = he.evaluator(obj_fn=he.pbias, simulations=simulations_vali, evaluation=observations_vali)[0]
    rmse_cali = he.evaluator(obj_fn=he.rmse, simulations=simulations_cali, evaluation=observations_cali)[0]
    rmse_vali = he.evaluator(obj_fn=he.rmse, simulations=simulations_vali, evaluation=observations_vali)[0]
    ve_cali = hr.ve(simulated_array=simulations_cali, observed_array=observations_cali)
    ve_vali = hr.ve(simulated_array=simulations_vali, observed_array=observations_vali)
    return kge_np_cali, rs_cali, alpha_cali, beta_cali, kge_np_vali, rs_vali, alpha_vali, beta_vali, pbias_cali, pbias_vali, rmse_cali, rmse_vali, ve_cali, ve_vali


# @profile
def write_to_file(rows):
    with open(Path(os.path.dirname(__file__), sys.argv[2]),
              'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(rows)


# @profile
def total_diag():
    simulations_cali, observations_cali, simulations_vali, observations_vali = csv_import()
    kge_np_cali, rs_cali, alpha_cali, beta_cali, kge_np_vali, rs_vali, alpha_vali, beta_vali, pbias_cali, pbias_vali, rmse_cali, rmse_vali, ve_cali, ve_vali = calc_metrics(
        simulations_cali, observations_cali, simulations_vali, observations_vali)

    fields = ['Run', 'KGE_NP', 'PBIAS', 'RMSE', 'VE', 'KGE_NP_Cost', 'PBIAS_Cost', 'rs', 'alpha', 'beta']
    row_cali = ["HYDROGRAPH_CALIBRATION", kge_np_cali[0], pbias_cali, rmse_cali, ve_cali, math.fabs((kge_np_cali - 1)),
                (math.fabs(pbias_cali)), rs_cali[0], alpha_cali[0], beta_cali[0]]
    row_vali = ["HYDROGRAPH_VALIDATION", kge_np_vali[0], pbias_vali, rmse_vali, ve_vali, math.fabs((kge_np_vali - 1)),
                (math.fabs(pbias_vali)), rs_vali[0], alpha_vali[0], beta_vali[0]]
    write_to_file([fields, row_cali, row_vali])


if __name__ == '__main__':
    # pr = cProfile.Profile()
    # pr.enable()

    # nse_cali = he.evaluator(he.nse, simulations_cali, observations_cali)
    # nse_vali = he.evaluator(he.nse, simulations_vali, observations_vali)
    # kge_np_cali, rs_cali, alpha_cali, beta_cali = he.kgenp(simulations=simulations_cali.to_numpy(),
    #                                                        evaluation=observations_cali.to_numpy())

    total_diag()

    # pr.disable()
    # after your program ends
    # pr.print_stats(sort="calls")
