"""
Tools to processing-process Raven output files.
"""
import glob
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

model_path = Path("/media/mainman/Work/RAVEN/models")
model_types = ["GR4J", "MOHYSE", "HYMOD"]

csv_dict = {}
perf_metrics = ['KGE_NP', 'PBIAS', 'RMSE', 'VE']


def test_function(model_name: str):
    file_list = glob.glob(f"{model_path}/CH-*/{model_name}/processor_0/model_best/CH-*_{model_name}_Diagnostics.csv")
    metrics_df = pd.DataFrame(columns=['Run', 'KGE_NP', 'PBIAS', 'RMSE', 'VE'])
    cali_df = pd.DataFrame(columns=perf_metrics)
    vali_df = pd.DataFrame(columns=perf_metrics)
    for i, f in enumerate(file_list):
        csv_pd = pd.read_csv(open(f), sep=",")
        cali = csv_pd.loc[csv_pd['Run'] == 'HYDROGRAPH_CALIBRATION']
        vali = csv_pd.loc[csv_pd['Run'] == 'HYDROGRAPH_VALIDATION']
        cali_df = pd.concat([cali_df, cali])
        vali_df = pd.concat([vali_df, vali])
    # for pm in perf_metrics:
    #     cali_list.plot.box(column=[pm])
    fig, axes = plt.subplots(ncols=8, dpi=300)
    fig.suptitle(model_name)

    for i, pm in enumerate(perf_metrics):
        cali_df.boxplot(column=pm, ax=axes[2 * i])
        axes[2 * i].set_title('Cali')
        vali_df.boxplot(column=pm, ax=axes[2 * i + 1])
        axes[2 * i + 1].set_title('Vali')
    # plt.show()

    plt.savefig(f"/home/mainman/Documents/Studium/UniBe/Master's Thesis/data/figures/{model_name}.png")


for m in model_types:
    test_function(m)

if __name__ == '__main__':
    pass
