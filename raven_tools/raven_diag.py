import sys
import hydroeval as he
import pandas as pd
import datetime as dt
import csv

if __name__ == '__main__':
    df = pd.read_csv(sys.argv[1], index_col=1)
    df = pd.read_csv("/home/sirian/Applications/Hydrology/RAVEN/testmodels/GR4J/processor_0/model/output/raven_broye_gr4j_Hydrographs.csv", index_col=1)
    # simulations = df.iloc[:"2000-01-01",:4]
    end = df.index.searchsorted("2000-01-01")
    simulations = df.iloc[:end,3]
    observations = df.iloc[:end,4]

    nse = he.evaluator(he.nse, simulations, observations)
    kge_orig, r_orig, alpha_orig, beta_orig = he.evaluator(he.kge, simulations, observations)
    kge_prime, r_prime, alpha_prime, beta_prime = he.evaluator(he.kgeprime, simulations, observations)
    kge_np, r_np, alpha_np, beta_np = he.evaluator(he.kgenp, simulations, observations)

    fields = ['Run', 'NSE', 'KGE_ORIG', 'KGE_PRIME', 'KGE_NP']
    row = ["HYDROGRAPH_CALIBRATION", nse[0],kge_orig[0],kge_prime[0],kge_np[0]]

    with open(sys.argv[2],'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerow(row)