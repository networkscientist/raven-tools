import glob
from pathlib import Path
from datetime import datetime
import time
import pandas as pd

writePath = "/home/sirian/Applications/Hydrology/RAVEN/data/"
meteo_path = "/home/sirian/Applications/Hydrology/RAVEN/data/forcings/"
resultPath = "/home/sirian/Applications/Hydrology/RAVEN/data/Discharge/result.csv"
original_csv_path = "/home/sirian/Applications/Hydrology/RAVEN/data/Discharge/original_from_hydromaps"
column_names = {
    "tre200hx" : "TEMP_MAX",
    "tre200hn" : "TEMP_MIN",
    "tre200h0" : "TEMP_AVE",
    "rre150h0" : "PRECIP"
}


## Script to generate rvt file from daily q value .csv files, separated in yearly files.
df_q = pd.concat([pd.read_csv(f, sep=',') for f in glob.glob(original_csv_path + "/*.csv")], ignore_index=True)
df_q['dt'] = pd.to_datetime(df_q['dt'])
df_q = df_q.sort_values(by='dt', ascending=True).dropna()
with open((writePath + "Discharge/discharge.rvt"), 'w') as f:
    f.write(":ObservationData\tHYDROGRAPH\t1\tm3/s\n")
    f.write("1920-01-01\t0:00:00\t1\t36160\n")
    dfAsString = df_q.to_string(justify="right", header=False, index=False, columns={'q'})
    f.write(dfAsString)
    f.write("\n:EndObservationData")


def subset_dataframe_time(dataframe, start_date, end_date):
    start_date = datetime.strptime(start_date,"%Y-%m-%d")
    end_date = datetime.strptime(end_date,"%Y-%m-%d")
    # mask = (dataframe['time'] >= start_date) & (dataframe['time']<= end_date)
    mask = dataframe['time'].between(start_date,end_date+pd.DateOffset(days=1), inclusive="left")
    dataframe = dataframe[mask]
    return dataframe


## Function to generate rvt file from daily meteo value .csv files from IDAWEB
def write_rvt_meteo(pet: bool = False, start_date: object = '2008-01-01', end_date: object = '2008-12-31') -> object:
    """
    Write the rvt file for the forcings data.
    :param start_date:
    :param end_date:
    :param pet: True if used for netCDF format, else False
    """
    rvt_filename = "GaugePAY_meteo"
    if pet:
        filename = f"{writePath}forcings/{rvt_filename}_nc.rvt"
    else:
        filename = f"{writePath}{rvt_filename}.rvt"
    df_meteo = pd.read_csv(meteo_path + "order_103168_PAY_tre200h0_1_data.txt", sep=";", usecols=[1, 2])
    meteo_files = Path(meteo_path).glob('*re*data.txt')
    for f in meteo_files:
        print(f)
        meteo_append = pd.read_csv(f, sep=";", usecols=[1, 2])
        df_meteo = pd.merge(df_meteo, meteo_append)
    df_meteo = df_meteo.rename(columns=column_names)
    df_meteo = df_meteo[['time','TEMP_AVE','TEMP_MIN','TEMP_MAX','PRECIP']]
    df_meteo['time'] = pd.to_datetime(df_meteo['time'], format="%Y%m%d%H")
    df_meteo = df_meteo.sort_values(by='time', ascending=True).dropna()
    df_meteo = subset_dataframe_time(df_meteo, start_date, end_date)
    with open(filename, 'w') as f:
        print(filename)
        f.write(f":MultiData\n{start_date}\t0:00:00\t0.041666667\t{len(df_meteo)}\n")
        f.write(":Parameters\tTEMP_AVE\tTEMP_MIN\tTEMP_MAX\n")
        f.write(":Units\tC\tC\tC\n")
        dfAsString = df_meteo.to_string(justify="right", header=False, index=False, columns={'TEMP_AVE','TEMP_MIN','TEMP_MAX'})
        f.write(dfAsString)
        f.write("\n:EndMultiData")


write_rvt_meteo(pet=True, start_date='2008-01-01', end_date='2008-12-31')

## Script to generate rvt file from monthly PET .csv files from IDAWEB
df_pet = pd.read_csv(meteo_path + "order_103168_PAY_ets150m0_1_data.txt", sep=";", usecols=[1, 2])

df_pet['time'] = pd.to_datetime(df_pet['time'], format="%Y%m")
df_pet = df_pet.sort_values(by='time', ascending=True).dropna()
with open((writePath + "forcings/GaugePAY_pet.rvt"), 'w') as f:
    f.write(":Data\n1981-01-01\t0:00:00\t0.041666667\t498\n")
    f.write(":Parameters\tPET_MONTH_AVE\n")
    f.write(":Units\tmm/d\n")
    dfAsString = df_pet.to_string(justify="right", header=False, index=False, columns={'ets150m0'})
    f.write(dfAsString)
    f.write("\n:EndData")


def generate_discharge_rvt(orig_path, orig_filename, out_path):
    df_q = pd.read_csv(orig_path + orig_filename, sep='\t')
    df_q['time'] = df_q.YYYY.map(str) + df_q.MM.map(str)
    df_q['DD'] = df_q['DD'].apply(lambda x: '{0:0>2}'.format(x))
    return df_q


out_path = "Discharge/Broye_Near_Payerne_Qobs_daily.rvt"
orig_filename = "Discharge/BroPay_Q_2034_hourly.asc"
df_q = generate_discharge_rvt(writePath, orig_filename, out_path)


def write_rvt_file(**data_type):
    if data_type["type"] == "Multi":
        print("Multi")
    else:
        print("Something else")


write_rvt_file(type="Multi")
