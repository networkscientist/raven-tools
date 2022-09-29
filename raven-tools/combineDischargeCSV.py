import glob
from datetime import datetime
from pathlib import Path

import pandas as pd

home_path = Path.home()
raven_path: str = "/media/mainman/Data/RAVEN"
data_path: str = f"{raven_path}/data/"
forcings_path: str = f"{raven_path}/data/forcings/"
result_discharge_Path: str = f"{raven_path}/data/Discharge/"
original_discharge_csv_path: str = result_discharge_Path
discharge_file_name: str = "BroPay_Q_2034_daily.rvt"

column_names = {
    "rre150h0": "PRECIP",
    "tre200h0": "TEMP_AVE",
    "tre200hx": "TEMP_MAX",
    "tre200hn": "TEMP_MIN",
}

time_column_names_asc = {
    "YYYY": "year",
    "MM": "month",
    "DD": "day",
    "HH": "hour",
}
# Date range (can be overriden in function call)
start_date = '2000-01-01'
end_date = '2000-12-31'


def write_rvt_pet():
    """Script to generate rvt file from monthly PET .csv files from IDAWEB

    :return:
    """
    # Read in PET data from CSV file
    df_pet: pd.DataFrame = pd.read_csv(forcings_path + "order_103168_PAY_ets150m0_1_data.txt", sep=";", usecols=[1, 2])
    # Convert time column to datetime format and sort according to it
    df_pet['time'] = pd.to_datetime(df_pet['time'], format="%Y%m")
    df_pet = df_pet.sort_values(by='time', ascending=True).dropna()
    df_pet = subset_dataframe_time(df_pet,start_date,end_date)
    # Export to RAVEN .rvt file

    with open(f"{forcings_path}GaugePAY_pet.rvt", 'w') as f:
        print("Writing PET data...")
        f.write(f":Data\n2000-01-01\t0:00:00\t30\t{len(df_pet)}\n")
        f.write(":Parameters\tPET_MONTH_AVE\n")
        f.write(":Units\tmm/d\n")
        dfAsString = df_pet.to_string(justify="right", header=False, index=False, columns={'ets150m0'})
        f.write(dfAsString)
        f.write("\n:EndData")
        print(" ... done")


# ## Script to generate rvt file from bi-hourly q value .asc file
def write_rvt_discharge_from_asc():
    """Function to generate rvt file from bi-hourly q value .asc file

    """


# ## Script to generate rvt file from daily q value .csv files, separated in yearly files.
def write_rvt_discharge_from_hydromaps():
    """Function to generate rvt file from daily q value .csv files, separated in yearly files.

    """
    # Read in the data from CSV files and concatenate into a single DataFrame
    df_q: pd.DataFrame = pd.concat([pd.read_csv(f, sep=',') for f in glob.glob(original_discharge_csv_path + "/*.csv")],
                                   ignore_index=True)
    # Convert dt column to datetime format and sort by it
    df_q['dt'] = pd.to_datetime(df_q['dt'])
    df_q = df_q.sort_values(by='dt', ascending=True).dropna()
    # Export to RAVEN .rvt file
    with open(f"{result_discharge_Path}BroPay_Q_2034_daily.rvt", 'w') as f:
        f.write(":ObservationData\tHYDROGRAPH\t1\tm3/s\n")
        f.write("1920-01-01\t0:00:00\t1\t36160\n")
        df_as_string = df_q.to_string(justify="right", header=False, index=False, columns={'q'})
        f.write(df_as_string)
        f.write("\n:EndObservationData")


def subset_dataframe_time(dataframe: pd.DataFrame, start_date: str, end_date: str) -> pd.DataFrame:
    """Subsetting a dataframe using a time interval.

        :param pd.DataFrame dataframe: Original DataFrame to subset.
        :param str start_date: Start date (inclusive)
        :param str end_date: End date (inclusive)
        :return pd.DataFrame subset_dataframe: Subset DataFrame
    """
    # Date to string conversion
    start_date = datetime.strptime(start_date, "%Y-%m-%d")
    end_date = datetime.strptime(end_date, "%Y-%m-%d")
    # Create data interval mask. Offsetting by 1 day to include end date.
    mask = dataframe['time'].between(start_date, end_date + pd.DateOffset(days=1), inclusive="left")
    # Apply the mask
    subset_dataframe = dataframe[mask]
    return subset_dataframe


def asc_to_rvt(ofile):
    df_meteo: pd.DataFrame = pd.read_csv(original_discharge_csv_path + ofile, sep="\t")
    # Convert time column to datetime format
    df_meteo = df_meteo.rename(columns=time_column_names_asc)

    df_meteo['time'] = pd.to_datetime(df_meteo[['year','month','day','hour']], format="%Y%m%d%H")
    # Sort by date
    # df_meteo_ordered = df_meteo_ordered.sort_values(by='time', ascending=True).dropna()
    # Subset according to start and end date
    start_date = "2000-01-01"
    start_time = "0:01:00"
    end_date = "2000-12-31"
    df_meteo = subset_dataframe_time(df_meteo, start_date, end_date)
    export_to_rvt_file(start_date, start_time, df_meteo)

def export_to_rvt_file(date, time,df):
    with open(result_discharge_Path + discharge_file_name, 'w') as f:
        # print(rvt_filename)
        f.write(f":ObservationData\tHYDROGRAPH\t1\tm3/s\n{date}\t{time}\t0.083333333\t{len(df)}\n")
        # For gauged precipitation data
        df_as_string = df.to_string(justify="right", header=False, index=False,
                                                  columns=['Q.BroPay'])
        f.write(df_as_string)
        # f.write(df_as_string)
        f.write("\n:EndObservationData")


asc_to_rvt("BroPay_Q_2034_hourly.asc")


## Function to generate rvt file from daily meteo value .csv files from IDAWEB
def write_rvt_meteo(netcdf: bool = False, start_date: str = '2008-01-01', end_date: str = '2008-12-31'):
    """Write the rvt file for the forcings data.

    This function is used to generate RAVEN .rvt files for gauged temperature data and either gauged or gridded
    precipitation data for a specified date range.

    :param bool netcdf: True if used for netCDF format, else False
    :param str start_date: Start date (inclusive). Either use value set in starting block or set here.
    :param str end_date: End date (inclusive). Either use value set in starting block or set here.
    """
    rvt_filename: str = "GaugePAY_meteo"
    # Check if netCDF is used for precipitation and set file name accordingly
    if netcdf:
        rvt_filename: str = f"{forcings_path}/{rvt_filename}_nc.rvt"
        print("Writing netCDF RVT file...")
    else:
        rvt_filename: str = f"{forcings_path}/{rvt_filename}.rvt"
        print("Writing RVT file using gauged precipitation data...")
    # Read average temperature values from csv file. Only use datetime and temperature columns.
    df_meteo: pd.DataFrame = pd.read_csv(forcings_path + "order_103168_PAY_tre200h0_1_data.txt", sep=";",
                                         usecols=[1, 2])
    # Read in the other data files using a RegEx selector.
    meteo_files = Path(forcings_path).glob('*re*data.txt')
    for f in meteo_files:
        # print(f)
        # Append the data to the df_meteo DataFrame to generate a single file
        meteo_append: pd.DataFrame = pd.read_csv(f, sep=";", usecols=[1, 2])
        df_meteo = pd.merge(df_meteo, meteo_append)
    # Rename and reorder columns
    df_meteo = df_meteo.rename(columns=column_names)
    df_meteo_ordered: pd.DataFrame = df_meteo[['time', 'TEMP_AVE', 'TEMP_MIN', 'TEMP_MAX', 'PRECIP']]

    # Convert time column to datetime format
    df_meteo_ordered['time'] = pd.to_datetime(df_meteo_ordered['time'], format="%Y%m%d%H")
    # Sort by date
    df_meteo_ordered = df_meteo_ordered.sort_values(by='time', ascending=True).dropna()
    # Subset according to start and end date
    df_meteo_ordered = subset_dataframe_time(df_meteo_ordered, start_date, end_date)

    # Export to RVT file
    with open(rvt_filename, 'w') as f:
        # print(rvt_filename)
        f.write(f":MultiData\n{start_date}\t0:00:00\t0.041666667\t{len(df_meteo_ordered)}\n")

        # For netCDF precipitation data
        if netcdf:
            f.write(":Parameters\tTEMP_AVE\tTEMP_MIN\tTEMP_MAX\n")
            f.write(":Units\tC\tC\tC\n")
            df_as_string = df_meteo_ordered.to_string(justify="right", header=False, index=False,
                                                      columns=['TEMP_AVE', 'TEMP_MIN', 'TEMP_MAX'])
        # For gauged precipitation data
        else:
            f.write(":Parameters\tPRECIP\tTEMP_AVE\tTEMP_MIN\tTEMP_MAX\n")
            f.write(":Units\tmm/d\tC\tC\tC\n")
            df_as_string = df_meteo_ordered.to_string(justify="right", header=False, index=False,
                                                      columns=['PRECIP', 'TEMP_AVE', 'TEMP_MIN',
                                                               'TEMP_MAX'])
        f.write(df_as_string)
        # f.write(df_as_string)
        f.write("\n:EndMultiData")
    print("... done")


# class RvtFile:
#     def __init__(self, start_date="2008-01-01", end_date="2008-12-31", rvt_filename="GaugePAY"):
#         self.StartDate = start_date
#         self.EndDate = end_date
#         self.RvtPath = Path(f"{forcings_path}/{rvt_filename}")
#     def set_start_date(self, new_date: str):
#         self.StartDate = new_date
#         return True


write_rvt_meteo(netcdf=False, start_date='2008-01-01', end_date='2008-12-31')
write_rvt_meteo(netcdf=True, start_date='2008-01-01', end_date='2008-12-31')
write_rvt_discharge_from_hydromaps()
write_rvt_pet()