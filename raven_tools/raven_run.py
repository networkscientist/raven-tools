"""
Tools to generate Raven .rv* files needed to run Raven models.
"""

import os
from pathlib import Path

import pandas
import yaml

with open("raven_tools/config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)


def forcing_block(start: int, end: int):
    """Create Dictionary of forcing data to write in RVT file.

    This function creates a Dictionary of forcing data to be written into an RVT file. From a start and end year,
    it creates the relevant directory Paths, which are based on Swiss gridded input data. Further parameters,
    as required by RAVEN are added and finally, the forcings data block is returned as a Dictionary.

    :param int start: Start year of forcings data files
    :param int end: End year of focings data files
    :return: The forcing data block
    :rtype: Dictionary

    """

    forcing_data = {
        'Rainfall': [
            ":GriddedForcing           Rainfall",
            "    :ForcingType          RAINFALL",
            f"    :FileNameNC           data_obs/RhiresD_v2.0_swiss.lv95/merged/RhiresD_ch01h.swiss.lv95_{start}01010000_{end}12310000_clipped.nc",
            "    :VarNameNC            RhiresD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            "    :RedirectToFile       data_obs/GridWeights.txt ",
            ":EndGriddedForcing"],
        'Average Temperature': [
            ":GriddedForcing           Average Temperature",
            "    :ForcingType          TEMP_AVE",
            f"    :FileNameNC           data_obs/TabsD_v2.0_swiss.lv95/merged/TabsD_ch01r.swiss.lv95_{start}01010000_{end}12310000_clipped.nc",
            "    :VarNameNC            TabsD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            "    :RedirectToFile       data_obs/GridWeights.txt ",
            ":EndGriddedForcing"],
        'Maximum Temperature': [
            ":GriddedForcing           Maximum Temperature",
            "    :ForcingType          TEMP_MAX",
            f"    :FileNameNC           data_obs/TmaxD_v2.0_swiss.lv95/merged/TmaxD_ch01r.swiss.lv95_{start}01010000_{end}12310000_clipped.nc",
            "    :VarNameNC            TmaxD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            "    :RedirectToFile       data_obs/GridWeights.txt ",
            ":EndGriddedForcing"],
        'Minimum Temperature': [
            ":GriddedForcing           Minimum Temperature",
            "    :ForcingType          TEMP_MIN",
            f"    :FileNameNC           data_obs/TminD_v2.0_swiss.lv95/merged/TminD_ch01r.swiss.lv95_{start}01010000_{end}12310000_clipped.nc",
            "    :VarNameNC            TminD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            "    :RedirectToFile       data_obs/GridWeights.txt ",
            ":EndGriddedForcing"
        ]
    }

    forcing_header = f"# Years: {start} to {end}\n"

    return forcing_data


def write_rvh(model_dir=config['ModelDir'], model_name=config['ModelName'], data_dir=config['DataDir']):
    csv_file = pandas.read_csv(data_dir + "Hydromap Attributes/CH-0057_attributes.csv", sep=",", skiprows=[8],
                               index_col='attribute_names')
    model_name = (model_name + ".rvh")
    file_name = Path(model_name + ".rvh")
    file_type = ":FileType          rvt ASCII Raven 3.5"
    author = ":WrittenBy         Peter Zweifel"
    creation_date = ":CreationDate      April 2022"
    description = "#\n# Emulation of GR4J simulation of Broye\n#------------------------------------------------------------------------\n"
    sub_basins = ":SubBasins\n" \
                 "  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,       GAUGED\n" \
                 "  :Units     ,          none,          none,   none,          km,         none\n" \
                 "            1,        Broye_Payerne,            -1,   NONE,       _AUTO,     1\n" \
                 ":EndSubBasins\n"

    hrus = ":HRUs\n" \
           "  :Attributes,  AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID,LAND_USE_CLASS, VEG_CLASS, SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT\n" \
           "  :Units     ,   km2,         m,      deg,       deg,     none,          none,      none,         none,            none,          none,   deg,    deg\n" \
           f"            1, {csv_file.loc['area_ch1903plus']['values']},     {csv_file.loc['a0401_eu_dem_v11_e40n20crp_chv1_0']['values']},   {csv_file.loc['lab_y']['values']},     {csv_file.loc['lab_x']['values']},        1,        LU_ALL,   VEG_ALL,    DEFAULT_P,          [NONE],        [NONE],   7.0,  217.0\n" \
           ":EndHRUs"

    with open((os.path.join(model_dir, model_name)), 'w') as ff:
        ff.write(f"#########################################################################\n"
                 f"{file_type}\n"
                 f"{author}\n"
                 f"{creation_date}\n")
        ff.write(f"{description}\n")
        ff.write(f"{sub_basins}\n")
        ff.write(f"{hrus}")


def write_rvt(start_year: int, end_year: int, model_dir=config['ModelDir'], model_name=config['ModelName']):
    """Write to Raven *.rvt file.

    :param int start_year: Start year of forcings data files
    :param int end_year: End year of focings data files
    :param Path model_dir: Root directory of *.rvX files
    :param Path model_name: Name of the .rvt file to be written

    """
    model_name = Path(model_name + ".rvt")
    file_type = ":FileType          rvt ASCII Raven 3.5"
    author = ":WrittenBy         Peter Zweifel"
    creation_date = ":CreationDate      April 2022"
    description = "#\n# Emulation of GR4J simulation of Broye\n#------------------------------------------------------------------------\n"
    gauge = [
        ":Gauge PYR2034\n",
        "  :Latitude    46.835913\n",
        "  :Longitude 6.9360708\n",
        "  :Elevation  441.0\n",
        ":EndGauge\n"
    ]

    flow_observation = [
        "# observed streamflow\n",
        ":RedirectToFile data_obs/BroPay_Q_2034_daily.rvt"
    ]
    with open((os.path.join(model_dir, model_name)), 'w') as ff:
        ff.write(f"#########################################################################\n"
                 f"{file_type}\n"
                 f"{author}\n"
                 f"{creation_date}\n")
        ff.write(f"{description}\n")
        ff.write(f"# meteorological forcings\n")
        for f in forcing_block(start_year, end_year).values():
            for t in f:
                ff.write(f"{t}\n")

        ff.writelines(gauge)
        ff.writelines(flow_observation)


if __name__ == '__main__':
    model_dir = Path(config['ModelDir'])
    model_name = Path("raven_broye_gr4j.rvt")
    start_year = 1981
    end_year = 2000
    write_rvt(1981, 1982, model_dir, model_name)
