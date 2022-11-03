"""
Tools to generate Raven .rv* files needed to run Raven models.
"""

import logging
import sys
from pathlib import Path

import pandas
import yaml

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stdout)
frm = logging.Formatter("{levelname}: {message} ({filename}/{funcName}/{lineno})",
                        style="{")
handler.setFormatter(frm)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.debug('Trying to read config.yaml file')
with open("raven_tools/config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
model_dir = config['ModelDir']
model_type = config['ModelName']
data_dir = config['DataDir']
project_dir = config['ProjectDir']
catchment = config['Catchment']
model_sub_dir = config['ModelSubDir']
supported_models = [
    "GR4J",
    "HYMOD",
    "HMETS",
    "HBV",
    "MOHYSE"
]

default_params = {
    "params": {
        "GR4J": {
            "GR4J_X1": "0.5",
            "GR4J_X2": "0.5",
            "GR4J_X3": "0.5",
            "GR4J_X4": "0.5",
            "Cemaneige_X1": "0.5",
            "Cemaneige_X2": "0.5",
            "Melt_Factor": "0.5",
            "Airsnow_Coeff": "0.5"
        },
        "HYMOD": {
            "HYMOD_Param_01": "0.5",
            "HYMOD_Param_02": "0.5",
            "HYMOD_Param_03": "0.5",
            "HYMOD_Param_04": "0.5",
            "HYMOD_Param_05": "0.5",
            "HYMOD_Param_06": "0.5",
            "HYMOD_Param_07": "0.5",
            "HYMOD_Param_08": "0.5",
            "HYMOD_Param_09": "0.5"
        },
        "HMETS": {
            "HMETS_Param_01": "0.5",
            "HMETS_Param_02": "0.5",
            "HMETS_Param_03": "0.5",
            "HMETS_Param_04": "0.5",
            "HMETS_Param_05": "0.5",
            "HMETS_Param_06": "0.5",
            "HMETS_Param_07": "0.5",
            "HMETS_Param_08": "0.5",
            "HMETS_Param_09": "0.5",
            "HMETS_Param_10": "0.5",
            "HMETS_Param_11": "0.5",
            "HMETS_Param_12": "0.5",
            "HMETS_Param_13": "0.5",
            "HMETS_Param_14": "0.5",
            "HMETS_Param_15": "0.5",
            "HMETS_Param_16": "0.5",
            "HMETS_Param_17": "0.5",
            "HMETS_Param_18": "0.5",
            "HMETS_Param_19": "0.5",
            "HMETS_Param_20": "0.5",
            "HMETS_Param_21": "0.5"
        },
        "HBV": {
            "HBV_Param_01": "0.5",
            "HBV_Param_02": "0.5",
            "HBV_Param_03": "0.5",
            "HBV_Param_04": "0.5",
            "HBV_Param_05": "0.5",
            "HBV_Param_06": "0.5",
            "HBV_Param_07": "0.5",
            "HBV_Param_08": "0.5",
            "HBV_Param_09": "0.5",
            "HBV_Param_10": "0.5",
            "HBV_Param_11": "0.5",
            "HBV_Param_12": "0.5",
            "HBV_Param_13": "0.5",
            "HBV_Param_14": "0.5",
            "HBV_Param_15": "0.5",
            "HBV_Param_16": "0.5",
            "HBV_Param_17": "0.5",
            "HBV_Param_18": "0.5",
            "HBV_Param_19": "0.5",
            "HBV_Param_20": "0.5",
            "HBV_Param_21": "0.5"
        },
        "MOHYSE": {
        }
    },
    "names": {
        "GR4J": {
            "GR4J_X1": "GR4J_x1",
            "GR4J_X2": "GR4J_x2",
            "GR4J_X3": "GR4J_x3",
            "GR4J_X4": "GR4J_x4",
            "Cemaneige_X1": "Cemaneige_x1",
            "Cemaneige_X2": "Cemaneige_x2",
            "Melt_Factor": "Melt_Factor",
            "Airsnow_Coeff": "Airsnow_Coeff"
        },
        "HYMOD": {
            "HYMOD_Param_01": "Res_Constant",
            "HYMOD_Param_02": "C_max",
            "HYMOD_Param_03": "T_s",
            "HYMOD_Param_04": "K_s",
            "HYMOD_Param_05": "Melt_Factor",
            "HYMOD_Param_06": "DD_Melt_Temp",
            "HYMOD_Param_07": "B_exp",
            "HYMOD_Param_08": "PET_Correction",
            "HYMOD_Param_09": "ALPHA"
        },
        "HMETS": {
            "HMETS_Param_01": "GAMMA_SHAPE",
            "HMETS_Param_02": "GAMMA_SCALE",
            "HMETS_Param_03": "GAMMA_SHAPE2",
            "HMETS_Param_04": "GAMMA_SCALE2",
            "HMETS_Param_05": "MIN_MELT_FACTOR",
            "HMETS_Param_06": "0.5",
            "HMETS_Param_07": "DD_MELT_TEMP",
            "HMETS_Param_08": "DD_AGGRADATION",
            "HMETS_Param_09": "SNOW_SWI_MIN",
            "HMETS_Param_10": "0.5",
            "HMETS_Param_11": "SWI_REDUCT_COEFF",
            "HMETS_Param_12": "DD_REFREEZE_TEMP",
            "HMETS_Param_13": "REFREEZE_FACTOR",
            "HMETS_Param_14": "REFREEZE_EXP",
            "HMETS_Param_15": "PET_CORRECTION_TOPSOIL",
            "HMETS_Param_16": "HMETS_RUNOFF_COEFF",
            "HMETS_Param_17": "PERC_COEFF_TOPSOIL",
            "HMETS_Param_18": "BASEFLOW_COEFF_TOPSOIL",
            "HMETS_Param_19": "BASEFLOW_COEFF_PHREATIC",
            "HMETS_Param_20": "0.5",
            "HMETS_Param_21": "0.5"
        },
        "HBV": {
            "HBV_Param_01": "RAINSNOW_TEMP",
            "HBV_Param_02": "MELT_FACTOR",
            "HBV_Param_03": "REFREEZE_FACTOR",
            "HBV_Param_04": "SNOW_SWI",
            "HBV_Param_05": "POROSITY",
            "HBV_Param_06": "FIELD_CAPACITY",
            "HBV_Param_07": "HBV_BETA",
            "HBV_Param_08": "MAX_PERC_RATE_FAST_RES",
            "HBV_Param_09": "BASEFLOW_COEFF_FAST_RES",
            "HBV_Param_10": "BASEFLOW_COEFF_SLOW_RES",
            "HBV_Param_11": "TIME_CONC_MAX_BAS",
            "HBV_Param_12": "PRECIP_LAPSE_PCALT",
            "HBV_Param_13": "ADIABATIC_LAPSE_TCALT",
            "HBV_Param_14": "SAT_WILT",
            "HBV_Param_15": "1_PLUS_ALPHA",
            "HBV_Param_16": "MAX_CAP_RISE_RATE",
            "HBV_Param_17": "THICKNESS_TOPSOIL",
            "HBV_Param_18": "HBV_MELT_FOR_CORR",
            "HBV_Param_19": "GLAC_STORAGE_COEFF",
            "HBV_Param_20": "RAIN_CORRECTION_RFCF",
            "HBV_Param_21": "SNOW_CORRECTION_SFCF"
        },
        "MOHYSE": {
        }
    }
}

header_line = "#########################################################################"
file_type = ":FileType          rvt ASCII Raven 3.5"
author = ":WrittenBy         Peter Zweifel"
creation_date = ":CreationDate      April 2022"
description = [
    "#",
    "# Emulation of GR4J simulation of Broye",
    "#------------------------------------------------------------------------"]
header = [header_line, file_type, author, creation_date, *description]

subsection_header_line: str = "#-----------------------------------------------------------------"
newline: str = "\n"


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


def write_rvh(model_dir=config['ModelDir'], model_type=config['ModelName'], data_dir=config['DataDir'],
              project_dir=config['ProjectDir'], catchment=config['Catchment'], model_sub_dir=config['ModelSubDir'],
              attribute_csv="CH-0057_attributes.csv"):
    assert model_type in supported_models, f"Model type {model_type} not in list of supported models."
    csv_file = pandas.read_csv(Path(project_dir, data_dir, "Hydromap Attributes", attribute_csv), sep=",",
                               skiprows=[8],
                               index_col='attribute_names')
    file_name: str = f"{catchment}_{model_type}.rvh"
    file_path: Path = Path(project_dir, model_dir, catchment, model_type, model_sub_dir, file_name)
    if model_type == "GR4J":
        try:
            sub_basins: list[str] = [
                ":SubBasins",
                "  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,       GAUGED",
                "  :Units     ,          none,          none,   none,          km,         none",
                "            1,        Broye_Payerne,            -1,   NONE,       _AUTO,     1",
                ":EndSubBasins"
            ]

            hrus: list[str] = [
                ":HRUs",
                "  :Attributes,  AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID,LAND_USE_CLASS, VEG_CLASS, SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT",
                "  :Units     ,   km2,         m,      deg,       deg,     none,          none,      none,         none,            none,          none,   deg,    deg",
                f"            1, {csv_file.loc['area_ch1903plus']['values']},     {csv_file.loc['a0401_eu_dem_v11_e40n20crp_chv1_0']['values']},   {csv_file.loc['lab_y']['values']},     {csv_file.loc['lab_x']['values']},        1,        LU_ALL,   VEG_ALL,    DEFAULT_P,          [NONE],        [NONE],   7.0,  217.0",
                ":EndHRUs"
            ]
        except:
            print("An error occured")

    elif model_type == "HYMOD":
        pass
    elif model_type == "HMETS":
        pass
    elif model_type == "HBV":
        pass
    elif model_type == "MOHYSE":
        pass
    else:
        print("Model not supported, please choose GR4J, HYMOD, HMETS, HBV or MOHYSE.")
    with open(file_path, 'w') as ff:
        ff.writelines(f"{line}{newline}" for line in header)
        ff.write(newline)
        ff.writelines(f"{line}{newline}" for line in subsection_header("Sub-Basins"))
        ff.writelines(f"{line}{newline}" for line in sub_basins)
        ff.write(newline)
        ff.writelines(f"{line}{newline}" for line in subsection_header("HRUs"))
        ff.writelines(f"{line}{newline}" for line in hrus)


def write_rvt(start_year: int, end_year: int, model_dir=model_dir, model_type=model_type,
              project_dir=project_dir, catchment=catchment, model_sub_dir=model_sub_dir):
    """Write to Raven *.rvt file.

    :param model_sub_dir:
    :param catchment:
    :param project_dir:
    :param int start_year: Start year of forcings data files
    :param int end_year: End year of focings data files
    :param Path model_dir: Root directory of *.rvX files
    :param Path model_type: Name of the .rvt file to be written

    """
    file_name: str = f"{catchment}_{model_type}.rvt"
    file_path: Path = Path(project_dir, model_dir, catchment, model_type, model_sub_dir, file_name)

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
    with open(file_path, 'w') as ff:
        ff.writelines(header)
        ff.write(f"# meteorological forcings\n")
        for f in forcing_block(start_year, end_year).values():
            for t in f:
                ff.write(f"{t}\n")

        ff.writelines(gauge)
        ff.writelines(flow_observation)


def generate_template(model_type, csv_file, params=default_params, param_or_name="params"):
    """Generates template text which can be written to .rvp file

    :param model_type:
    :param csv_file:
    :param params:
    :param param_or_name:
    :return:
    """

    assert model_type in supported_models, f"model_type expected GR4J, HYMOD, HMETS, HBV or MOHYSE, got {model_type} instead "
    logger.debug("Starting if-tree for model_type...")
    if model_type == "GR4J":
        logger.debug(f"Model type is {model_type}.")
        logger.debug(f"Building rvp dictionary...")
        rvp = {
            "Soil Classes":
                [
                    ":SoilClasses",
                    "   :Attributes",
                    "   :Units",
                    "       SOIL_PROD",
                    "       SOIL_ROUT",
                    "       SOIL_TEMP",
                    "       SOIL_GW",
                    ":EndSoilClasses"
                ],
            "Soil Profiles":
                [
                    "#     name,#horizons,{soiltype,thickness}x{#horizons}",
                    "#     GR4J_X1 is thickness of first layer (SOIL_PROD), here 0.529",
                    ":SoilProfiles",
                    f"   DEFAULT_P, 4, SOIL_PROD, {params[param_or_name][model_type]['GR4J_X1']}, SOIL_ROUT, 0.300, SOIL_TEMP, 1.000, SOIL_GW, 1.000,",
                    ":EndSoilProfiles"
                ],
            "Vegetation Classes":
                [
                    ":VegetationClasses",
                    "   :Attributes, MAX_HT, MAX_LAI, MAX_LEAF_COND",
                    "   :Units, m, none, mm_per_s",
                    "   VEG_ALL, 0.0, 0.0, 0.0",
                    ":EndVegetationClasses"
                ],
            "Land Use Classes":
                [
                    ":LandUseClasses",
                    "   :Attributes, IMPERM, FOREST_COV",
                    "   :Units, frac, frac",
                    f"   LU_ALL, {int(csv_file.loc['a0425_clc18_5_aurbv1_0']['values']) / 100}, {int(csv_file.loc['a0418_clc18_5_afrtv1_0']['values']) / 100}",
                    ":EndLandUseClasses"
                ],
            "Global Parameters":
                [
                    ":GlobalParameter RAINSNOW_TEMP       0.0",
                    ":GlobalParameter RAINSNOW_DELTA      1.0",
                    f":GlobalParameter AIRSNOW_COEFF     {params[param_or_name][model_type]['Airsnow_Coeff']} # [1/d] = 1.0 - CEMANEIGE_X2 = 1.0 - x6",
                    ":GlobalParameter AVG_ANNUAL_SNOW    16.9 # [mm]  =       CEMANEIGE_X1 =       x5",
                    "#:GlobalParameter PRECIP_LAPSE     0.0004 I assume not necessary for gridded data",
                    "#:GlobalParameter ADIABATIC_LAPSE  0.0065 not necessary for gridded data"
                ],
            "Soil Parameters":
                [
                    ":SoilParameterList",
                    "   :Parameters, POROSITY, GR4J_X3, GR4J_X2",
                    "   :Units, none, mm, mm / d",
                    f"   [DEFAULT], 1.0, {params[param_or_name][model_type]['GR4J_X3']}, {params[param_or_name][model_type]['GR4J_X2']}",
                    ":EndSoilParameterList"
                ],
            "Land Use Parameters":
                [
                    ":LandUseParameterList",
                    "   :Parameters, GR4J_X4, MELT_FACTOR",
                    "   :Units, d, mm / d / C",
                    f"   [DEFAULT], {params[param_or_name][model_type]['GR4J_X4']}, {params[param_or_name][model_type]['Melt_Factor']}",
                    ":EndLandUseParameterList"
                ]
        }
        logger.debug(f"rvp dictionary built.")
        logger.debug(f"Returning rvp dictionary...")
        return rvp
    elif model_type == "HYMOD":
        logger.debug(f"Model type is {model_type}.")
        logger.debug(f"Building rvp dictionary...")
        rvp = {
            "Soil Classes":
                [
                    ":SoilClasses",
                    "   :Attributes",
                    "   :Units",
                    "       TOPSOIL",
                    "       GWSOIL",
                    ":EndSoilClasses"
                ],
            "Soil Profiles":
                [
                    "#     name,#horizons,{soiltype,thickness}x{#horizons}",
                    "# ",
                    ":SoilProfiles",
                    "   LAKE, 0,",
                    "   ROCK, 0,",
                    "   # DEFAULT_P,      2, TOPSOIL,  HYMOD_PARA_2, GWSOIL, 10.0",
                    f"   DEFAULT_P, 2, TOPSOIL, {params[param_or_name][model_type]['HYMOD_Param_02']}, GWSOIL, 10.0",
                    ":EndSoilProfiles"
                ],
            "Land Use Classes":
                [
                    ":LandUseClasses",
                    "   :Attributes, IMPERM, FOREST_COV,",
                    "   :Units, frac, frac,",
                    "       LU_ALL, 0.0, 1.0",
                    ":EndLandUseClasses"
                ],
            "Vegetation Classes":
                [
                    ":VegetationClasses,",
                    "   :Attributes, MAX_HT, MAX_LAI, MAX_LEAF_COND,",
                    "   :Units, m, none, mm_per_s,",
                    "       VEG_ALL, 0.0, 0.0, 0.0",
                    ":EndVegetationClasses"
                ],
            "Global Parameters":
                [
                    f":GlobalParameter RAINSNOW_TEMP {params[param_or_name][model_type]['HYMOD_Param_03']}",
                    "   #:GlobalParameter      RAINSNOW_TEMP    HYMOD_PARA_3"
                ],
            "Soil Parameters":
                [
                    ":SoilParameterList",
                    "   :Parameters, POROSITY, PET_CORRECTION, BASEFLOW_COEFF,",
                    "   :Units, -, -, 1 / d,",
                    "       # TOPSOIL,            1.0 ,    HYMOD_PARA_8,               0.0,",
                    "       #  GWSOIL,            1.0 ,             1.0,   HYMOD_PARA_4=Ks,",
                    f"       TOPSOIL, 1.0, {params[param_or_name][model_type]['HYMOD_Param_08']}, 0.0,",
                    f"       GWSOIL, 1.0, 1.0, {params[param_or_name][model_type]['HYMOD_Param_04']},",
                    ":EndSoilParameterList"
                ],
            "Land Use Parameters":
                [
                    ":LandUseParameterList",
                    "   :Parameters, MELT_FACTOR, DD_MELT_TEMP, PDM_B,",
                    "   :Units, mm / d / K, degC, -,",
                    "       # [DEFAULT],    HYMOD_PARA_5,    HYMOD_PARA_6,  HYMOD_PARA_7=Bexp,",
                    f"       [DEFAULT], {params[param_or_name][model_type]['HYMOD_Param_05']}, {params[param_or_name][model_type]['HYMOD_Param_06']}, {params[param_or_name][model_type]['HYMOD_Param_07']},",
                    ":EndLandUseParameterList"
                ]
        }
        logger.debug(f"rvp dictionary built.")
        logger.debug(f"Returning rvp dictionary...")
        return rvp
    elif model_type == "HBV":
        logger.debug(f"Model type is {model_type}.")
        logger.debug(f"Building rvp dictionary...")
        rvp = {
            "Soil Classes":
                [
                    ":SoilClasses",
                    "   :Attributes,",
                    "   :Units,",
                    "       TOPSOIL,      1.0,    0.0,       0",
                    "       SLOW_RES,     1.0,    0.0,       0",
                    "       FAST_RES,     1.0,    0.0,       0",
                    ":EndSoilClasses"
                ],
            "Soil Profiles":
                [
                    "#     name,#horizons,{soiltype,thickness}x{#horizons}",
                    "# ",
                    f"   DEFAULT_P,      3,    TOPSOIL,            {params[param_or_name][model_type]['HBV_Param_17']},   FAST_RES,    100.0, SLOW_RES,    100.0",
                    ":EndSoilProfiles"
                ],
            "Vegetation Classes":
                [
                    ":VegetationClasses",
                    "   :Attributes, MAX_HT, MAX_LAI, MAX_LEAF_COND,",
                    "   :Units, m, none, mm_per_s,",
                    "       VEG_ALL, 0.0, 0.0, 0.0",
                    ":EndVegetationClasses"
                ],
            "Vegetation Parameters":
                [
                    ":VegetationParameterList",
                    "   :Parameters,  MAX_CAPACITY, MAX_SNOW_CAPACITY,  TFRAIN,  TFSNOW,",
                    "   :Units,                 mm,                mm,    frac,    frac,",
                    "       VEG_ALL,             10000,             10000,    0.88,    0.88,",
                    ":EndVegetationParameterList"
                ],
            "Land Use Classes":
                [
                    ":LandUseClasses",
                    "   :Attributes, IMPERM, FOREST_COV,",
                    "   :Units, frac, frac,",
                    "       LU_ALL, 0.0, 1.0",
                    ":EndLandUseClasses"
                ],
            "Global Parameters":
                [
                    f":GlobalParameter RAINSNOW_TEMP       {params[param_or_name][model_type]['HBV_Param_01']}",
                    ":GlobalParameter RAINSNOW_DELTA      1.0 #constant",
                    f"#:GlobalParameter PRECIP_LAPSE     {params[param_or_name][model_type]['HBV_Param_12']} # I assume not necessary for gridded data, HBV_PARA_12=PCALT",
                    f"#:GlobalParameter ADIABATIC_LAPSE  {params[param_or_name][model_type]['HBV_Param_13']} # not necessary for gridded data, HBV_PARA_13=TCALT",
                    f":GlobalParameter SNOW_SWI  {params[param_or_name][model_type]['HBV_Param_04']} #HBV_PARA_04"
                ],
            "Land Use Parameters":
                [
                    ":LandUseParameterList",
                    "  :Parameters,   MELT_FACTOR, MIN_MELT_FACTOR,   HBV_MELT_FOR_CORR, REFREEZE_FACTOR, HBV_MELT_ASP_CORR",
                    "  :Units     ,        mm/d/K,          mm/d/K,                none,          mm/d/K,              none",
                    "  #              HBV_PARA_02,        CONSTANT,         HBV_PARA_18,     HBV_PARA_03,          CONSTANT",
                    f"    [DEFAULT],  {params[param_or_name][model_type]['HBV_Param_02']},             2.2,        {params[param_or_name][model_type]['HBV_Param_18']},    {params[param_or_name][model_type]['HBV_Param_03']},              0.48",
                    ":EndLandUseParameterList",
                    "",
                    ":LandUseParameterList",
                    " :Parameters, HBV_MELT_GLACIER_CORR,   HBV_GLACIER_KMIN, GLAC_STORAGE_COEFF, HBV_GLACIER_AG",
                    " :Units     ,                  none,                1/d,                1/d,           1/mm",
                    "   #                       CONSTANT,           CONSTANT,        HBV_PARA_19,       CONSTANT,",
                    f"   [DEFAULT],                  1.64,               0.05,       {params[param_or_name][model_type]['HBV_Param_19']},           0.05",
                    ":EndLandUseParameterList"

                ],
            "Soil Parameters":
                [
                    ":SoilParameterList",
                    "  :Parameters,                POROSITY,FIELD_CAPACITY,     SAT_WILT,     HBV_BETA, MAX_CAP_RISE_RATE,  MAX_PERC_RATE,  BASEFLOW_COEFF,            BASEFLOW_N",
                    "  :Units     ,                    none,          none,         none,         none,              mm/d,           mm/d,             1/d,                  none",
                    "  #                        HBV_PARA_05,   HBV_PARA_06,  HBV_PARA_14,  HBV_PARA_07,       HBV_PARA_16,       CONSTANT,        CONSTANT,              CONSTANT,",
                    f"    [DEFAULT],            {params[param_or_name][model_type]['HBV_Param_05']},  {params[param_or_name][model_type]['HBV_Param_06']}, {params[param_or_name][model_type]['HBV_Param_14']}, {params[param_or_name][model_type]['HBV_Param_07']},      {params[param_or_name][model_type]['HBV_Param_16']},            0.0,             0.0,                   0.0",
                    "  #                                                        CONSTANT,                                     HBV_PARA_08,     HBV_PARA_09, 1+HBV_PARA_15=1+ALPHA,",
                    f"     FAST_RES,                _DEFAULT,      _DEFAULT,          0.0,     _DEFAULT,          _DEFAULT,   {params[param_or_name][model_type]['HBV_Param_08']},    {params[param_or_name][model_type]['HBV_Param_09']},              1.877607",
                    "  #                                                        CONSTANT,                                                      HBV_PARA_10,              CONSTANT,",
                    f"     SLOW_RES,                _DEFAULT,      _DEFAULT,          0.0,     _DEFAULT,          _DEFAULT,       _DEFAULT,    {params[param_or_name][model_type]['HBV_Param_10']},                   1.0",
                    ":EndSoilParameterList"
                ]
        }
        logger.debug(f"rvp dictionary built.")
        logger.debug(f"Returning rvp dictionary...")
        return rvp
    elif model_type == "HMETS":
        logger.debug(f"Model type is {model_type}.")
        logger.debug(f"Building rvp dictionary...")
        rvp = {
            "Soil Classes":
                [
                    ":SoilClasses",
                    "   :Attributes,",
                    "   :Units,",
                    "       TOPSOIL,",
                    "       PHREATIC,",
                    ":EndSoilClasses"
                ],
            "Soil Profiles":
                [
                    "#     name,#horizons,{soiltype,thickness}x{#horizons}",
                    "   LAKE, 0",
                    "   ROCK, 0",
                    "   # DEFAULT_P, 2, TOPSOIL,          x(20)/1000, PHREATIC,         x(21)/1000,",
                    f"  DEFAULT_P, 2, TOPSOIL,     {params[param_or_name][model_type]['HMETS_Param_20b']}, PHREATIC,     {params[param_or_name][model_type]['HMETS_Param_21b']},",
                    ":EndSoilProfiles"
                ],
            "Vegetation Classes":
                [
                    ":VegetationClasses",
                    "   :Attributes, MAX_HT, MAX_LAI, MAX_LEAF_COND,",
                    "   :Units, m, none, mm_per_s,",
                    "       FOREST,             4,             5,             5,",
                    ":EndVegetationClasses"
                ],
            "Land Use Classes":
                [
                    ":LandUseClasses",
                    "   :Attributes, IMPERM, FOREST_COV,",
                    "   :Units, frac, frac,",
                    "       FOREST, 0.0, 1.0",
                    ":EndLandUseClasses"
                ],
            "Global Parameters":
                [
                    f":GlobalParameter  SNOW_SWI_MIN {params[param_or_name][model_type]['HMETS_Param_09']} # x(9)",
                    f":GlobalParameter  SNOW_SWI_MAX {params[param_or_name][model_type]['HMETS_Param_09b']} # x(9)+x(10)",
                    f":GlobalParameter  SWI_REDUCT_COEFF {params[param_or_name][model_type]['HMETS_Param_11']} # x(11)",
                    ":GlobalParameter SNOW_SWI 0.05 #not sure why/if needed"
                ],
            "Vegetation Parameters":
                [
                    ":VegetationParameterList",
                    "   :Parameters,  RAIN_ICEPT_PCT,  SNOW_ICEPT_PCT,",
                    "   :Units,               -,               -,",
                    "       [DEFAULT],             0.0,             0.0,",
                    ":EndVegetationParameterList"
                ],
            "Land Use Parameters":
                [
                    ":LandUseParameterList",
                    "   :Parameters, MIN_MELT_FACTOR, MAX_MELT_FACTOR,    DD_MELT_TEMP,  DD_AGGRADATION, REFREEZE_FACTOR,    REFREEZE_EXP, DD_REFREEZE_TEMP, HMETS_RUNOFF_COEFF,",
                    "   :Units,          mm/d/C,          mm/d/C,               C,            1/mm,          mm/d/C,               -,                C,                  -,",
                    f"      [DEFAULT],  {params[param_or_name][model_type]['HMETS_Param_05']}, {params[param_or_name][model_type]['HMETS_Param_05b']},  {params[param_or_name][model_type]['HMETS_Param_07']},  {params[param_or_name][model_type]['HMETS_Param_08']},  {params[param_or_name][model_type]['HMETS_Param_13']},  {params[param_or_name][model_type]['HMETS_Param_14']},   {params[param_or_name][model_type]['HMETS_Param_12']},     {params[param_or_name][model_type]['HMETS_Param_16']},",
                    "#      x(5),       x(5)+x(6),            x(7),            x(8),           x(13),           x(14),            x(12),              x(16),",
                    ":EndLandUseParameterList",
                    "",
                    ":LandUseParameterList",
                    "   :Parameters,     GAMMA_SHAPE,     GAMMA_SCALE,    GAMMA_SHAPE2,    GAMMA_SCALE2,",
                    "   :Units,               -,             1/d,               -,             1/d,",
                    f"      [DEFAULT],  {params[param_or_name][model_type]['HMETS_Param_01']},  {params[param_or_name][model_type]['HMETS_Param_02']},  {params[param_or_name][model_type]['HMETS_Param_03']},  {params[param_or_name][model_type]['HMETS_Param_04']},",
                    "#      x(1),            x(2),            x(3),            x(4),",
                    ":EndLandUseParameterList"
                ],
            "Soil Parameters":
                [
                    ":SoilParameterList",
                    "   :Parameters,        POROSITY,      PERC_COEFF,  PET_CORRECTION, BASEFLOW_COEFF",
                    "   :Units,               -,             1/d,               -,            1/d",
                    f"      TOPSOIL,             1.0,  {params[param_or_name][model_type]['HMETS_Param_17']},  {params[param_or_name][model_type]['HMETS_Param_15']}, {params[param_or_name][model_type]['HMETS_Param_18']}",
                    f"      PHREATIC,             1.0,             0.0,             0.0, {params[param_or_name][model_type]['HMETS_Param_19']}",
                    "#      TOPSOIL,             1.0,           x(17),           x(15),          x(18)",
                    "#      PHREATIC,             1.0,             0.0,             0.0,          x(19)",
                    ":EndSoilParameterList"
                ]
        }
        logger.debug(f"rvp dictionary built.")
        logger.debug(f"Returning rvp dictionary...")
        return rvp
    elif model_type == "MOHYSE":
        logger.debug(f"Model type is {model_type}.")
        logger.debug(f"Building rvp dictionary...")
        rvp = {
            "Soil Classes":
                [
                    ":SoilClasses",
                    "   :Attributes,",
                    "   :Units,",
                    "       TOPSOIL",
                    "       GWSOIL",
                    ":EndSoilClasses"
                ],
            "Soil Profiles":
                [
                    "#  name,#horizons,{soiltype,thickness}x{#horizons}",
                    "# ",
                    "   LAKE, 0",
                    "   ROCK, 0",
                    "#  DEFAULT_P,      2, TOPSOIL, MOHYSE_PARA_5, GWSOIL, 10.0",
                    f"  DEFAULT_P,      2, TOPSOIL,     {params[param_or_name][model_type]['MOHYSE_Param_05']}, GWSOIL, 10.0",
                    ":EndSoilProfiles"
                ],
            "Vegetation Classes":
                [
                    ":VegetationClasses",
                    "   :Attributes, MAX_HT, MAX_LAI, MAX_LEAF_COND,",
                    "   :Units, m, none, mm_per_s,",
                    "       VEG_ALL, 0.0, 0.0, 0.0",
                    ":EndVegetationClasses"
                ],
            "Vegetation Parameters":
                [
                    ":VegetationParameterList",
                    "   :Parameters,    SAI_HT_RATIO,  RAIN_ICEPT_PCT,  SNOW_ICEPT_PCT,",
                    "   :Units,               -,               -,               -, ",
                    "       [DEFAULT],             0.0,             0.0,             0.0,   ",
                    ":EndVegetationParameterList"
                ],
            "Land Use Classes":
                [
                    ":LandUseClasses",
                    "   :Attributes, IMPERM, FOREST_COV,",
                    "   :Units, frac, frac,",
                    "       LU_ALL, 0.0, 1.0",
                    ":EndLandUseClasses"
                ],
            "Global Parameters":
                [
                    "#:GlobalParameter      RAINSNOW_TEMP              -2.0",
                    ":GlobalParameter       TOC_MULTIPLIER              1.0",
                    f"# :GlobalParameter     MOHYSE_PET_COEFF  {params[param_or_name][model_type]['MOHYSE_Param_01']}",
                    ":GlobalParameter       MOHYSE_PET_COEFF         1.0000"
                ],
            "Land Use Parameters":
                [
                    ":LandUseParameterList",
                    "   :Parameters,     MELT_FACTOR,       AET_COEFF, FOREST_SPARSENESS, DD_MELT_TEMP,",
                    "   :Units,          mm/d/K,            mm/d,                 -,         degC,",
                    "#      [DEFAULT],   MOHYSE_PARA_3,   MOHYSE_PARA_2,               0.0,MOHYSE_PARA_4, ",
                    f"      [DEFAULT],          {params[param_or_name][model_type]['MOHYSE_Param_03']},          {params[param_or_name][model_type]['MOHYSE_Param_02']},               0.0,       {params[param_or_name][model_type]['MOHYSE_Param_04']},",
                    ":EndLandUseParameterList"
                ],
            "Soil Parameters":
                [
                    ":SoilParameterList",
                    "   :Parameters,        POROSITY,  PET_CORRECTION,        HBV_BETA,  BASEFLOW_COEFF,      PERC_COEFF, ",
                    "   :Units,               -,               -,               -,             1/d,             1/d, ",
                    "#      TOPSOIL,            1.0 ,             1.0,             1.0,   MOHYSE_PARA_7,   MOHYSE_PARA_6,",
                    "#      GWSOIL,            1.0 ,             1.0,             1.0,   MOHYSE_PARA_8,             0.0,",
                    f"      TOPSOIL,            1.0 ,             1.0,             1.0,          {params[param_or_name][model_type]['MOHYSE_Param_07']},          {params[param_or_name][model_type]['MOHYSE_Param_06']},",
                    f"      GWSOIL,            1.0 ,             1.0,             1.0,          {params[param_or_name][model_type]['MOHYSE_Param_08']},             0.0,",
                    ":EndSoilParameterList"
                ]
        }
        logger.debug(f"rvp dictionary built.")
        logger.debug(f"Returning rvp dictionary...")
        return rvp
    else:
        logger.debug(f"None of if/elif statements evaluated to True, since model type is: {model_type}")
        print(f"Model type {model_type} currently not supported.")


# TODO: is 'melt_factor' actually a parameter of GR4J? In the tutorial file, 'annual snow' is called 'CEMANEIGE_X1'...


def write_rvp(model_dir=model_dir,
              model_type=model_type,
              data_dir=data_dir,
              project_dir=project_dir,
              catchment=catchment,
              model_sub_dir=model_sub_dir,
              params=default_params,
              template_type="Raven",
              attribute_csv_name="CH-0057_attributes.csv",
              attribute_csv_dir="Hydromap Attributes"):
    """Writes .rvp file, either as an Ostrich or Raven template.

    :param model_dir:
    :param model_type:
    :param data_dir:
    :param project_dir:
    :param catchment:
    :param model_sub_dir:
    :param params:
    :param raven_template:
    :param ostrich_template:
    """

    logger.debug("Arrived in function write_rvp().")
    assert model_type in supported_models, f"Got model type: {model_type}, which is not supported, check variable \"" \
                                           f"supported_models."
    logger.debug(f"Trying to read catchment attribute CSV file {Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name)}...")
    csv_file = pandas.read_csv(Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name), sep=",",
                               skiprows=[8],
                               index_col='attribute_names')
    logger.debug("Attribute catchment attribute CSV file read.")
    file_name: str = f"{catchment}_{model_type}.rvp"
    logger.debug(f".rvp filename set to {file_name}.")
    file_path: Path = Path(project_dir, model_dir, catchment, model_type, model_sub_dir, file_name)
    logger.debug(f".rvp file path set to {file_path}.")
    template_sections = {}
    logger.debug("Empty dict template_sections created.")
    logger.debug("Entering if-tree for template type...")
    if template_type == "Raven":
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Trying to generate .rvp template sections with function generate_template()...")
        template_sections = generate_template(model_type=model_type, csv_file=csv_file, params=params,
                                              param_or_name="params")
        logger.debug("Wrote .rvp template sections generated by generate_template() to dict template_sections")
    if template_type == "Ostrich":
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Adding .tpl suffix to file_path for Ostrich template file...")
        file_path = Path((str(file_path) + ".tpl"))
        logger.debug(f"New file_path: {file_path}")
        logger.debug(f"Trying to generate .rvp.tpl template sections with function generate_template()...")
        template_sections = generate_template(model_type=model_type, csv_file=csv_file, params=params,
                                              param_or_name="names")
        logger.debug("Wrote .rvp.tpl template sections generated by generate_template() to dict template_sections")
    logger.debug(f"Trying to write to file {file_path}")
    with open(file_path, 'w') as ff:
        # print("Header:")
        ff.writelines(f"{line}{newline}" for line in header)
        logger.debug("Header lines written.")
        ff.write(newline)
        logger.debug("Entering template_sections for-loop...")
        for section in template_sections:
            # print(f"Section Header: {section}")
            ff.writelines(f"{line}\n" for line in subsection_header(section))
            logger.debug("Subsection header written.")
            # for subsection in raven_sections[section]:
            # print(f"Subsection: {template_sections[section]}")
            ff.writelines(f"{lin}\n" for lin in template_sections[section])
            logger.debug("Template section written")
            ff.write(newline)
        logger.debug("template_sections for-loop finished.")
        # print(template_sections)
        logger.debug("Function write_rvp() finished.")


def subsection_header(title: str) -> list[str]:
    """Generates subsection header for a .rvX file from given title.

    :param str title: Title to be used in subsection header.
    :return subsection_head: List of header line as strings.
    :rtype subsection_head: list[str]
    """
    logger.debug("Arrived in function subsection_header().")
    subsection_head: list[str] = [
        subsection_header_line,
        f"# ----{title}--------------------------------------------",
        subsection_header_line
    ]
    logger.debug("Subsection header written.")
    logger.debug("Returning subsection header...")
    logger.debug("Leaving function subsection_header() afterwards...")
    return subsection_head


if __name__ == '__main__':
    # model_dir = Path(config['ModelDir'])
    # model_type = Path("raven_broye_gr4j.rvt")
    start_year = 1981
    end_year = 2000
    write_rvt(1981, 1982)
