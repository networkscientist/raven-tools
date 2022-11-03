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
            "MOHYSE_Param_01": "0.5",
            "MOHYSE_Param_02": "0.5",
            "MOHYSE_Param_03": "0.5",
            "MOHYSE_Param_04": "0.5",
            "MOHYSE_Param_05": "0.5",
            "MOHYSE_Param_06": "0.5",
            "MOHYSE_Param_07": "0.5",
            "MOHYSE_Param_08": "0.5",
            "MOHYSE_Param_09": "0.5",
            "MOHYSE_Param_10": "0.5"
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
            "MOHYSE_Param_01": "MOHYSE_PET_COEFF",
            "MOHYSE_Param_02": "AET_COEFF",
            "MOHYSE_Param_03": "MELT_FACTOR",
            "MOHYSE_Param_04": "DD_MELT_TEMP",
            "MOHYSE_Param_05": "Thickness_TOPSOIL",
            "MOHYSE_Param_06": "PERC_COEFF_TOPSOIL",
            "MOHYSE_Param_07": "BASEFLOW_COEFF_TOPSOIL",
            "MOHYSE_Param_08": "BASEFLOW_COEFF_GWSOIL",
            "MOHYSE_Param_09": "GAMMA_SHAPE",
            "MOHYSE_Param_10": "TIME_CONC"

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


def generate_template(model_type, csv_file, params=default_params, param_or_name="params", generate_rvp=False,
                      generate_rvh=False, generate_rvi=False, generate_rvc=False, generate_rvt=False):
    """Generates template text which can be written to .rvp file

    :param model_type:
    :param csv_file:
    :param params:
    :param param_or_name:
    :return:
    """

    assert model_type in supported_models, f"model_type expected GR4J, HYMOD, HMETS, HBV or MOHYSE, got {model_type} instead "
    logger.debug("Starting if-tree for file_type...")
    if generate_rvp is True:
        logger.debug("generate_rvp evaluated to True...")
        logger.debug("Starting if-tree for model_type...")
        if model_type == "GR4J":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvp dictionary...")
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
            logger.debug(f"{model_type} rvp dictionary built.")
            logger.debug(f"Returning {model_type} rvp dictionary...")
            return rvp
        elif model_type == "HYMOD":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvp dictionary...")
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
            logger.debug(f"{model_type} rvp dictionary built.")
            logger.debug(f"Returning {model_type} rvp dictionary...")
            return rvp
        elif model_type == "HBV":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvp dictionary...")
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
            logger.debug(f"{model_type} rvp dictionary built.")
            logger.debug(f"Returning {model_type} rvp dictionary...")
            return rvp
        elif model_type == "HMETS":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvp dictionary...")
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
            logger.debug(f"{model_type} rvp dictionary built.")
            logger.debug(f"Returning {model_type} rvp dictionary...")
            return rvp
        elif model_type == "MOHYSE":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvp dictionary...")
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
            logger.debug(f"{model_type} rvp dictionary built.")
            logger.debug(f"Returning {model_type} rvp dictionary...")
            return rvp
        else:
            logger.debug(f"None of if/elif statements evaluated to True, since model type is: {model_type}")
            print(f"Model type {model_type} currently not supported.")

    if generate_rvh is True:
        logger.debug("generate_rvh evaluated to True...")
        logger.debug("Starting if-tree for model_type...")
        if model_type == "GR4J":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvh dictionary...")
            rvh = {
                "Subbasins":
                    [
                        ":SubBasins",
                        "  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,       GAUGED",
                        "  :Units     ,          none,          none,   none,          km,         none",
                        "            1,        Broye_Payerne,            -1,   NONE,       _AUTO,     1",
                        ":EndSubBasins"
                    ],
                "HRUs":
                    [
                        ":HRUs",
                        "  :Attributes,  AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID,LAND_USE_CLASS, VEG_CLASS, SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT",
                        "  :Units     ,   km2,         m,      deg,       deg,     none,          none,      none,         none,            none,          none,   deg,    deg",
                        f"            1, {csv_file.loc['area_ch1903plus']['values']},     {csv_file.loc['a0401_eu_dem_v11_e40n20crp_chv1_0']['values']},   {csv_file.loc['lab_y']['values']},     {csv_file.loc['lab_x']['values']},        1,        LU_ALL,   VEG_ALL,    DEFAULT_P,          [NONE],        [NONE],   {csv_file.loc['a0404_eu_dem_v11_e40n20_slp8v1_0']['values']},  {csv_file.loc['a0407_eu_dem_v11_asp8sm_maskv1_0']['values']}",
                        ":EndHRUs"
                    ]
            }
            logger.debug(f"{model_type} rvh dictionary built.")
            logger.debug(f"Returning {model_type} rvh dictionary...")
            return rvh
        elif model_type == "HYMOD":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvh dictionary...")
            rvh = {
                "Subbasins":
                    [
                        ":SubBasins",
                        "  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,       GAUGED",
                        "  :Units     ,          none,          none,   none,          km,         none",
                        "            1,        Broye_Payerne,            -1,   NONE,       _AUTO,     1",
                        ":EndSubBasins"
                    ],
                "HRUs":
                    [
                        ":HRUs",
                        "  :Attributes,  AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID,LAND_USE_CLASS, VEG_CLASS, SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT",
                        "  :Units     ,   km2,         m,      deg,       deg,     none,          none,      none,         none,            none,          none,   deg,    deg",
                        f"            1, {csv_file.loc['area_ch1903plus']['values']},     {csv_file.loc['a0401_eu_dem_v11_e40n20crp_chv1_0']['values']},   {csv_file.loc['lab_y']['values']},     {csv_file.loc['lab_x']['values']},        1,        LU_ALL,   VEG_ALL,    DEFAULT_P,          [NONE],        [NONE],   {csv_file.loc['a0404_eu_dem_v11_e40n20_slp8v1_0']['values']},  {csv_file.loc['a0407_eu_dem_v11_asp8sm_maskv1_0']['values']}",
                        ":EndHRUs"
                    ],
                "Subbasin Properties":
                    [
                        ":SubBasinProperties",
                        "#                         HYMOD_PARA_1,                  3,",
                        "   :Parameters,           RES_CONSTANT,     NUM_RESERVOIRS,",
                        "   :Units,                         1/d,                  -,",
                        f"              1,          {params[param_or_name][model_type]['HYMOD_Param_01']},                  3,",
                        ":EndSubBasinProperties"
                    ]
            }
            logger.debug(f"{model_type} rvh dictionary built.")
            logger.debug(f"Returning {model_type} rvh dictionary...")
            return rvh
        elif model_type == "HBV":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvh dictionary...")
            rvh = {
                "Subbasins":
                    [
                        ":SubBasins",
                        "  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,       GAUGED",
                        "  :Units     ,          none,          none,   none,          km,         none",
                        "            1,        Broye_Payerne,            -1,   NONE,       _AUTO,     1",
                        ":EndSubBasins"
                    ],
                "HRUs":
                    [
                        ":HRUs",
                        "  :Attributes,  AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID,LAND_USE_CLASS, VEG_CLASS, SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT",
                        "  :Units     ,   km2,         m,      deg,       deg,     none,          none,      none,         none,            none,          none,   deg,    deg",
                        f"            1, {csv_file.loc['area_ch1903plus']['values']},     {csv_file.loc['a0401_eu_dem_v11_e40n20crp_chv1_0']['values']},   {csv_file.loc['lab_y']['values']},     {csv_file.loc['lab_x']['values']},        1,        LU_ALL,   VEG_ALL,    DEFAULT_P,          [NONE],        [NONE],   {csv_file.loc['a0404_eu_dem_v11_e40n20_slp8v1_0']['values']},  {csv_file.loc['a0407_eu_dem_v11_asp8sm_maskv1_0']['values']}",
                        ":EndHRUs"
                    ],
                "Subbasin Properties":
                    [
                        ":SubBasinProperties",
                        "#                       HBV_PARA_11, DERIVED FROM HBV_PARA_11,",
                        "#                            MAXBAS,                 MAXBAS/2,"
                        "   :Parameters,           TIME_CONC,             TIME_TO_PEAK,",
                        "   :Units,                        d,                        d,",
                        f"              1,          {params[param_or_name][model_type]['HBV_Param_11']},                  {params[param_or_name][model_type]['HBV_Param_11']},",
                        ":EndSubBasinProperties"
                    ]
            }
            # TODO: Actually it is param 11b... Change!
            logger.debug(f"{model_type} rvh dictionary built.")
            logger.debug(f"Returning {model_type} rvh dictionary...")
            return rvh
        elif model_type == "HMETS":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvp dictionary...")
            rvh = {
                "Subbasins":
                    [
                        ":SubBasins",
                        "  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,       GAUGED",
                        "  :Units     ,          none,          none,   none,          km,         none",
                        "            1,        Broye_Payerne,            -1,   NONE,       _AUTO,     1",
                        ":EndSubBasins"
                    ],
                "HRUs":
                    [
                        ":HRUs",
                        "  :Attributes,  AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID,LAND_USE_CLASS, VEG_CLASS, SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT",
                        "  :Units     ,   km2,         m,      deg,       deg,     none,          none,      none,         none,            none,          none,   deg,    deg",
                        f"            1, {csv_file.loc['area_ch1903plus']['values']},     {csv_file.loc['a0401_eu_dem_v11_e40n20crp_chv1_0']['values']},   {csv_file.loc['lab_y']['values']},     {csv_file.loc['lab_x']['values']},        1,        LU_ALL,   VEG_ALL,    DEFAULT_P,          [NONE],        [NONE],   {csv_file.loc['a0404_eu_dem_v11_e40n20_slp8v1_0']['values']},  {csv_file.loc['a0407_eu_dem_v11_asp8sm_maskv1_0']['values']}",
                        ":EndHRUs"
                    ]
            }
            logger.debug(f"{model_type} rvh dictionary built.")
            logger.debug(f"Returning {model_type} rvh dictionary...")
            return rvh
        elif model_type == "MOHYSE":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvh dictionary...")
            rvh = {
                "Subbasins":
                    [
                        ":SubBasins",
                        "  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,       GAUGED",
                        "  :Units     ,          none,          none,   none,          km,         none",
                        "            1,        Broye_Payerne,            -1,   NONE,       _AUTO,     1",
                        ":EndSubBasins"
                    ],
                "HRUs":
                    [
                        ":HRUs",
                        "  :Attributes,  AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID,LAND_USE_CLASS, VEG_CLASS, SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT",
                        "  :Units     ,   km2,         m,      deg,       deg,     none,          none,      none,         none,            none,          none,   deg,    deg",
                        f"            1, {csv_file.loc['area_ch1903plus']['values']},     {csv_file.loc['a0401_eu_dem_v11_e40n20crp_chv1_0']['values']},   {csv_file.loc['lab_y']['values']},     {csv_file.loc['lab_x']['values']},        1,        LU_ALL,   VEG_ALL,    DEFAULT_P,          [NONE],        [NONE],   {csv_file.loc['a0404_eu_dem_v11_e40n20_slp8v1_0']['values']},  {csv_file.loc['a0407_eu_dem_v11_asp8sm_maskv1_0']['values']}",
                        ":EndHRUs"
                    ],
                "Subbasin Properties":
                    [
                        ":SubBasinProperties",
                        "#          1.0 / MOHYSE_PARA_10,   MOHYSE_PARA_9",
                        "   :Parameters,     GAMMA_SCALE,     GAMMA_SHAPE,",
                        "   :Units,                  1/d,               -",
                        f"              1,          {params[param_or_name][model_type]['MOHYSE_Param_10']},                  {params[param_or_name][model_type]['MOHYSE_Param_09']},",
                        ":EndSubBasinProperties"
                    ]
            }
            logger.debug(f"{model_type} rvh dictionary built.")
            logger.debug(f"Returning {model_type} rvh dictionary...")
            return rvh
        else:
            logger.debug(f"None of if/elif statements evaluated to True, since model type is: {model_type}")
            print(f"Model type {model_type} currently not supported.")

    if generate_rvi is True:
        logger.debug("generate_rvi evaluated to True...")
        logger.debug("Starting if-tree for model_type...")
        if model_type == "GR4J":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvi dictionary...")
            rvi = {
                "Model Organisation":
                    [
                        ":StartDate             1981-01-01 00:00:00",
                        ":EndDate               2019-09-30 00:00:00",
                        ":TimeStep              1.0",
                        ":Method                ORDERED_SERIES",
                        ":RunName               Broye_GR4J",
                        ":EvaluationPeriod CALIBRATION 1981-01-01 2000-01-01",
                        ":EvaluationPeriod VALIDATION 2000-01-02 2019-09-30"
                    ],
                "Model Options":
                    [
                        ":SoilModel             SOIL_MULTILAYER  4",
                        ":Routing               ROUTE_NONE",
                        ":CatchmentRoute        ROUTE_DUMP",
                        ":Evaporation           PET_HAMON",
                        ":RainSnowFraction      RAINSNOW_DINGMAN",
                        ":PotentialMeltMethod   POTMELT_DEGREE_DAY",
                        ":OroTempCorrect        OROCORR_SIMPLELAPSE",
                        ":OroPrecipCorrect      OROCORR_SIMPLELAPSE"
                    ],
                "Soil Layer Alias Definitions":
                    [
                        ":Alias PRODUCT_STORE      SOIL[0]",
                        ":Alias ROUTING_STORE      SOIL[1]",
                        ":Alias TEMP_STORE         SOIL[2]",
                        ":Alias GW_STORE           SOIL[3]"
                    ],
                "Hydrologic Process Order":
                    [
                        ":HydrologicProcesses",
                        "   :Precipitation            PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE",
                        "   :SnowTempEvolve           SNOTEMP_NEWTONS    SNOW_TEMP",
                        "   :SnowBalance              SNOBAL_CEMA_NIEGE  SNOW            PONDED_WATER",
                        "   :OpenWaterEvaporation     OPEN_WATER_EVAP    PONDED_WATER    ATMOSPHERE     	 # Pn",
                        "   :Infiltration             INF_GR4J           PONDED_WATER    MULTIPLE       	 # Ps-",
                        "   :SoilEvaporation          SOILEVAP_GR4J      PRODUCT_STORE   ATMOSPHERE     	 # Es",
                        "   :Percolation              PERC_GR4J          PRODUCT_STORE   TEMP_STORE     	 # Perc",
                        "       :Flush                    RAVEN_DEFAULT      SURFACE_WATER   TEMP_STORE     	 # Pn-Ps",
                        "   :Split                    RAVEN_DEFAULT      TEMP_STORE      CONVOLUTION[0] CONVOLUTION[1] 0.9  # Split Pr",
                        "   :Convolve                 CONVOL_GR4J_1      CONVOLUTION[0]  ROUTING_STORE  	 # Q9",
                        "   :Convolve                 CONVOL_GR4J_2      CONVOLUTION[1]  TEMP_STORE     	 # Q1",
                        "   :Percolation              PERC_GR4JEXCH      ROUTING_STORE   GW_STORE       	 # F(x1)",
                        "   :Percolation              PERC_GR4JEXCH2     TEMP_STORE      GW_STORE       	 # F(x1)",
                        "       :Flush                    RAVEN_DEFAULT      TEMP_STORE      SURFACE_WATER  	 # Qd",
                        "   :Baseflow                 BASE_GR4J          ROUTING_STORE   SURFACE_WATER  	 # Qr",
                        ":EndHydrologicProcesses"
                    ],
                "Output Options":
                    [
                        "#  :EvaluationMetrics NASH_SUTCLIFFE KLING_GUPTA"
                    ]
            }
            logger.debug(f"{model_type} rvi dictionary built.")
            logger.debug(f"Returning {model_type} rvi dictionary...")
            return rvi
        elif model_type == "HYMOD":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvi dictionary...")
            rvi = {
                "Model Organisation":
                    [
                        ":StartDate          1981-01-01 00:00:00",
                        ":EndDate            2019-09-30 00:00:00",
                        ":TimeStep           1.0",
                        ":Method             ORDERED_SERIES",
                        ":RunName            Broye_HYMOD"
                    ],
                "Model Options":
                    [
                        ":Routing             ROUTE_NONE",
                        ":CatchmentRoute      ROUTE_RESERVOIR_SERIES",
                        ":Evaporation         PET_HAMON",
                        ":OW_Evaporation      PET_HAMON",
                        ":SWRadiationMethod   SW_RAD_NONE",
                        ":LWRadiationMethod   LW_RAD_NONE",
                        ":CloudCoverMethod    CLOUDCOV_NONE",
                        ":RainSnowFraction    RAINSNOW_THRESHOLD",
                        ":PotentialMeltMethod POTMELT_DEGREE_DAY",
                        ":PrecipIceptFract    PRECIP_ICEPT_NONE",
                        ":SoilModel           SOIL_MULTILAYER 2",
                        ":EvaluationPeriod   CALIBRATION   1981-01-01   2000-01-01",
                        ":EvaluationPeriod   VALIDATION    2000-01-02   2019-09-30"
                    ],
                "Hydrologic Process Order":
                    [
                        ":HydrologicProcesses",
                        "   :Precipitation     PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE",
                        "   :SnowBalance       SNOBAL_SIMPLE_MELT SNOW            PONDED_WATER",
                        "   :Infiltration      INF_PDM            PONDED_WATER    MULTIPLE",
                        "#  :Flush            RAVEN_DEFAULT      SURFACE_WATER   SOIL[1]   HYMOD_PARAM_9=ALPHA",
                        f"  :Flush             RAVEN_DEFAULT      SURFACE_WATER   SOIL[1]          {params[param_or_name][model_type]['HYMOD_Param_09']}",
                        "   :SoilEvaporation   SOILEVAP_PDM       SOIL[0]         ATMOSPHERE",
                        "   :Baseflow          BASE_LINEAR        SOIL[1]         SURFACE_WATER",
                        ":EndHydrologicProcesses"
                    ],
                "Output Options":
                    [
                        "#  :EvaluationMetrics NASH_SUTCLIFFE KLING_GUPTA"
                    ]
            }
            logger.debug(f"{model_type} rvi dictionary built.")
            logger.debug(f"Returning {model_type} rvi dictionary...")
            return rvi
        elif model_type == "HBV":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvi dictionary...")
            rvi = {
                "Model Organisation":
                    [
                        ":StartDate             1981-01-01 00:00:00",
                        ":EndDate               2019-09-30 00:00:00",
                        ":TimeStep              1.0",
                        ":RunName               Broye_HBV",
                        ":EvaluationPeriod CALIBRATION 1981-01-01 2000-01-01",
                        ":EvaluationPeriod VALIDATION 2000-01-02 2019-09-30",
                    ],
                "Model Options":
                    [
                        ":Routing             	    ROUTE_NONE",
                        ":CatchmentRoute      	    TRIANGULAR_UH",
                        ":Evaporation         	    PET_FROMMONTHLY",
                        ":OW_Evaporation      	    PET_FROMMONTHLY",
                        ":SWRadiationMethod   	    SW_RAD_DEFAULT",
                        ":SWCloudCorrect      	    SW_CLOUD_CORR_NONE",
                        ":SWCanopyCorrect     	    SW_CANOPY_CORR_NONE",
                        ":LWRadiationMethod   	    LW_RAD_DEFAULT",
                        ":RainSnowFraction    	    RAINSNOW_HBV",
                        ":PotentialMeltMethod 	    POTMELT_HBV",
                        ":OroTempCorrect      	    OROCORR_HBV",
                        ":OroPrecipCorrect    	    OROCORR_HBV",
                        ":OroPETCorrect       	    OROCORR_HBV",
                        ":CloudCoverMethod    	    CLOUDCOV_NONE",
                        ":PrecipIceptFract    	    PRECIP_ICEPT_USER",
                        ":MonthlyInterpolationMethod MONTHINT_LINEAR_21",
                        ":SoilModel                  SOIL_MULTILAYER 3"
                    ],
                "Soil Alias Layer Definitions":
                    [
                        ":Alias       FAST_RESERVOIR SOIL[1]",
                        ":Alias       SLOW_RESERVOIR SOIL[2]",
                        ":LakeStorage SLOW_RESERVOIR"
                    ],
                "Hydrologic Process Order":
                    [
                        ":HydrologicProcesses",
                        "   :SnowRefreeze      FREEZE_DEGREE_DAY  SNOW_LIQ        SNOW",
                        "   :Precipitation     PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE",
                        "   :CanopyEvaporation CANEVP_ALL         CANOPY          ATMOSPHERE",
                        ":  CanopySnowEvap    CANEVP_ALL         CANOPY_SNOW     ATMOSPHERE",
                        "   :SnowBalance       SNOBAL_SIMPLE_MELT SNOW            SNOW_LIQ",
                        "       :-->Overflow     RAVEN_DEFAULT      SNOW_LIQ        PONDED_WATER",
                        "   :Flush             RAVEN_DEFAULT      PONDED_WATER    GLACIER",
                        "       :-->Conditional HRU_TYPE IS GLACIER",
                        "   :GlacierMelt       GMELT_HBV          GLACIER_ICE     GLACIER",
                        "   :GlacierRelease    GRELEASE_HBV_EC    GLACIER         SURFACE_WATER",
                        "   :Infiltration      INF_HBV            PONDED_WATER    MULTIPLE",
                        "   :Flush             RAVEN_DEFAULT      SURFACE_WATER   FAST_RESERVOIR",
                        "       :-->Conditional HRU_TYPE IS_NOT GLACIER",
                        "   :SoilEvaporation   SOILEVAP_HBV       SOIL[0]         ATMOSPHERE",
                        "   :CapillaryRise     RISE_HBV           FAST_RESERVOIR 	SOIL[0]",
                        "   :LakeEvaporation   LAKE_EVAP_BASIC    SLOW_RESERVOIR  ATMOSPHERE",
                        "   :Percolation       PERC_CONSTANT      FAST_RESERVOIR 	SLOW_RESERVOIR",
                        "   :Baseflow          BASE_POWER_LAW     FAST_RESERVOIR  SURFACE_WATER",
                        "   :Baseflow          BASE_LINEAR        SLOW_RESERVOIR  SURFACE_WATER",
                        ":EndHydrologicProcesses"
                    ],
                "Output Options":
                    [
                        "#:EvaluationMetrics NASH_SUTCLIFFE"
                    ]
            }
            logger.debug(f"{model_type} rvi dictionary built.")
            logger.debug(f"Returning {model_type} rvi dictionary...")
            return rvi
        elif model_type == "HMETS":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvi dictionary...")
            rvi = {
                "Model Organisation":
                    [
                        ":StartDate               1954-01-01 00:00:00",
                        ":Duration                2010-12-31 00:00:00",
                        ":TimeStep                1.0",
                        ":Method                  ORDERED_SERIES",
                        ":RunName                 raven_broye_hmets",
                        ":EvaluationPeriod CALIBRATION 1981-01-01 2000-01-01",
                        ":EvaluationPeriod VALIDATION 2000-01-02 2019-09-30"
                    ],
                "Model Options":
                    [
                        ":PotentialMeltMethod     POTMELT_HMETS",
                        ":RainSnowFraction        RAINSNOW_DATA",
                        ":Evaporation             PET_DATA",
                        "#:Evaporation            PET_OUDIN",
                        ":CatchmentRoute          ROUTE_DUMP",
                        ":Routing                 ROUTE_NONE",
                        ":SoilModel               SOIL_TWO_LAYER"
                    ],
                "Alias Definitions":
                    [
                        ":Alias DELAYED_RUNOFF CONVOLUTION[1]"
                    ],
                "Hydrologic Process Order":
                    [
                        ":HydrologicProcesses",
                        "   :SnowBalance     SNOBAL_HMETS    MULTIPLE     MULTIPLE",
                        "   :Precipitation   RAVEN_DEFAULT   ATMOS_PRECIP MULTIPLE",
                        "   :Infiltration    INF_HMETS       PONDED_WATER MULTIPLE",
                        "   :Overflow      OVERFLOW_RAVEN  SOIL[0]      DELAYED_RUNOFF",
                        "   :Baseflow        BASE_LINEAR     SOIL[0]      SURFACE_WATER   # interflow, really",
                        "   :Percolation     PERC_LINEAR     SOIL[0]      SOIL[1]         # recharge",
                        "   :Overflow      OVERFLOW_RAVEN  SOIL[1]      DELAYED_RUNOFF",
                        "   :SoilEvaporation SOILEVAP_ALL    SOIL[0]      ATMOSPHERE      # AET",
                        "   :Convolve        CONVOL_GAMMA    CONVOLUTION[0] SURFACE_WATER #'surface runoff'",
                        "   :Convolve        CONVOL_GAMMA_2  DELAYED_RUNOFF SURFACE_WATER #'delayed runoff'",
                        "   :Baseflow        BASE_LINEAR     SOIL[1]      SURFACE_WATER",
                        ":EndHydrologicProcesses"
                    ],
                "Output Options":
                    [
                        "#:EvaluationMetrics NASH_SUTCLIFFE"
                    ]
            }
            logger.debug(f"{model_type} rvp dictionary built.")
            logger.debug(f"Returning {model_type} rvi dictionary...")
            return rvi
        elif model_type == "MOHYSE":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvi dictionary...")
            rvi = {
                "Model Organisation":
                    [
                        ":StartDate               1954-01-01 00:00:00",
                        ":Duration                2010-12-31 00:00:00",
                        ":TimeStep                1.0",
                        ":Method                  ORDERED_SERIES",
                        ":RunName                 Broye_MOHYSE",
                        ":EvaluationPeriod CALIBRATION 1981-01-01 2000-01-01",
                        ":EvaluationPeriod VALIDATION 2000-01-02 2019-09-30"
                    ],
                "Model Options":
                    [
                        ":SoilModel             SOIL_TWO_LAYER",
                        ":PotentialMeltMethod   POTMELT_DEGREE_DAY",
                        ":Routing               ROUTE_NONE",
                        ":CatchmentRoute        ROUTE_GAMMA_CONVOLUTION",
                        ":Evaporation           PET_MOHYSE",
                        ":DirectEvaporation",
                        ":RainSnowFraction      RAINSNOW_DATA"
                    ],
                "Alias Definitions":
                    [
                        "# :Alias MOHYSE_PARA_1      1.5589    # :GlobalParameter         MOHYSE_PET_COEFF",
                        "# :Alias MOHYSE_PARA_2	    0.9991    # LandUseParameterList --> AET_COEFF",
                        "# :Alias MOHYSE_PARA_3	    2.1511    # LandUseParameterList --> MELT_FACTOR",
                        "# :Alias MOHYSE_PARA_4	   -1.6101    # LandUseParameterList --> DD_MELT_TEMP",
                        "# :Alias MOHYSE_PARA_5	    0.5000    # SoilProfiles         --> thickness of TOPSOIL (in mm????? must be m!!!)",
                        "# :Alias MOHYSE_PARA_6	    0.1050    # SoilParameterList    --> PERC_COEFF (TOPSOIL)",
                        "# :Alias MOHYSE_PARA_7	    0.0533    # SoilParameterList    --> BASEFLOW_COEFF (TOPSOIL)",
                        "# :Alias MOHYSE_PARA_8	    0.0132    # SoilParameterList    --> BASEFLOW_COEFF (GWSOIL)",
                        "# :Alias MOHYSE_PARA_9	    1.0474    # :SubBasinProperties  --> GAMMA_SHAPE",
                        "# :Alias MOHYSE_PARA_10	    7.9628    # :SubBasinProperties  --> TIME_CONC = MOHYSE_PARA_10 / 0.3 = 26.542666666"
                    ],
                "Hydrologic Process Order":
                    [
                        ":HydrologicProcesses",
                        "   :SoilEvaporation  SOILEVAP_LINEAR    SOIL[0]            ATMOSPHERE",
                        "   :SnowBalance      SNOBAL_SIMPLE_MELT SNOW PONDED_WATER",
                        "   :Precipitation    RAVEN_DEFAULT      ATMOS_PRECIP       MULTIPLE",
                        "   :Infiltration     INF_HBV            PONDED_WATER       SOIL[0]",
                        "   :Baseflow         BASE_LINEAR        SOIL[0]            SURFACE_WATER",
                        "   :Percolation      PERC_LINEAR        SOIL[0]            SOIL[1]",
                        "   :Baseflow         BASE_LINEAR        SOIL[1]            SURFACE_WATER",
                        ":EndHydrologicProcesses"
                    ],
                "Output Options":
                    [
                        "#:EvaluationMetrics NASH_SUTCLIFFE"
                    ]
            }
            logger.debug(f"{model_type} rvi dictionary built.")
            logger.debug(f"Returning {model_type} rvi dictionary...")
            return rvi
        else:
            logger.debug(f"None of if/elif statements evaluated to True, since model type is: {model_type}")
            print(f"Model type {model_type} currently not supported.")

    if generate_rvc is True:
        logger.debug("generate_rvc evaluated to True...")
        logger.debug("Starting if-tree for model_type...")
        if model_type == "GR4J":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvc dictionary...")
            rvc = {
                "Soil Profiles":
                    [
                        f"# SOIL[0] = {params[param_or_name][model_type]['GR4J_X1']} * 1000. / 2.0 (initialize to 1/2 full)",
                        "# SOIL[1] = 0.3m * 1000. / 2.0   (initialize to 1/2 full)"
                    ],
                "HRU States":
                    [
                        ":HRUStateVariableTable",
                        "   :Attributes SOIL[0] SOIL[1]",
                        "   :Units      mm      mm",
                        "   1           264.5   15.0",
                        ":EndHRUStateVariableTable"
                    ]
            }
            logger.debug(f"{model_type} rvc dictionary built.")
            logger.debug(f"Returning {model_type} rvc dictionary...")
            return rvc
        elif model_type == "HYMOD":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvc dictionary...")
            rvc = {
                "Empty":
                    [
                        "# Nothing to set here."
                    ]
            }
            logger.debug(f"{model_type} rvc dictionary built.")
            logger.debug(f"Returning {model_type} rvc dictionary...")
            return rvc
        elif model_type == "HBV":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvc dictionary...")
            # TODO: Which param is HBV_PARA_???
            rvc = {
                "Basin":
                    [
                        ":BasinInitialConditions",
                        ":Attributes, ID,              Q",
                        ":Units,      none,         m3/s",
                        "#                  HBV_PARA_???",
                        "1,             1.0",
                        ":EndBasinInitialConditions"
                    ],
                "Lower Groundwater Storage":
                    [
                        "# Initial Lower groundwater storage - for each HRU",
                        "",
                        ":InitialConditions SOIL[2]",
                        "# derived from thickness: HBV_PARA_17 [m] * 1000.0 / 2.0",
                        f"{params[param_or_name][model_type]['HBV_Param_17']}",
                        ":EndInitialConditions"
                    ]
            }
            logger.debug(f"{model_type} rvc dictionary built.")
            logger.debug(f"Returning {model_type} rvc dictionary...")
            return rvc
        elif model_type == "HMETS":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvc dictionary...")
            rvc = {
                "Initial Storage":
                    [
                        "# initialize to 1/2 full",
                        "# x(20)/2",
                        ":UniformInitialConditions SOIL[0] 155.36055",
                        "# x(21)/2",
                        ":UniformInitialConditions SOIL[1] 458.09735"
                    ],
                "HRUs":
                    [
                        ":HRUStateVariableTable",
                        ":Attributes SOIL[0] SOIL[1]",
                        ":Units mm mm",
                        "1 155.36055 458.09735",
                        ":EndHRUStateVariableTable",
                    ]
            }
            logger.debug(f"{model_type} rvc dictionary built.")
            logger.debug(f"Returning {model_type} rvc dictionary...")
            return rvc
        elif model_type == "MOHYSE":
            logger.debug(f"Model type is {model_type}.")
            logger.debug(f"Building {model_type} rvc dictionary...")
            rvc = {
                "Empty":
                    [
                        "# Nothing to set here."
                    ]
            }
            logger.debug(f"{model_type} rvc dictionary built.")
            logger.debug(f"Returning {model_type} rvc dictionary...")
            return rvc
        else:
            logger.debug(f"None of if/elif statements evaluated to True, since model type is: {model_type}")
            print(f"Model type {model_type} currently not supported.")

    if generate_rvt is True:
        pass



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
    logger.debug(
        f"Trying to read catchment attribute CSV file {Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name)}...")
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
                                              param_or_name="params", generate_rvp=True)
        logger.debug("Wrote .rvp template sections generated by generate_template() to dict template_sections")
    if template_type == "Ostrich":
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Adding .tpl suffix to file_path for Ostrich template file...")
        file_path = Path((str(file_path) + ".tpl"))
        logger.debug(f"New file_path: {file_path}")
        logger.debug(f"Trying to generate .rvp.tpl template sections with function generate_template()...")
        template_sections = generate_template(model_type=model_type, csv_file=csv_file, params=params,
                                              param_or_name="names", generate_rvp=True)
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


# TODO: is 'melt_factor' actually a parameter of GR4J? In the tutorial file, 'annual snow' is called 'CEMANEIGE_X1'...


def write_rvh(model_dir=model_dir, model_type=model_type, data_dir=data_dir,
              project_dir=project_dir, catchment=catchment, model_sub_dir=model_sub_dir, params=default_params,
              template_type="Raven",
              attribute_csv_name="CH-0057_attributes.csv",
              attribute_csv_dir="Hydromap Attributes"):
    logger.debug("Arrived in function write_rvh().")
    assert model_type in supported_models, f"Got model type: {model_type}, which is not supported, check variable \"" \
                                           f"supported_models."
    logger.debug(
        f"Trying to read catchment attribute CSV file {Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name)}...")
    csv_file = pandas.read_csv(Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name), sep=",",
                               skiprows=[8],
                               index_col='attribute_names')
    logger.debug("Attribute catchment attribute CSV file read.")
    file_name: str = f"{catchment}_{model_type}.rvh"
    logger.debug(f".rvh filename set to {file_name}.")
    file_path: Path = Path(project_dir, model_dir, catchment, model_type, model_sub_dir, file_name)
    logger.debug(f".rvh file path set to {file_path}.")
    template_sections = {}
    logger.debug("Empty dict template_sections created.")
    logger.debug("Entering if-tree for template type...")
    if template_type == "Raven":
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Trying to generate .rvh template sections with function generate_template()...")
        template_sections = generate_template(model_type=model_type, csv_file=csv_file, params=params,
                                              param_or_name="params", generate_rvh=True)
        logger.debug("Wrote .rvh template sections generated by generate_template() to dict template_sections")
    if template_type == "Ostrich":
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Adding .tpl suffix to file_path for Ostrich template file...")
        file_path = Path((str(file_path) + ".tpl"))
        logger.debug(f"New file_path: {file_path}")
        logger.debug(f"Trying to generate .rvh.tpl template sections with function generate_template()...")
        template_sections = generate_template(model_type=model_type, csv_file=csv_file, params=params,
                                              param_or_name="names", generate_rvh=True)
        logger.debug("Wrote .rvh.tpl template sections generated by generate_template() to dict template_sections")

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
        logger.debug("Function write_rvh() finished.")


def write_rvi(model_dir=model_dir, model_type=model_type, data_dir=data_dir,
              project_dir=project_dir, catchment=catchment, model_sub_dir=model_sub_dir, params=default_params,
              template_type="Raven",
              attribute_csv_name="CH-0057_attributes.csv",
              attribute_csv_dir="Hydromap Attributes"):
    logger.debug("Arrived in function write_rvi().")
    assert model_type in supported_models, f"Got model type: {model_type}, which is not supported, check variable \"" \
                                           f"supported_models."
    logger.debug(
        f"Trying to read catchment attribute CSV file {Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name)}...")
    csv_file = pandas.read_csv(Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name), sep=",",
                               skiprows=[8],
                               index_col='attribute_names')
    logger.debug("Attribute catchment attribute CSV file read.")
    file_name: str = f"{catchment}_{model_type}.rvi"
    logger.debug(f".rvi filename set to {file_name}.")
    file_path: Path = Path(project_dir, model_dir, catchment, model_type, model_sub_dir, file_name)
    logger.debug(f".rvi file path set to {file_path}.")
    template_sections = {}
    logger.debug("Empty dict template_sections created.")
    logger.debug("Entering if-tree for template type...")
    if template_type == "Raven":
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Trying to generate .rvi template sections with function generate_template()...")
        template_sections = generate_template(model_type=model_type, csv_file=csv_file, params=params,
                                              param_or_name="params", generate_rvi=True)
        logger.debug("Wrote .rvi template sections generated by generate_template() to dict template_sections")
    if template_type == "Ostrich":
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Adding .tpl suffix to file_path for Ostrich template file...")
        file_path = Path((str(file_path) + ".tpl"))
        logger.debug(f"New file_path: {file_path}")
        logger.debug(f"Trying to generate .rvi.tpl template sections with function generate_template()...")
        template_sections = generate_template(model_type=model_type, csv_file=csv_file, params=params,
                                              param_or_name="names", generate_rvi=True)
        logger.debug("Wrote .rvi.tpl template sections generated by generate_template() to dict template_sections")

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
        logger.debug("Function write_rvi() finished.")


def write_rvc(model_dir: str=model_dir, model_type=model_type, data_dir=data_dir,
              project_dir=project_dir, catchment=catchment, model_sub_dir=model_sub_dir, params=default_params,
              template_type="Raven",
              attribute_csv_name="CH-0057_attributes.csv",
              attribute_csv_dir="Hydromap Attributes"):
    """

    :param str model_dir: Directory name of current model
    :param model_type:
    :param data_dir:
    :param project_dir:
    :param catchment:
    :param model_sub_dir:
    :param params:
    :param template_type:
    :param attribute_csv_name:
    :param attribute_csv_dir:
    """
    logger.debug("Arrived in function write_rvc().")
    assert model_type in supported_models, f"Got model type: {model_type}, which is not supported, check variable \"" \
                                           f"supported_models."
    logger.debug(
        f"Trying to read catchment attribute CSV file {Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name)}...")
    csv_file = pandas.read_csv(Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name), sep=",",
                               skiprows=[8],
                               index_col='attribute_names')
    logger.debug("Attribute catchment attribute CSV file read.")
    file_name: str = f"{catchment}_{model_type}.rvc"
    logger.debug(f".rvc filename set to {file_name}.")
    file_path: Path = Path(project_dir, model_dir, catchment, model_type, model_sub_dir, file_name)
    logger.debug(f".rvc file path set to {file_path}.")
    template_sections = {}
    logger.debug("Empty dict template_sections created.")
    logger.debug("Entering if-tree for template type...")
    if template_type == "Raven":
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Trying to generate .rvc template sections with function generate_template()...")
        template_sections = generate_template(model_type=model_type, csv_file=csv_file, params=params,
                                              param_or_name="params", generate_rvc=True)
        logger.debug("Wrote .rvc template sections generated by generate_template() to dict template_sections")
    if template_type == "Ostrich":
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Adding .tpl suffix to file_path for Ostrich template file...")
        file_path = Path((str(file_path) + ".tpl"))
        logger.debug(f"New file_path: {file_path}")
        logger.debug(f"Trying to generate .rvc.tpl template sections with function generate_template()...")
        template_sections = generate_template(model_type=model_type, csv_file=csv_file, params=params,
                                              param_or_name="names", generate_rvc=True)
        logger.debug("Wrote .rvc.tpl template sections generated by generate_template() to dict template_sections")

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
        logger.debug("Function write_rvc() finished.")


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
