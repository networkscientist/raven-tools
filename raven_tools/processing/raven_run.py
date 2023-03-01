"""
Tools to generate Raven .rv* files needed to run Raven models.
"""
import logging
import os
from pathlib import Path

import pandas

from raven_tools import config

try:
    logger = logging.getLogger(__name__)
    logger.debug("TEEEST")

    logger.debug(f"CWD: {os.getcwd()}")
    logger.debug('Trying to read project_config.yaml file')
except:
    logger.exception("Error loading logger.")

try:
    conf = config.variables.project_config
    model_dir = conf['ModelDir']
    model_type = conf['ModelName']
    data_dir = conf['DataDir']
    project_dir = conf['ProjectDir']
    catchment = conf['Catchment']
    catchment_id = conf['CatchmentID']
    gauge_lat = conf['GaugeLat']
    gauge_lon = conf['GaugeLon']
    model_sub_dir = conf['ModelSubDir']
    grid_weights_file = conf['GridWeights']
    author = conf['Author']
    generation_date = conf['Date']
    poetry_location = conf['PoetryEnvLocation']
except:
    logger.exception("Error getting project_config from __init__.py!")

try:
    default_params = config.variables.default_params
except:
    logger.exception("Error getting default_params from config.py")

subsection_header_line: str = "#-----------------------------------------------------------------"
newline: str = "\n"


def get_catchment_info(csv_file):
    pandas.read_csv(csv_file)


def create_header(author=conf['Author'], catchment=catchment, model=model_type):
    try:
        header_line = "#########################################################################"
        file_type = ":FileType          rvt ASCII Raven 3.5"
        author_line = f":WrittenBy         {author}"
        creation_date = ":CreationDate      April 2022"
        description = [
            "#",
            f"# Emulation of {model} simulation of {catchment}",
            "#------------------------------------------------------------------------"]
        header = [header_line, file_type, author_line, creation_date, *description]
        return header
    except NameError:
        print("Probably project_config file could not be found...")
        pass


def forcing_block(start_year: int, end_year: int, catchment: str):
    """Create Dictionary of forcing data to write in RVT file.

    This function creates a Dictionary of forcing data to be written into an RVT file. From a start and end year,
    it creates the relevant directory Paths, which are based on Swiss gridded input data. Further parameters,
    as required by RAVEN are added and finally, the forcings data block is returned as a Dictionary.

    Args:
        start_year : int
            Start year of forcings data files
        end_year : int
            End year of forcings data files

    Returns:
        forcing_data: dict[str, list[str]]
            The forcing data block

    """

    forcing_data: dict[str, list[str]] = {
        'Rainfall': [
            ":GriddedForcing           Rainfall",
            "    :ForcingType          RAINFALL",
            f"    :FileNameNC           data_obs/RhiresD_v2.0_swiss.lv95/merged/RhiresD_ch01h.swiss.lv95_{start_year}01010000_{end_year}12310000_{catchment}_clipped.nc",
            "    :VarNameNC            RhiresD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            f"    :RedirectToFile       data_obs/{grid_weights_file}",
            ":EndGriddedForcing"],
        'Average Temperature': [
            ":GriddedForcing           Average Temperature",
            "    :ForcingType          TEMP_AVE",
            f"    :FileNameNC           data_obs/TabsD_v2.0_swiss.lv95/merged/TabsD_ch01r.swiss.lv95_{start_year}01010000_{end_year}12310000_{catchment}_clipped.nc",
            "    :VarNameNC            TabsD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            f"    :RedirectToFile       data_obs/{grid_weights_file}",
            ":EndGriddedForcing"],
        'Maximum Temperature': [
            ":GriddedForcing           Maximum Temperature",
            "    :ForcingType          TEMP_MAX",
            f"    :FileNameNC           data_obs/TmaxD_v2.0_swiss.lv95/merged/TmaxD_ch01r.swiss.lv95_{start_year}01010000_{end_year}12310000_{catchment}_clipped.nc",
            "    :VarNameNC            TmaxD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            f"    :RedirectToFile       data_obs/{grid_weights_file}",
            ":EndGriddedForcing"],
        'Minimum Temperature': [
            ":GriddedForcing           Minimum Temperature",
            "    :ForcingType          TEMP_MIN",
            f"    :FileNameNC           data_obs/TminD_v2.0_swiss.lv95/merged/TminD_ch01r.swiss.lv95_{start_year}01010000_{end_year}12310000_{catchment}_clipped.nc",
            "    :VarNameNC            TminD",
            "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
            f"    :RedirectToFile       data_obs/{grid_weights_file}",
            ":EndGriddedForcing"
        ]
    }

    # forcing_header = f"# Years: {start_year} to {end_year}\n"

    return forcing_data


def write_rvt(start_year: int,
              end_year: int,
              model_dir=model_dir,
              model_type=model_type,
              project_dir=project_dir,
              catchment=catchment,
              catchment_id=catchment_id,
              gauge_lat=gauge_lat,
              gauge_lon=gauge_lon,
              model_sub_dir=model_sub_dir,
              author=conf['Author'],
              gauge_short_code="DefAult",
              station_elevation="100"):
    """Write to Raven .rvt file.

    Args:
        model_sub_dir:
        catchment:
        project_dir:
        start_year : int
            Start year of forcings data files
        end_year : int
            End year of forcings data files
        model_dir: Root directory of .rvX files
        model_type: Name of the .rvt file to be written

    """
    logger.debug("Entered write_rvt function.")
    file_name: str = f"{catchment}_{model_type}.rvt"
    file_path: Path = Path(model_dir, model_sub_dir, file_name)
    logger.debug(f"file_name = {file_name}")
    logger.debug(f"file_path = {file_path}")

    gauge = [
        f":Gauge {gauge_short_code}\n",
        f"  :Latitude    {gauge_lat}\n",
        f"  :Longitude {gauge_lon}\n",
        f"  :Elevation  {station_elevation}\n",
        ":EndGauge\n"
    ]

    flow_observation = [
        "# observed streamflow\n",
        ":RedirectToFile data_obs/BroPay_Q_2034_daily.rvt"
    ]
    with open(file_path, 'w') as ff:
        ff.writelines(f"{line}{newline}" for line in create_header(author, catchment, model_type))
        ff.write(f"# meteorological forcings\n")
        for f in forcing_block(start_year, end_year, catchment).values():
            for t in f:
                ff.write(f"{t}\n")

        ff.writelines(gauge)
        ff.writelines(flow_observation)


def generate_template_rvx(csv_file=None, model_type=model_type, params=default_params, param_or_name="names") -> dict:
    """Generates template text which can be written to .rvX file.

        Args:
            model_type : str
                Name of the model type
            csv_file : str
                File path of the csv file
            params : dict
                Parameters for the model
            param_or_name : str
                Name of the parameter or parameter value

        Returns:
            rvx_params : dict
                Dictionary containing the parameters for the .rvX file
    """

    logger.debug("Entered generate_template() function.")
    assert model_type in config.variables.supported_models, f"model_type expected GR4J, HYMOD, HMETS, HBV or MOHYSE, got {model_type} instead "
    logger.debug("model_type is in the list of supported models.")

    logger.debug("Trying to create rvx_params dictionary...")
    rvx_params = \
        {
            "GR4J":
                {
                    "rvp":
                        {
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
                                    f"   DEFAULT_P, 4, SOIL_PROD, {params[param_or_name]['GR4J']['GR4J_X1']}, SOIL_ROUT, 0.300, SOIL_TEMP, 1.000, SOIL_GW, 1.000,",
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
                                    f":GlobalParameter AIRSNOW_COEFF     {params[param_or_name]['GR4J']['Airsnow_Coeff']} # [1/d] = 1.0 - CEMANEIGE_X2 = 1.0 - x6",
                                    ":GlobalParameter AVG_ANNUAL_SNOW    16.9 # [mm]  =       CEMANEIGE_X1 =       x5",
                                    "#:GlobalParameter PRECIP_LAPSE     0.0004 I assume not necessary for gridded data",
                                    "#:GlobalParameter ADIABATIC_LAPSE  0.0065 not necessary for gridded data"
                                ],
                            "Soil Parameters":
                                [
                                    ":SoilParameterList",
                                    "   :Parameters, POROSITY, GR4J_X3, GR4J_X2",
                                    "   :Units, none, mm, mm / d",
                                    f"   [DEFAULT], 1.0, {params[param_or_name]['GR4J']['GR4J_X3']}, {params[param_or_name]['GR4J']['GR4J_X2']}",
                                    ":EndSoilParameterList"
                                ],
                            "Land Use Parameters":
                                [
                                    ":LandUseParameterList",
                                    "   :Parameters, GR4J_X4, MELT_FACTOR",
                                    "   :Units, d, mm / d / C",
                                    f"   [DEFAULT], {params[param_or_name]['GR4J']['GR4J_X4']}, {params[param_or_name]['GR4J']['Melt_Factor']}",
                                    ":EndLandUseParameterList"
                                ]
                        },
                    "rvh":
                        {"Subbasins":
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
                                ]},
                    "rvi":
                        {"Model Organisation":
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
                        },
                    "rvc":
                        {
                            "Soil Profiles":
                                [
                                    f"# SOIL[0] = {params[param_or_name]['GR4J']['GR4J_X1']} * 1000. / 2.0 (initialize to 1/2 full)",
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
                        },
                    "rvt":
                        {}
                },
            "HYMOD":
                {
                    "rvp":
                        {
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
                                    f"   DEFAULT_P, 2, TOPSOIL, {params[param_or_name]['HYMOD']['HYMOD_Param_02']}, GWSOIL, 10.0",
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
                                    f":GlobalParameter RAINSNOW_TEMP {params[param_or_name]['HYMOD']['HYMOD_Param_03']}",
                                    "   #:GlobalParameter      RAINSNOW_TEMP    HYMOD_PARA_3"
                                ],
                            "Soil Parameters":
                                [
                                    ":SoilParameterList",
                                    "   :Parameters, POROSITY, PET_CORRECTION, BASEFLOW_COEFF,",
                                    "   :Units, -, -, 1 / d,",
                                    "       # TOPSOIL,            1.0 ,    HYMOD_PARA_8,               0.0,",
                                    "       #  GWSOIL,            1.0 ,             1.0,   HYMOD_PARA_4=Ks,",
                                    f"       TOPSOIL, 1.0, {params[param_or_name]['HYMOD']['HYMOD_Param_08']}, 0.0,",
                                    f"       GWSOIL, 1.0, 1.0, {params[param_or_name]['HYMOD']['HYMOD_Param_04']},",
                                    ":EndSoilParameterList"
                                ],
                            "Land Use Parameters":
                                [
                                    ":LandUseParameterList",
                                    "   :Parameters, MELT_FACTOR, DD_MELT_TEMP, PDM_B,",
                                    "   :Units, mm / d / K, degC, -,",
                                    "       # [DEFAULT],    HYMOD_PARA_5,    HYMOD_PARA_6,  HYMOD_PARA_7=Bexp,",
                                    f"       [DEFAULT], {params[param_or_name]['HYMOD']['HYMOD_Param_05']}, {params[param_or_name]['HYMOD']['HYMOD_Param_06']}, {params[param_or_name]['HYMOD']['HYMOD_Param_07']},",
                                    ":EndLandUseParameterList"
                                ]
                        },
                    "rvh":
                        {"Subbasins":
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
                                    f"              1,          {params[param_or_name]['HYMOD']['HYMOD_Param_01']},                  3,",
                                    ":EndSubBasinProperties"
                                ]
                        },
                    "rvi":
                        {"Model Organisation":
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
                                    "\t:Precipitation     PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE",
                                    "   :SnowBalance       SNOBAL_SIMPLE_MELT SNOW            PONDED_WATER",
                                    "   :Infiltration      INF_PDM            PONDED_WATER    MULTIPLE",
                                    "#  :Flush            RAVEN_DEFAULT      SURFACE_WATER   SOIL[1]   HYMOD_PARAM_9=ALPHA",
                                    f"  :Flush             RAVEN_DEFAULT      SURFACE_WATER   SOIL[1]          {params[param_or_name]['HYMOD']['HYMOD_Param_09']}",
                                    "   :SoilEvaporation   SOILEVAP_PDM       SOIL[0]         ATMOSPHERE",
                                    "   :Baseflow          BASE_LINEAR        SOIL[1]         SURFACE_WATER",
                                    ":EndHydrologicProcesses"
                                ],
                            "Output Options":
                                [
                                    "#  :EvaluationMetrics NASH_SUTCLIFFE KLING_GUPTA"
                                ]
                        },
                    "rvc":
                        {"Empty":
                            [
                                "# Nothing to set here."
                            ]
                        },
                    "rvt":
                        {}

                },
            "HMETS":
                {
                    "rvp":
                        {
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
                                    ":SoilProfiles",
                                    "   LAKE, 0",
                                    "   ROCK, 0",
                                    "   # DEFAULT_P, 2, TOPSOIL,          x(20)/1000, PHREATIC,         x(21)/1000,",
                                    f"  DEFAULT_P, 2, TOPSOIL,     {params[param_or_name]['HMETS']['HMETS_Param_20b']}, PHREATIC,     {params[param_or_name]['HMETS']['HMETS_Param_21b']},",
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
                                    f":GlobalParameter  SNOW_SWI_MIN {params[param_or_name]['HMETS']['HMETS_Param_09']} # x(9)",
                                    f":GlobalParameter  SNOW_SWI_MAX {params[param_or_name]['HMETS']['HMETS_Param_09b']} # x(9)+x(10)",
                                    f":GlobalParameter  SWI_REDUCT_COEFF {params[param_or_name]['HMETS']['HMETS_Param_11']} # x(11)",
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
                                    f"      [DEFAULT],  {params[param_or_name]['HMETS']['HMETS_Param_05']}, {params[param_or_name]['HMETS']['HMETS_Param_05b']},  {params[param_or_name]['HMETS']['HMETS_Param_07']},  {params[param_or_name]['HMETS']['HMETS_Param_08']},  {params[param_or_name]['HMETS']['HMETS_Param_13']},  {params[param_or_name]['HMETS']['HMETS_Param_14']},   {params[param_or_name]['HMETS']['HMETS_Param_12']},     {params[param_or_name]['HMETS']['HMETS_Param_16']},",
                                    "#      x(5),       x(5)+x(6),            x(7),            x(8),           x(13),           x(14),            x(12),              x(16),",
                                    ":EndLandUseParameterList",
                                    "",
                                    ":LandUseParameterList",
                                    "   :Parameters,     GAMMA_SHAPE,     GAMMA_SCALE,    GAMMA_SHAPE2,    GAMMA_SCALE2,",
                                    "   :Units,               -,             1/d,               -,             1/d,",
                                    f"      [DEFAULT],  {params[param_or_name]['HMETS']['HMETS_Param_01']},  {params[param_or_name]['HMETS']['HMETS_Param_02']},  {params[param_or_name]['HMETS']['HMETS_Param_03']},  {params[param_or_name]['HMETS']['HMETS_Param_04']},",
                                    "#      x(1),            x(2),            x(3),            x(4),",
                                    ":EndLandUseParameterList"
                                ],
                            "Soil Parameters":
                                [
                                    ":SoilParameterList",
                                    "   :Parameters,        POROSITY,      PERC_COEFF,  PET_CORRECTION, BASEFLOW_COEFF",
                                    "   :Units,               -,             1/d,               -,            1/d",
                                    f"      TOPSOIL,             1.0,  {params[param_or_name]['HMETS']['HMETS_Param_17']},  {params[param_or_name]['HMETS']['HMETS_Param_15']}, {params[param_or_name]['HMETS']['HMETS_Param_18']}",
                                    f"      PHREATIC,             1.0,             0.0,             0.0, {params[param_or_name]['HMETS']['HMETS_Param_19']}",
                                    "#      TOPSOIL,             1.0,           x(17),           x(15),          x(18)",
                                    "#      PHREATIC,             1.0,             0.0,             0.0,          x(19)",
                                    ":EndSoilParameterList"
                                ]
                        },
                    "rvh":
                        {"Subbasins":
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
                        },
                    "rvi":
                        {"Model Organisation":
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
                        },
                    "rvc":
                        {"Initial Storage":
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
                        },
                    "rvt":
                        {}
                },
            "HBV":
                {
                    "rvp":
                        {
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
                                    ":SoilProfiles",
                                    f"   DEFAULT_P,      3,    TOPSOIL,            {params[param_or_name]['HBV']['HBV_Param_17']},   FAST_RES,    100.0, SLOW_RES,    100.0",
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
                                    f":GlobalParameter RAINSNOW_TEMP       {params[param_or_name]['HBV']['HBV_Param_01']}",
                                    ":GlobalParameter RAINSNOW_DELTA      1.0 #constant",
                                    f"#:GlobalParameter PRECIP_LAPSE     {params[param_or_name]['HBV']['HBV_Param_12']} # I assume not necessary for gridded data, HBV_PARA_12=PCALT",
                                    f"#:GlobalParameter ADIABATIC_LAPSE  {params[param_or_name]['HBV']['HBV_Param_13']} # not necessary for gridded data, HBV_PARA_13=TCALT",
                                    f":GlobalParameter SNOW_SWI  {params[param_or_name]['HBV']['HBV_Param_04']} #HBV_PARA_04"
                                ],
                            "Land Use Parameters":
                                [
                                    ":LandUseParameterList",
                                    "  :Parameters,   MELT_FACTOR, MIN_MELT_FACTOR,   HBV_MELT_FOR_CORR, REFREEZE_FACTOR, HBV_MELT_ASP_CORR",
                                    "  :Units     ,        mm/d/K,          mm/d/K,                none,          mm/d/K,              none",
                                    "  #              HBV_PARA_02,        CONSTANT,         HBV_PARA_18,     HBV_PARA_03,          CONSTANT",
                                    f"    [DEFAULT],  {params[param_or_name]['HBV']['HBV_Param_02']},             2.2,        {params[param_or_name]['HBV']['HBV_Param_18']},    {params[param_or_name]['HBV']['HBV_Param_03']},              0.48",
                                    ":EndLandUseParameterList",
                                    "",
                                    ":LandUseParameterList",
                                    " :Parameters, HBV_MELT_GLACIER_CORR,   HBV_GLACIER_KMIN, GLAC_STORAGE_COEFF, HBV_GLACIER_AG",
                                    " :Units     ,                  none,                1/d,                1/d,           1/mm",
                                    "   #                       CONSTANT,           CONSTANT,        HBV_PARA_19,       CONSTANT,",
                                    f"   [DEFAULT],                  1.64,               0.05,       {params[param_or_name]['HBV']['HBV_Param_19']},           0.05",
                                    ":EndLandUseParameterList"

                                ],
                            "Soil Parameters":
                                [
                                    ":SoilParameterList",
                                    "  :Parameters,                POROSITY,FIELD_CAPACITY,     SAT_WILT,     HBV_BETA, MAX_CAP_RISE_RATE,  MAX_PERC_RATE,  BASEFLOW_COEFF,            BASEFLOW_N",
                                    "  :Units     ,                    none,          none,         none,         none,              mm/d,           mm/d,             1/d,                  none",
                                    "  #                        HBV_PARA_05,   HBV_PARA_06,  HBV_PARA_14,  HBV_PARA_07,       HBV_PARA_16,       CONSTANT,        CONSTANT,              CONSTANT,",
                                    f"    [DEFAULT],            {params[param_or_name]['HBV']['HBV_Param_05']},  {params[param_or_name]['HBV']['HBV_Param_06']}, {params[param_or_name]['HBV']['HBV_Param_14']}, {params[param_or_name]['HBV']['HBV_Param_07']},      {params[param_or_name]['HBV']['HBV_Param_16']},            0.0,             0.0,                   0.0",
                                    "  #                                                        CONSTANT,                                     HBV_PARA_08,     HBV_PARA_09, 1+HBV_PARA_15=1+ALPHA,",
                                    f"     FAST_RES,                _DEFAULT,      _DEFAULT,          0.0,     _DEFAULT,          _DEFAULT,   {params[param_or_name]['HBV']['HBV_Param_08']},    {params[param_or_name]['HBV']['HBV_Param_09']},              1.877607",
                                    "  #                                                        CONSTANT,                                                      HBV_PARA_10,              CONSTANT,",
                                    f"     SLOW_RES,                _DEFAULT,      _DEFAULT,          0.0,     _DEFAULT,          _DEFAULT,       _DEFAULT,    {params[param_or_name]['HBV']['HBV_Param_10']},                   1.0",
                                    ":EndSoilParameterList"
                                ]
                        },
                    "rvh":
                        {"Subbasins":
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
                                    f"              1,          {params[param_or_name]['HBV']['HBV_Param_11']},                  {params[param_or_name]['HBV']['HBV_Param_11']},",
                                    ":EndSubBasinProperties"
                                ]
                        },
                    "rvi":
                        {"Model Organisation":
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
                        },
                    "rvc":
                        {"Basin":
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
                                    f"{params[param_or_name]['HBV']['HBV_Param_17']}",
                                    ":EndInitialConditions"
                                ]
                        },
                    "rvt":
                        {}
                },
            "MOHYSE":
                {
                    "rvp":
                        {
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
                                    ":SoilProfiles",
                                    "   LAKE, 0",
                                    "   ROCK, 0",
                                    "#  DEFAULT_P,      2, TOPSOIL, MOHYSE_PARA_5, GWSOIL, 10.0",
                                    f"   DEFAULT_P,      2, TOPSOIL,     {params[param_or_name]['MOHYSE']['MOHYSE_Param_05']}, GWSOIL, 10.0",
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
                                    f"# :GlobalParameter     MOHYSE_PET_COEFF  {params[param_or_name]['MOHYSE']['MOHYSE_Param_01']}",
                                    ":GlobalParameter       MOHYSE_PET_COEFF         1.0000"
                                ],
                            "Land Use Parameters":
                                [
                                    ":LandUseParameterList",
                                    "   :Parameters,     MELT_FACTOR,       AET_COEFF, FOREST_SPARSENESS, DD_MELT_TEMP,",
                                    "   :Units,          mm/d/K,            mm/d,                 -,         degC,",
                                    "#      [DEFAULT],   MOHYSE_PARA_3,   MOHYSE_PARA_2,               0.0,MOHYSE_PARA_4, ",
                                    f"      [DEFAULT],          {params[param_or_name]['MOHYSE']['MOHYSE_Param_03']},          {params[param_or_name]['MOHYSE']['MOHYSE_Param_02']},               0.0,       {params[param_or_name]['MOHYSE']['MOHYSE_Param_04']},",
                                    ":EndLandUseParameterList"
                                ],
                            "Soil Parameters":
                                [
                                    ":SoilParameterList",
                                    "   :Parameters,        POROSITY,  PET_CORRECTION,        HBV_BETA,  BASEFLOW_COEFF,      PERC_COEFF, ",
                                    "   :Units,               -,               -,               -,             1/d,             1/d, ",
                                    "#      TOPSOIL,            1.0 ,             1.0,             1.0,   MOHYSE_PARA_7,   MOHYSE_PARA_6,",
                                    "#      GWSOIL,            1.0 ,             1.0,             1.0,   MOHYSE_PARA_8,             0.0,",
                                    f"      TOPSOIL,            1.0 ,             1.0,             1.0,          {params[param_or_name]['MOHYSE']['MOHYSE_Param_07']},          {params[param_or_name]['MOHYSE']['MOHYSE_Param_06']},",
                                    f"      GWSOIL,            1.0 ,             1.0,             1.0,          {params[param_or_name]['MOHYSE']['MOHYSE_Param_08']},             0.0,",
                                    ":EndSoilParameterList"
                                ]
                        },
                    "rvh":
                        {"Subbasins":
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
                                    f"              1,          {params[param_or_name]['MOHYSE']['MOHYSE_Param_10']},                  {params[param_or_name]['MOHYSE']['MOHYSE_Param_09']},",
                                    ":EndSubBasinProperties"
                                ]
                        },
                    "rvi":
                        {"Model Organisation":
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
                        },
                    "rvc":
                        {"Empty":
                            [
                                "# Nothing to set here."
                            ]
                        },
                    "rvt":
                        {}
                }
        }
    logger.debug("rvx_params dictionary created.")
    return rvx_params[model_type]


# TODO: Check the parameter list for those with suffixes 'b'!

# TODO: Which param is HBV_PARA_???


def generate_template_ostrich(model_type: str = model_type,
                              params: dict = default_params,
                              catchment: str = catchment,
                              author: str = author,
                              generation_date: str = generation_date) -> dict:
    """
    Generates template text which can be written to .rvp file

    Args:
        catchment : str
            Catchment name
        model_type : str
            Name of model type
        params : dict
            Dictionary containing parameters values

    Return:
        ost_params[model_type] : dict
            Dictionary with Ostrich parameters
    """

    file_name = f"{catchment}_{model_type}"
    param_or_name = "names"
    logger.debug("Entered generate_template() function.")
    assert model_type in config.variables.supported_models, f"model_type expected GR4J, HYMOD, HMETS, HBV or MOHYSE, got {model_type} instead "
    logger.debug("model_type is in the list of supported models.")
    module_root_dir: Path = Path().resolve()
    ost_params = {
        "GR4J":
            {
                "ost_in":
                    {
                        "Model Info":
                            [
                                f"# Model Type: {model_type}",
                                f"# Catchment Name: {catchment}",
                                f"# Catchment ID: {catchment_id}",
                                f"# Author: {author}",
                                f"# Generation Date: {generation_date}"
                            ],
                        "General Options":
                            [
                                f"ProgramType  	    ParaPADDS",
                                f"ObjectiveFunction   GCOP",
                                f"ModelExecutable     ./Ost-RAVEN.sh",
                                f"PreserveBestModel   ./save_best.sh",
                                f"",
                                f"ModelSubdir processor_",
                                f"",
                                f"# OstrichWarmStart yes",
                                f"# CheckSensitivities yes"
                            ],
                        "Extra Directories":
                            [
                                f"BeginExtraDirs",
                                f"model",
                                f"EndExtraDirs"
                            ],
                        "File Pairs":
                            [
                                f"BeginFilePairs",
                                f"{file_name}.rvp.tpl;	{file_name}.rvp",
                                f"#can be multiple (.rvh, .rvi)",
                                f"EndFilePairs"
                            ],
                        "Parameter Specification":
                            [
                                f"#Parameter/DV Specification",
                                f"#name,initial value, lower bound, upper bound, input, output, internal transformations",
                                f"#name exactly as in *.tpl",
                                f"BeginParams",
                                f"#parameter	   init.	 low		high	tx_in  tx_ost tx_out",
                                f"{params[param_or_name]['GR4J']['GR4J_X1']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['GR4J']['GR4J_X2']}  	random	 	-15		10	none   none 	none",
                                f"{params[param_or_name]['GR4J']['GR4J_X3']}  	random		10		700	none   none 	none",
                                f"{params[param_or_name]['GR4J']['GR4J_X4']}  	random		0		7 	none   none	none",
                                f"{params[param_or_name]['GR4J']['Melt_Factor']}  	random		1		30	none   none	none",
                                f"{params[param_or_name]['GR4J']['Airsnow_Coeff']}  	random		0		1	none   none 	none",
                                f"EndParams"
                            ],
                        "Response Variables":
                            [
                                f"# Reads the Nash-Sutcliffe value from a csv file. Semicolon is a filename separator",
                                f"BeginResponseVars",
                                f"#name	  filename			        keyword		line	col	token                               augmented?",
                                f"NSE      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	3	',' yes",
                                f"KGE      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	4	',' yes",
                                f"EndResponseVars",
                            ],
                        "Tied Response Variables":
                            [
                                f"#Negative Nash-Sutcliffe efficiency",
                                f"BeginTiedRespVars",
                                f"NegNSE 1 NSE wsum -1.00",
                                f"NegKGE 1 KGE wsum -1.00",
                                f"EndTiedRespVars",
                            ],
                        "GCOP Options":
                            [
                                f"BeginGCOP",
                                f"CostFunction NegNSE",
                                f"CostFunction NegKGE",
                                f"PenaltyFunction APM",
                                f"EndGCOP",
                            ],
                        "Constraints":
                            [
                                f"BeginConstraints",
                                f"# not needed when no constraints, but PenaltyFunction statement above is required",
                                f"# name     type     penalty    lwr   upr   resp.var",
                                f"EndConstraints",
                            ],
                        "Random Seed Control":
                            [
                                f"# Randomsed control added",
                                f"RandomSeed 3333",
                            ],
                        "Algorithm Settings":
                            [
                                f"#Algorithm should be last in this file (see p51 for APDDS):",
                                f"",
                                f"BeginParallelPADDSAlg",
                                f"PerturbationValue 0.20",
                                f"MaxIterations 50",
                                f"#	UseRandomParamValues",
                                f"SelectionMetric ExactHyperVolumeContribution",
                                f"EnableDebugging",
                                f"# UseInitialParamValues",
                                f"# Note: above intializes DDS to parameter values IN the initial",
                                f"#       model input files IF 'extract' option used in BeginParams",
                                f"#       block (column 'init')",
                                f"EndParallelPADDSAlg",
                                f"",
                                f"#Begin_PADDSAU_Alg",
                                f"#	PerturbationValue 0.20",
                                f"#    NumSearches 25",
                                f"#    MinItersPerSearch 20",
                                f"#    MaxItersPerSearch 30",
                                f"#    ParallelSearches yes",
                                f"#    Threshold 500",
                                f"#    Randomize no",
                                f"#    ReviseAU yes",
                                f"#	# UseInitialParamValues",
                                f"#	# Note: above intializes DDS to parameter values IN the initial",
                                f"#	#       model input files IF 'extract' option used in BeginParams",
                                f"#	#       block (column 'init')",
                                f"#End_DDSAU_Alg",
                            ]
                    },
                "save_best":
                    {
                        "Save Best":
                            [
                                f"#!/bin/bash{newline}",
                                f"set -e{newline}",
                                f"echo \"saving input files for the best solution found...\"{newline}"
                                f"if [ ! -e model_best ] ; then",
                                f"\tmkdir model_best",
                                f"fi{newline}",
                                f"cp model/{file_name}.rvp                    model_best/{file_name}.rvp",
                                f"cp model/output/{file_name}_Diagnostics.csv model_best/{file_name}_Diagnostics.csv",
                                f"cp model/output/{file_name}_Hydrographs.csv model_best/{file_name}_Hydrographs.csv{newline}",
                                f"exit 0",
                            ]
                    },
                "ost_raven":
                    {
                        "Ost-Raven":
                            [
                                f"#!/bin/bash{newline}",
                                f"set -e{newline}",
                                f"# Get the latest version of the diagnostics script and copy it to the model folder",
                                f"cp {str(module_root_dir)}/raven_tools/processing/raven_diag.py model/output/raven_diag.py{newline}",
                                f"# Copy the latest model files to the model folder",
                                f"cp ./{file_name}.rvp model/{file_name}.rvp{newline}",
                                f"## cd into the model folder",
                                f"cd model{newline}",
                                f"# Run Raven.exe",
                                f"./Raven.exe {file_name} -o output/",
                                f"cd output",
                                f"DIAG_FILE=$(pwd)/{file_name}_Diagnostics.csv",
                                f"HYDROGRAPH_FILE=$(pwd)/{file_name}_Hydrographs.csv",
                                f"source {poetry_location}",
                                f"# shellcheck disable=SC2086",
                                f"python ./raven_diag.py \"$HYDROGRAPH_FILE\" \"$DIAG_FILE\"{newline}",
                                f"exit 0",
                            ]
                    }
            },
        "HYMOD":
            {
                "ost_in":
                    {
                        "Model Info":
                            [
                                f"# Model Type: {model_type}",
                                f"# Catchment Name: {catchment}",
                                f"# Catchment ID: {catchment_id}",
                                f"# Author: {author}",
                                f"# Generation Date: {generation_date}"
                            ],
                        "General Options":
                            [
                                f"ProgramType  	    ParaPADDS",
                                f"ObjectiveFunction   GCOP",
                                f"ModelExecutable     ./Ost-RAVEN.sh",
                                f"PreserveBestModel   ./save_best.sh",
                                f"",
                                f"ModelSubdir processor_",
                                f"",
                                f"# OstrichWarmStart yes",
                                f"# CheckSensitivities yes"
                            ],
                        "Extra Directories":
                            [
                                f"BeginExtraDirs",
                                f"model",
                                f"EndExtraDirs"
                            ],
                        "File Pairs":
                            [
                                f"BeginFilePairs",
                                f"{file_name}.rvh.tpl;	{file_name}.rvh",
                                f"{file_name}.rvp.tpl;	{file_name}.rvp",
                                f"{file_name}.rvi.tpl;	{file_name}.rvi",
                                f"#can be multiple (.rvh, .rvi)",
                                f"EndFilePairs"
                            ],
                        "Parameter Specification":
                            [
                                f"#Parameter/DV Specification",
                                f"#name,initial value, lower bound, upper bound, input, output, internal transformations",
                                f"#name exactly as in *.tpl",
                                f"BeginParams",
                                f"#parameter	   init.	 low		high	tx_in  tx_ost tx_out",
                                f"{params[param_or_name]['HYMOD']['HYMOD_Param_01']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HYMOD']['HYMOD_Param_02']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HYMOD']['HYMOD_Param_03']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HYMOD']['HYMOD_Param_04']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HYMOD']['HYMOD_Param_05']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HYMOD']['HYMOD_Param_06']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HYMOD']['HYMOD_Param_07']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HYMOD']['HYMOD_Param_08']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HYMOD']['HYMOD_Param_09']}		random		0.01		2.5	none   none 	none",
                                f"EndParams"
                            ],
                        "Response Variables":
                            [
                                f"# Reads the Nash-Sutcliffe value from a csv file. Semicolon is a filename separator",
                                f"BeginResponseVars",
                                f"#name	  filename			        keyword		line	col	token                               augmented?",
                                f"NSE      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	3	',' yes",
                                f"KGE      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	4	',' yes",
                                f"EndResponseVars",
                            ],
                        "Tied Response Variables":
                            [
                                f"#Negative Nash-Sutcliffe efficiency",
                                f"BeginTiedRespVars",
                                f"NegNSE 1 NSE wsum -1.00",
                                f"NegKGE 1 KGE wsum -1.00",
                                f"EndTiedRespVars",
                            ],
                        "GCOP Options":
                            [
                                f"BeginGCOP",
                                f"CostFunction NegNSE",
                                f"CostFunction NegKGE",
                                f"PenaltyFunction APM",
                                f"EndGCOP",
                            ],
                        "Constraints":
                            [
                                f"BeginConstraints",
                                f"# not needed when no constraints, but PenaltyFunction statement above is required",
                                f"# name     type     penalty    lwr   upr   resp.var",
                                f"EndConstraints",
                            ],
                        "Random Seed Control":
                            [
                                f"# Randomsed control added",
                                f"RandomSeed 3333",
                            ],
                        "Algorithm Settings":
                            [
                                f"#Algorithm should be last in this file (see p51 for APDDS):",
                                f"",
                                f"BeginParallelPADDSAlg",
                                f"PerturbationValue 0.20",
                                f"MaxIterations 50",
                                f"#	UseRandomParamValues",
                                f"SelectionMetric ExactHyperVolumeContribution",
                                f"EnableDebugging",
                                f"# UseInitialParamValues",
                                f"# Note: above intializes DDS to parameter values IN the initial",
                                f"#       model input files IF 'extract' option used in BeginParams",
                                f"#       block (column 'init')",
                                f"EndParallelPADDSAlg",
                                f"",
                                f"#Begin_PADDSAU_Alg",
                                f"#	PerturbationValue 0.20",
                                f"#    NumSearches 25",
                                f"#    MinItersPerSearch 20",
                                f"#    MaxItersPerSearch 30",
                                f"#    ParallelSearches yes",
                                f"#    Threshold 500",
                                f"#    Randomize no",
                                f"#    ReviseAU yes",
                                f"#	# UseInitialParamValues",
                                f"#	# Note: above intializes DDS to parameter values IN the initial",
                                f"#	#       model input files IF 'extract' option used in BeginParams",
                                f"#	#       block (column 'init')",
                                f"#End_DDSAU_Alg",
                            ]
                    },
                "save_best":
                    {
                        "Save Best":
                            [
                                f"#!/bin/bash{newline}",
                                f"set -e{newline}",
                                f"echo \"saving input files for the best solution found...\"{newline}"
                                f"if [ ! -e model_best ] ; then",
                                f"\tmkdir model_best",
                                f"fi{newline}",
                                f"cp model/{file_name}.rvi                    model_best/{file_name}.rvi",
                                f"cp model/{file_name}.rvh                    model_best/{file_name}.rvh",
                                f"cp model/{file_name}.rvp                    model_best/{file_name}.rvp",
                                f"cp model/output/{file_name}_Diagnostics.csv model_best/{file_name}_Diagnostics.csv",
                                f"cp model/output/{file_name}_Hydrographs.csv model_best/{file_name}_Hydrographs.csv{newline}",
                                f"exit 0",
                            ]
                    },
                "ost_raven":
                    {
                        "Ost-Raven":
                            [
                                f"#!/bin/bash{newline}",
                                f"set -e{newline}",
                                f"# Get the latest version of the diagnostics script and copy it to the model folder",
                                f"cp {str(module_root_dir)}/raven_tools/processing/raven_diag.py model/output/raven_diag.py{newline}",
                                f"# Copy the latest model files to the model folder",
                                f"cp ./{file_name}.rvi model/{file_name}.rvi",
                                f"cp ./{file_name}.rvh model/{file_name}.rvh",
                                f"cp ./{file_name}.rvp model/{file_name}.rvp{newline}",
                                f"## cd into the model folder",
                                f"cd model{newline}",
                                f"# Run Raven.exe",
                                f"./Raven.exe {file_name} -o output/",
                                f"cd output",
                                f"DIAG_FILE=$(pwd)/{file_name}_Diagnostics.csv",
                                f"HYDROGRAPH_FILE=$(pwd)/{file_name}_Hydrographs.csv",
                                f"source {poetry_location}",
                                f"# shellcheck disable=SC2086",
                                f"python ./raven_diag.py \"$HYDROGRAPH_FILE\" \"$DIAG_FILE\"{newline}",
                                f"exit 0",
                            ]
                    }
            },
        "HMETS":
            {
                "ost_in":
                    {
                        "Model Info":
                            [
                                f"# Model Type: {model_type}",
                                f"# Catchment Name: {catchment}",
                                f"# Catchment ID: {catchment_id}",
                                f"# Author: {author}",
                                f"# Generation Date: {generation_date}"
                            ],
                        "General Options":
                            [
                                f"ProgramType  	    ParaPADDS",
                                f"ObjectiveFunction   GCOP",
                                f"ModelExecutable     ./Ost-RAVEN.sh",
                                f"PreserveBestModel   ./save_best.sh",
                                f"",
                                f"ModelSubdir processor_",
                                f"",
                                f"# OstrichWarmStart yes",
                                f"# CheckSensitivities yes"
                            ],
                        "Extra Directories":
                            [
                                f"BeginExtraDirs",
                                f"model",
                                f"EndExtraDirs"
                            ],
                        "File Pairs":
                            [
                                f"BeginFilePairs",
                                f"{file_name}.rvp.tpl;	{file_name}.rvp",
                                f"#can be multiple (.rvh, .rvi)",
                                f"EndFilePairs"
                            ],
                        "Parameter Specification":
                            [
                                f"#Parameter/DV Specification",
                                f"#name,initial value, lower bound, upper bound, input, output, internal transformations",
                                f"#name exactly as in *.tpl",
                                f"BeginParams",
                                f"#parameter	   init.	 low		high	tx_in  tx_ost tx_out",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_01']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_02']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_03']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_04']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_05']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_05b']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_06']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_07']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_08']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_09']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_09b']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_10']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_11']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_12']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_13']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_14']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_15']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_16']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_17']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_18']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_19']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_20']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_20b']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_21']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_21b']}		random		0.01		2.5	none   none 	none",
                                f"EndParams"
                            ],
                        "Tied Parameters":
                            [
                                f"BeginTiedParams",
                                f"# 1-parameter linear (TLIN = 2*XVAL) ",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_05b']} #  1 {params[param_or_name]['HMETS']['HMETS_Param_11']} linear 0.5 0.00 free # SNOW_SWI_MAX ",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_09b']} #  1 {params[param_or_name]['HMETS']['HMETS_Param_11']} linear 0.5 0.00 free # MAX_MELT_FACTOR ",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_20b']} #  1 {params[param_or_name]['HMETS']['HMETS_Param_11']} linear 0.5 0.00 free # :UniformInitialConditions SOIL[0] = HMETS_Param_20/2 ",
                                f"{params[param_or_name]['HMETS']['HMETS_Param_21b']} #  1 {params[param_or_name]['HMETS']['HMETS_Param_11']} linear 0.5 0.00 free # :UniformInitialConditions SOIL[1] = HMETS_Param_21/2 ",
                                f"EndTiedParams"
                            ],
                        "Response Variables":
                            [
                                f"# Reads the Nash-Sutcliffe value from a csv file. Semicolon is a filename separator",
                                f"BeginResponseVars",
                                f"#name	  filename			        keyword		line	col	token                               augmented?",
                                f"NSE      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	3	',' yes",
                                f"KGE      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	4	',' yes",
                                f"EndResponseVars",
                            ],
                        "Tied Response Variables":
                            [
                                f"#Negative Nash-Sutcliffe efficiency",
                                f"BeginTiedRespVars",
                                f"NegNSE 1 NSE wsum -1.00",
                                f"NegKGE 1 KGE wsum -1.00",
                                f"EndTiedRespVars",
                            ],
                        "GCOP Options":
                            [
                                f"BeginGCOP",
                                f"CostFunction NegNSE",
                                f"CostFunction NegKGE",
                                f"PenaltyFunction APM",
                                f"EndGCOP",
                            ],
                        "Constraints":
                            [
                                f"BeginConstraints",
                                f"# not needed when no constraints, but PenaltyFunction statement above is required",
                                f"# name     type     penalty    lwr   upr   resp.var",
                                f"EndConstraints",
                            ],
                        "Random Seed Control":
                            [
                                f"# Randomsed control added",
                                f"RandomSeed 3333",
                            ],
                        "Algorithm Settings":
                            [
                                f"#Algorithm should be last in this file (see p51 for APDDS):",
                                f"",
                                f"BeginParallelPADDSAlg",
                                f"PerturbationValue 0.20",
                                f"MaxIterations 50",
                                f"#	UseRandomParamValues",
                                f"SelectionMetric ExactHyperVolumeContribution",
                                f"EnableDebugging",
                                f"# UseInitialParamValues",
                                f"# Note: above intializes DDS to parameter values IN the initial",
                                f"#       model input files IF 'extract' option used in BeginParams",
                                f"#       block (column 'init')",
                                f"EndParallelPADDSAlg",
                                f"",
                                f"#Begin_PADDSAU_Alg",
                                f"#	PerturbationValue 0.20",
                                f"#    NumSearches 25",
                                f"#    MinItersPerSearch 20",
                                f"#    MaxItersPerSearch 30",
                                f"#    ParallelSearches yes",
                                f"#    Threshold 500",
                                f"#    Randomize no",
                                f"#    ReviseAU yes",
                                f"#	# UseInitialParamValues",
                                f"#	# Note: above intializes DDS to parameter values IN the initial",
                                f"#	#       model input files IF 'extract' option used in BeginParams",
                                f"#	#       block (column 'init')",
                                f"#End_DDSAU_Alg",
                            ]
                    },
                "save_best":
                    {
                        "Save Best":
                            [
                                f"#!/bin/bash{newline}",
                                f"set -e{newline}",
                                f"echo \"saving input files for the best solution found...\"{newline}",
                                f"if [ ! -e model_best ] ; then",
                                f"\tmkdir model_best",
                                f"fi{newline}",
                                f"cp model/{file_name}.rvp                    model_best/{file_name}.rvp",
                                f"cp model/output/{file_name}_Diagnostics.csv model_best/{file_name}_Diagnostics.csv",
                                f"cp model/output/{file_name}_Hydrographs.csv model_best/{file_name}_Hydrographs.csv{newline}",
                                f"exit 0",
                            ]
                    },
                "ost_raven":
                    {
                        "Ost-Raven":

                            [
                                f"#!/bin/bash{newline}",
                                f"set -e{newline}",
                                f"# Get the latest version of the diagnostics script and copy it to the model folder",
                                f"cp {str(module_root_dir)}/raven_tools/processing/raven_diag.py model/output/raven_diag.py{newline}",
                                f"# Copy the latest model files to the model folder",
                                f"cp ./{file_name}.rvc model/{file_name}.rvc",
                                f"cp ./{file_name}.rvp model/{file_name}.rvp{newline}",
                                f"## cd into the model folder",
                                f"cd model{newline}",
                                f"# Run Raven.exe",
                                f"./Raven.exe {file_name} -o output/",
                                f"cd output",
                                f"DIAG_FILE=$(pwd)/{file_name}_Diagnostics.csv",
                                f"HYDROGRAPH_FILE=$(pwd)/{file_name}_Hydrographs.csv",
                                f"source {poetry_location}",
                                f"# shellcheck disable=SC2086",
                                f"python ./raven_diag.py \"$HYDROGRAPH_FILE\" \"$DIAG_FILE\"{newline}",
                                f"exit 0",
                            ]
                    }
            },
        "HBV":
            {
                "ost_in":
                    {
                        "Model Info":
                            [
                                f"# Model Type: {model_type}",
                                f"# Catchment Name: {catchment}",
                                f"# Catchment ID: {catchment_id}",
                                f"# Author: {author}",
                                f"# Generation Date: {generation_date}"
                            ],
                        "General Options":
                            [
                                f"ProgramType  	    ParaPADDS",
                                f"ObjectiveFunction   GCOP",
                                f"ModelExecutable     ./Ost-RAVEN.sh",
                                f"PreserveBestModel   ./save_best.sh",
                                f"",
                                f"ModelSubdir processor_",
                                f"",
                                f"# OstrichWarmStart yes",
                                f"# CheckSensitivities yes"
                            ],
                        "Extra Directories":
                            [
                                f"BeginExtraDirs",
                                f"model",
                                f"EndExtraDirs"
                            ],
                        "File Pairs":
                            [
                                f"BeginFilePairs",
                                f"{file_name}.rvp.tpl;	{file_name}.rvp",
                                f"#can be multiple (.rvh, .rvi)",
                                f"EndFilePairs"
                            ],
                        "Parameter Specification":
                            [
                                f"#Parameter/DV Specification",
                                f"#name,initial value, lower bound, upper bound, input, output, internal transformations",
                                f"#name exactly as in *.tpl",
                                f"BeginParams",
                                f"#parameter	   init.	 low		high	tx_in  tx_ost tx_out",
                                f"{params[param_or_name]['HBV']['HBV_Param_01']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_02']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_03']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_04']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_05']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_06']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_07']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_08']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_09']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_10']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_11']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_12']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_13']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_14']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_15']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_16']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_17']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_18']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_19']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_20']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['HBV']['HBV_Param_21']}		random		0.01		2.5	none   none 	none",

                                f"EndParams"
                            ],
                        "Tied Parameters":
                            [
                                f"BeginTiedParams",
                                f"# 1-parameter linear (TLIN = 2*XVAL) ",
                                f"{params[param_or_name]['HBV']['HBV_Param_11b']} 1 {params[param_or_name]['HBV']['HBV_Param_11b']} linear 0.5 0.00 free",
                                f"{params[param_or_name]['HBV']['HBV_Param_17b']} 1 {params[param_or_name]['HBV']['HBV_Param_17']} linear 0.5 0.00 free",
                                f"EndTiedParams"
                            ],
                        "Response Variables":
                            [
                                f"# Reads the Nash-Sutcliffe value from a csv file. Semicolon is a filename separator",
                                f"BeginResponseVars",
                                f"#name	  filename			        keyword		line	col	token                               augmented?",
                                f"NSE      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	3	',' yes",
                                f"KGE      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	4	',' yes",
                                f"EndResponseVars",
                            ],
                        "Tied Response Variables":
                            [
                                f"#Negative Nash-Sutcliffe efficiency",
                                f"BeginTiedRespVars",
                                f"NegNSE 1 NSE wsum -1.00",
                                f"NegKGE 1 KGE wsum -1.00",
                                f"EndTiedRespVars",
                            ],
                        "GCOP Options":
                            [
                                f"BeginGCOP",
                                f"CostFunction NegNSE",
                                f"CostFunction NegKGE",
                                f"PenaltyFunction APM",
                                f"EndGCOP",
                            ],
                        "Constraints":
                            [
                                f"BeginConstraints",
                                f"# not needed when no constraints, but PenaltyFunction statement above is required",
                                f"# name     type     penalty    lwr   upr   resp.var",
                                f"EndConstraints",
                            ],
                        "Random Seed Control":
                            [
                                f"# Randomsed control added",
                                f"RandomSeed 3333",
                            ],
                        "Algorithm Settings":
                            [
                                f"#Algorithm should be last in this file (see p51 for APDDS):",
                                f"",
                                f"BeginParallelPADDSAlg",
                                f"PerturbationValue 0.20",
                                f"MaxIterations 50",
                                f"#	UseRandomParamValues",
                                f"SelectionMetric ExactHyperVolumeContribution",
                                f"EnableDebugging",
                                f"# UseInitialParamValues",
                                f"# Note: above intializes DDS to parameter values IN the initial",
                                f"#       model input files IF 'extract' option used in BeginParams",
                                f"#       block (column 'init')",
                                f"EndParallelPADDSAlg",
                                f"",
                                f"#Begin_PADDSAU_Alg",
                                f"#	PerturbationValue 0.20",
                                f"#    NumSearches 25",
                                f"#    MinItersPerSearch 20",
                                f"#    MaxItersPerSearch 30",
                                f"#    ParallelSearches yes",
                                f"#    Threshold 500",
                                f"#    Randomize no",
                                f"#    ReviseAU yes",
                                f"#	# UseInitialParamValues",
                                f"#	# Note: above intializes DDS to parameter values IN the initial",
                                f"#	#       model input files IF 'extract' option used in BeginParams",
                                f"#	#       block (column 'init')",
                                f"#End_DDSAU_Alg",
                            ]
                    },
                "save_best":
                    {
                        "Save Best":
                            [
                                f"#!/bin/bash{newline}",
                                f"set -e{newline}",
                                f"echo \"saving input files for the best solution found...\"{newline}"
                                f"if [ ! -e model_best ] ; then",
                                f"\tmkdir model_best",
                                f"fi{newline}",
                                f"cp model/{file_name}.rvi                    model_best/{file_name}.rvi",
                                f"cp model/{file_name}.rvh                    model_best/{file_name}.rvh",
                                f"cp model/{file_name}.rvp                    model_best/{file_name}.rvp",
                                f"cp model/output/{file_name}_Diagnostics.csv model_best/{file_name}_Diagnostics.csv",
                                f"cp model/output/{file_name}_Hydrographs.csv model_best/{file_name}_Hydrographs.csv{newline}",
                                f"exit 0",
                            ]
                    },
                "ost_raven":
                    {
                        "Ost-Raven":
                            [
                                f"#!/bin/bash{newline}",
                                f"set -e{newline}",
                                f"# Get the latest version of the diagnostics script and copy it to the model folder",
                                f"cp {str(module_root_dir)}/raven_tools/processing/raven_diag.py model/output/raven_diag.py{newline}",
                                f"# Copy the latest model files to the model folder",
                                f"cp ./{file_name}.rvi model/{file_name}.rvi",
                                f"cp ./{file_name}.rvh model/{file_name}.rvh",
                                f"cp ./{file_name}.rvc model/{file_name}.rvc",
                                f"cp ./{file_name}.rvp model/{file_name}.rvp{newline}",
                                f"## cd into the model folder",
                                f"cd model{newline}",
                                f"# Run Raven.exe",
                                f"./Raven.exe {file_name} -o output/",
                                f"cd output",
                                f"DIAG_FILE=$(pwd)/{file_name}_Diagnostics.csv",
                                f"HYDROGRAPH_FILE=$(pwd)/{file_name}_Hydrographs.csv",
                                f"source {poetry_location}",
                                f"# shellcheck disable=SC2086",
                                f"python ./raven_diag.py \"$HYDROGRAPH_FILE\" \"$DIAG_FILE\"{newline}",
                                f"exit 0",
                            ]
                    }
            },
        "MOHYSE":
            {
                "ost_in":
                    {
                        "Model Info":
                            [
                                f"# Model Type: {model_type}",
                                f"# Catchment Name: {catchment}",
                                f"# Catchment ID: {catchment_id}",
                                f"# Author: {author}",
                                f"# Generation Date: {generation_date}"
                            ],
                        "General Options":
                            [
                                f"ProgramType  	    ParaPADDS",
                                f"ObjectiveFunction   GCOP",
                                f"ModelExecutable     ./Ost-RAVEN.sh",
                                f"PreserveBestModel   ./save_best.sh",
                                f"",
                                f"ModelSubdir processor_",
                                f"",
                                f"# OstrichWarmStart yes",
                                f"# CheckSensitivities yes"
                            ],
                        "Extra Directories":
                            [
                                f"BeginExtraDirs",
                                f"model",
                                f"EndExtraDirs"
                            ],
                        "File Pairs":
                            [
                                f"BeginFilePairs",
                                f"{file_name}.rvp.tpl;	{file_name}.rvp",
                                f"{file_name}.rvc.tpl;	{file_name}.rvc",
                                f"#can be multiple (.rvh, .rvi)",
                                f"EndFilePairs"
                            ],
                        "Parameter Specification":
                            [
                                f"#Parameter/DV Specification",
                                f"#name,initial value, lower bound, upper bound, input, output, internal transformations",
                                f"#name exactly as in *.tpl",
                                f"BeginParams",
                                f"#parameter	   init.	 low		high	tx_in  tx_ost tx_out",
                                f"{params[param_or_name]['MOHYSE']['MOHYSE_Param_01']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['MOHYSE']['MOHYSE_Param_02']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['MOHYSE']['MOHYSE_Param_03']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['MOHYSE']['MOHYSE_Param_04']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['MOHYSE']['MOHYSE_Param_05']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['MOHYSE']['MOHYSE_Param_06']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['MOHYSE']['MOHYSE_Param_07']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['MOHYSE']['MOHYSE_Param_08']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['MOHYSE']['MOHYSE_Param_09']}		random		0.01		2.5	none   none 	none",
                                f"{params[param_or_name]['MOHYSE']['MOHYSE_Param_10']}		random		0.01		2.5	none   none 	none",
                                f"EndParams"
                            ],
                        "Response Variables":
                            [
                                f"# Reads the Nash-Sutcliffe value from a csv file. Semicolon is a filename separator",
                                f"BeginResponseVars",
                                f"#name	  filename			        keyword		line	col	token                               augmented?",
                                f"NSE      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	3	',' yes",
                                f"KGE      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	4	',' yes",
                                f"EndResponseVars",
                            ],
                        "Tied Response Variables":
                            [
                                f"#Negative Nash-Sutcliffe efficiency",
                                f"BeginTiedRespVars",
                                f"NegNSE 1 NSE wsum -1.00",
                                f"NegKGE 1 KGE wsum -1.00",
                                f"EndTiedRespVars",
                            ],
                        "GCOP Options":
                            [
                                f"BeginGCOP",
                                f"CostFunction NegNSE",
                                f"CostFunction NegKGE",
                                f"PenaltyFunction APM",
                                f"EndGCOP",
                            ],
                        "Constraints":
                            [
                                f"BeginConstraints",
                                f"# not needed when no constraints, but PenaltyFunction statement above is required",
                                f"# name     type     penalty    lwr   upr   resp.var",
                                f"EndConstraints",
                            ],
                        "Random Seed Control":
                            [
                                f"# Randomsed control added",
                                f"RandomSeed 3333",
                            ],
                        "Algorithm Settings":
                            [
                                f"#Algorithm should be last in this file (see p51 for APDDS):",
                                f"",
                                f"BeginParallelPADDSAlg",
                                f"PerturbationValue 0.20",
                                f"MaxIterations 50",
                                f"#	UseRandomParamValues",
                                f"SelectionMetric ExactHyperVolumeContribution",
                                f"EnableDebugging",
                                f"# UseInitialParamValues",
                                f"# Note: above intializes DDS to parameter values IN the initial",
                                f"#       model input files IF 'extract' option used in BeginParams",
                                f"#       block (column 'init')",
                                f"EndParallelPADDSAlg",
                                f"",
                                f"#Begin_PADDSAU_Alg",
                                f"#	PerturbationValue 0.20",
                                f"#    NumSearches 25",
                                f"#    MinItersPerSearch 20",
                                f"#    MaxItersPerSearch 30",
                                f"#    ParallelSearches yes",
                                f"#    Threshold 500",
                                f"#    Randomize no",
                                f"#    ReviseAU yes",
                                f"#	# UseInitialParamValues",
                                f"#	# Note: above intializes DDS to parameter values IN the initial",
                                f"#	#       model input files IF 'extract' option used in BeginParams",
                                f"#	#       block (column 'init')",
                                f"#End_DDSAU_Alg",
                            ]
                    },
                "save_best":
                    {
                        "Save Best":
                            [
                                f"#!/bin/bash{newline}",
                                f"set -e{newline}",
                                f"echo \"saving input files for the best solution found...\"{newline}"
                                f"if [ ! -e model_best ] ; then",
                                f"\tmkdir model_best",
                                f"fi{newline}",
                                f"cp model/{file_name}.rvp                    model_best/{file_name}.rvp",
                                f"cp model/output/{file_name}_Diagnostics.csv model_best/{file_name}_Diagnostics.csv",
                                f"cp model/output/{file_name}_Hydrographs.csv model_best/{file_name}_Hydrographs.csv{newline}",
                                f"exit 0",
                            ]
                    },
                "ost_raven":
                    {
                        "Ost-Raven":
                            [
                                f"#!/bin/bash{newline}",
                                f"set -e{newline}",
                                f"# Get the latest version of the diagnostics script and copy it to the model folder",
                                f"cp {str(module_root_dir)}/raven_tools/processing/raven_diag.py model/output/raven_diag.py{newline}",
                                f"# Copy the latest model files to the model folder",
                                f"cp ./{file_name}.rvi model/{file_name}.rvi",
                                f"cp ./{file_name}.rvh model/{file_name}.rvh",
                                f"cp ./{file_name}.rvp model/{file_name}.rvp{newline}",
                                f"## cd into the model folder",
                                f"cd model{newline}",
                                f"# Run Raven.exe",
                                f"./Raven.exe {file_name} -o output/",
                                f"cd output",
                                f"DIAG_FILE=$(pwd)/{file_name}_Diagnostics.csv",
                                f"HYDROGRAPH_FILE=$(pwd)/{file_name}_Hydrographs.csv",
                                f"source {poetry_location}",
                                f"# shellcheck disable=SC2086",
                                f"python ./raven_diag.py \"$HYDROGRAPH_FILE\" \"$DIAG_FILE\"{newline}",
                                f"exit 0",
                            ]
                    }
            },
    }
    return ost_params[model_type]


# TODO: is 'melt_factor' actually a parameter of GR4J? In the tutorial file, 'annual snow' is called 'CEMANEIGE_X1'...


def subsection_header(title: str) -> list[str]:
    """Generates subsection header for a .rvX file from given title.
    Args:
        title (str):
            Title to be used in subsection header.
    Returns:
        subsection_head: list[str]
            List of header line as strings.
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


def write_rvx(model_dir: str = model_dir,
              model_type: str = model_type,
              data_dir: str = data_dir,
              project_dir: Path = project_dir,
              catchment: str = catchment,
              model_sub_dir: str = model_sub_dir,
              params: dict = default_params,
              template_type: str = "Raven",
              attribute_csv_name: str = "CH-0057_attributes.csv",
              attribute_csv_dir: str = "Hydromap Attributes",
              rvx_type: str = "rvi",
              author=conf['Author']):
    """Writes .rvX file(s), either as an Ostrich or Raven template.
    Args:
        model_dir (str): The directory where the model files are stored. Default is "model_dir".
        model_type (str): The type of model to use. Default is "model_type".
        data_dir (str): The directory where the input data is stored. Default is "data_dir".
        project_dir (Path): The directory of the project. Default is "project_dir".
        catchment (str): The name of the catchment for which the RVX file should be created. Default is "catchment".
        model_sub_dir (str): The sub-directory where the model files are stored. Default is "model_sub_dir".
        params (dict): A dictionary of model parameters. Default is "default_params".
        template_type (str): The type of template to use. Default is "Raven".
        attribute_csv_name (str): The name of the attribute CSV file. Default is "CH-0057_attributes.csv".
        attribute_csv_dir (str): The directory where the attribute CSV file is stored. Default is "Hydromap Attributes".
        rvx_type (str): The type of RVX file to create. Default is "rvi".

        This function writes the RVX file for the given catchment, model type, and parameters.
    """

    logger.debug("Arrived in function write_rvx().")
    assert model_type in config.variables.supported_models, f"Got model type: {model_type}, which is not supported, check variable \"" \
                                                            f"supported_models."
    try:
        logger.debug(
            f"Trying to read catchment attribute CSV file {Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name)}...")
        csv_file = pandas.read_csv(Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name), sep=",",
                                   skiprows=[8],
                                   index_col='attribute_names', usecols=[0, 1])
    except:
        logger.exception()
    logger.debug("Attribute catchment attribute CSV file read.")
    file_name: str = f"{catchment}_{model_type}.{rvx_type}"
    logger.debug(f".{rvx_type} filename set to {file_name}.")
    file_path: Path = Path(project_dir, model_dir, catchment, model_type, model_sub_dir, file_name)
    logger.debug(f".{rvx_type} file path set to {file_path}.")
    template_sections = {}
    logger.debug("Empty dict template_sections created.")

    logger.debug("Entering if-tree for template type...")
    if template_type == "Raven":
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Trying to generate .{rvx_type} template sections with function generate_template()...")
        template_sections = generate_template_rvx(model_type=model_type, csv_file=csv_file, params=params,
                                                  param_or_name="params")
        logger.debug(f"Wrote .{rvx_type} template sections generated by generate_template() to dict template_sections")
    if template_type == "Ostrich":
        file_path: Path = Path(project_dir, model_dir, catchment, model_type, file_name)
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Adding .tpl suffix to file_path for Ostrich template file...")
        file_path = Path((str(file_path) + ".tpl"))
        logger.debug(f"New file_path: {file_path}")
        logger.debug(f"Trying to generate .{rvx_type}.tpl template sections with function generate_template()...")
        template_sections = generate_template_rvx(model_type=model_type, csv_file=csv_file, params=params,
                                                  param_or_name="names")
        logger.debug(
            f"Wrote .{rvx_type}.tpl template sections generated by generate_template() to dict template_sections")

    logger.debug(f"Trying to write to file {file_path}")
    with open(file_path, 'w') as ff:
        ff.writelines(f"{line}{newline}" for line in create_header(author, catchment, model_type))
        logger.debug("Header lines written.")
        ff.write(newline)
        logger.debug("Entering template_sections for-loop...")
        logger.debug(f"template_sections dictionary:\n {template_sections}")
        for section in template_sections[rvx_type]:
            logger.debug(f"Current section: {section}")
            ff.writelines(f"{line}\n" for line in subsection_header(section))
            logger.debug("Subsection header written.")
            ff.writelines(f"{lin}\n" for lin in template_sections[rvx_type][section])
            logger.debug("Template section written")
            ff.write(newline)
        logger.debug("template_sections for-loop finished.")
        logger.debug("Function write_rvx() finished.")


def write_ostrich(
        model_dir: str = model_dir,
        model_type: str = model_type,
        project_dir: Path = project_dir,
        catchment: str = catchment,
        params: dict = default_params,
        ost_in: bool = True,
        save_best: bool = True,
        ost_raven: bool = True
):
    """Writes Ostrich input files ostIn.txt, save_best.sh and Ost-RAVEN.sh
    Args:
        model_dir (str): The directory where the model files are stored. Default is "model_dir".
        model_type (str): The type of model to use. Default is "model_type".
        project_dir (Path): The directory of the project. Default is "project_dir".
        catchment (str): The name of the catchment for which the RVX file should be created. Default is "catchment".
        params (dict): A dictionary of model parameters. Default is "default_params".
        ost_in (bool): Flag to indicate if Ostrich input file should be created. Default is True.
        save_best (bool): Flag to indicate if save_best.sh file should be created. Default is True.
        ost_raven (bool): Flag to indicate if Ost-RAVEN.sh file should be created. Default is True.

    This function writes Ostrich input files ostIn.txt, save_best.sh and Ost-RAVEN.sh
    """
    logger.debug("Arrived in function write_ostrich().")
    assert model_type in config.variables.supported_models, f"Got model type: {model_type}, which is not supported, check variable \"" \
                                                            f"supported_models."

    ost_in_file_name: str = f"ostIn.txt"
    save_best_file_name: str = f"save_best.sh"
    ost_raven_file_name: str = f"Ost-RAVEN.sh"
    logger.debug(f"filename w/o suffix set to {ost_in_file_name}.")
    template_sections = generate_template_ostrich(model_type=model_type, params=params, catchment=catchment)
    logger.debug(f"Dictionary template_sections created.")
    logger.debug(f"Variable ost_in evaluated to {ost_in}")
    if ost_in:
        logger.debug(
            f"Trying to write to file: {Path(project_dir, model_dir, catchment, model_type, ost_in_file_name)}")
        with open(Path(project_dir, model_dir, catchment, model_type, ost_in_file_name), 'w') as ff:
            logger.debug("Entering template_sections for-loop...")
            logger.debug(f"template_sections dictionary: {newline} {template_sections}")
            for section in template_sections["ost_in"]:
                logger.debug(f"Current section: {section}")
                ff.writelines(f"{line}\n" for line in subsection_header(section))
                logger.debug("Subsection header written.")
                ff.writelines(f"{lin}\n" for lin in template_sections["ost_in"][section])
                logger.debug("Template section written.")
                ff.write(newline)
            logger.debug("template_sections for-loop finished.")
    logger.debug(f"Variable save_best evaluated to {save_best}")
    if save_best:
        logger.debug(
            f"Trying to write to file: {Path(project_dir, model_dir, catchment, model_type, save_best_file_name)}")
        with open(Path(project_dir, model_dir, catchment, model_type, save_best_file_name), 'w') as ff:
            logger.debug("Entering template_sections for-loop...")
            logger.debug(f"template_sections dictionary: {newline} {template_sections}")
            for section in template_sections["save_best"]:
                ff.writelines(f"{lin}\n" for lin in template_sections["save_best"][section])
                logger.debug("Template section written.")
                ff.write(newline)
            logger.debug("template_sections for-loop finished.")
    logger.debug(f"Variable ost_raven evaluated to {ost_raven}")
    if ost_raven:
        logger.debug(
            f"Trying to write to file: {Path(project_dir, model_dir, catchment, model_type, ost_raven_file_name)}")
        with open(Path(project_dir, model_dir, catchment, model_type, ost_raven_file_name), 'w') as ff:
            logger.debug("Entering template_sections for-loop...")
            logger.debug(f"template_sections dictionary: {newline} {template_sections}")
            for section in template_sections["ost_raven"]:
                ff.writelines(f"{lin}\n" for lin in template_sections["ost_raven"][section])
                logger.debug("Template section written.")
                ff.write(newline)
            logger.debug("template_sections for-loop finished.")


if __name__ == '__main__':
    # model_dir = Path(project_config['ModelDir'])
    # model_type = Path("raven_broye_gr4j.rvt")
    start_year = 1981
    end_year = 2000
    # write_rvt(1981, 1982)
