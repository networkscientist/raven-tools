"""
Tools to generate Raven .rv* files needed to run Raven models.
"""
import logging
import math
import os
from pathlib import Path

import pandas
import pandas as pd

import raven_tools.config.variables
import raven_tools.processing.raven_preprocess
from raven_tools import config

try:
    logger = logging.getLogger(__name__)
    logger.debug(f"CWD: {os.getcwd()}")
    logger.debug('Trying to read project_config.yaml file')
except:
    logger.exception("Error loading logger.")

try:
    conf = config.variables.project_config
    model_dir = conf['ModelDir']
    model_type = conf['ModelName']
    project_dir = conf['ProjectDir']
    catchment_name = conf['Catchment']
    catchment_id = conf['CatchmentID']
    gauge_lat = conf['GaugeLat']
    gauge_lon = conf['GaugeLon']
    model_sub_dir = conf['ModelSubDir']
    grid_weights_file = conf['GridWeights']
    author = conf['Author']
    generation_date = conf['Date']
    poetry_location = conf['PoetryEnvLocation']
    start_year = conf['StartYear']
    end_year = conf['EndYear']
    cali_end_year = conf['CaliEndYear']
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


def create_header(catchment_ch_id: str, author=conf['Author'], creation_date=generation_date, model=model_type,
                  rvx_type: str = "rvi"):
    """Creates header info for .rvX files

    Args:
        catchment_ch_id: str
            Catchment id
        author: str
            Author name
        creation_date: str
            Date of file generation
        model: str
            Model type, e.g. 'GR4J'
        rvx_type: str
            .rvX file type

    Returns:
        header: list
            List with header info
    """
    try:
        header_line = "#########################################################################"
        file_type = f":FileType          {rvx_type} ASCII Raven 3.5"
        author_line = f":WrittenBy         {author}"
        creation_date = f":CreationDate      {creation_date}"
        description = [
            "#",
            f"# Emulation of {model} simulation of {catchment_ch_id}",
            "#------------------------------------------------------------------------"]
        header = [header_line, file_type, author_line, creation_date, *description]
        return header
    except NameError:
        print("Probably project_config file could not be found...")
        pass


def forcing_block(start_year: int, end_year: int, catchment_ch_id: str, model_type: str):
    """Create Dictionary of forcing data to write in RVT file.

    This function creates a Dictionary of forcing data to be written into an RVT file. From a start and end year,
    it creates the relevant directory Paths, which are based on Swiss gridded input data. Further parameters,
    as required by RAVEN are added and finally, the forcings data block is returned as a Dictionary.

    Args:
        start_year : int
            Start year of forcings data files
        end_year : int
            End year of forcings data files
        catchment_ch_id: str
            Catchment id
        model_type: str
            Model type, e.g. 'GR4J'

    Returns:
        forcing_data: dict[str, list[str]]
            The forcing data block

    """
    if model_type == "HBV":
        grid_weights_file_path = f"data_obs/RhiresD_v2.0_swiss.lv95/out/grid_weights_{catchment_ch_id}_hbv.txt"
    else:
        grid_weights_file_path = f"data_obs/RhiresD_v2.0_swiss.lv95/out/grid_weights_{catchment_ch_id}.txt"
    forcing_rainfall = [
        ":GriddedForcing           Rainfall",
        "    :ForcingType          RAINFALL",
        f"    :FileNameNC           data_obs/RhiresD_v2.0_swiss.lv95/out/RhiresD_v2.0_swiss.lv95_{start_year}01010000_{end_year}12310000_{catchment_ch_id}_clipped.nc",
        "    :VarNameNC            RhiresD",
        "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
        f"    :RedirectToFile       {grid_weights_file_path}",
        ":EndGriddedForcing"]
    forcing_temp_ave = [
        ":GriddedForcing           Average Temperature",
        "    :ForcingType          TEMP_AVE",
        f"    :FileNameNC           data_obs/TabsD_v2.0_swiss.lv95/out/TabsD_v2.0_swiss.lv95_{start_year}01010000_{end_year}12310000_{catchment_ch_id}_clipped.nc",
        "    :VarNameNC            TabsD",
        "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
        f"    :RedirectToFile       {grid_weights_file_path}",
        ":EndGriddedForcing"]
    forcing_temp_max = [
        ":GriddedForcing           Maximum Temperature",
        "    :ForcingType          TEMP_MAX",
        f"    :FileNameNC           data_obs/TmaxD_v2.0_swiss.lv95/out/TmaxD_v2.0_swiss.lv95_{start_year}01010000_{end_year}12310000_{catchment_ch_id}_clipped.nc",
        "    :VarNameNC            TmaxD",
        "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
        f"    :RedirectToFile       {grid_weights_file_path}",
        ":EndGriddedForcing"]
    forcing_temp_min = [
        ":GriddedForcing           Minimum Temperature",
        "    :ForcingType          TEMP_MIN",
        f"    :FileNameNC           data_obs/TminD_v2.0_swiss.lv95/out/TminD_v2.0_swiss.lv95_{start_year}01010000_{end_year}12310000_{catchment_ch_id}_clipped.nc",
        "    :VarNameNC            TminD",
        "    :DimNamesNC           E N time     # must be in the order of (x,y,t) ",
        f"    :RedirectToFile       {grid_weights_file_path}",
        ":EndGriddedForcing"
    ]

    forcing_data = {
        'Rainfall':
            forcing_rainfall,
        'Average Temperature':
            forcing_temp_ave,
        'Maximum Temperature':
            forcing_temp_max,
        'Minimum Temperature':
            forcing_temp_min
    }

    return forcing_data


def write_rvt(start_year: int,
              end_year: int,
              catchment_ch_id: str,
              hru_info: dict,
              model_dir=model_dir,
              model_type=model_type,
              project_dir=project_dir,
              gauge_lat=gauge_lat,
              gauge_lon=gauge_lon,
              model_sub_dir=model_sub_dir,
              author=conf['Author'],
              gauge_short_code="DefAult",
              station_elevation="100",
              catchment_gauge_id="1000",
              params=default_params,
              param_or_name: str = "names",
              template_type: str = "Raven",
              data_dir: Path = Path(conf['ProjectDir'], conf['DataDir'])):
    """Write to Raven .rvt file.

    Args:
        start_year : int
            Start year of forcings data files
        end_year : int
            End year of forcings data files
        catchment_ch_id: str
            Catchment id
        hru_info: dict
            Information on the HRU
        model_dir: str
            Root directory of .rvX files
        model_type: str
            Model type, e.g. 'GR4J'
        project_dir: str
            Root directory of project
        gauge_lat: str
            Latitude of gauge location
        gauge_lon: str
            Longitude of gauge location
        model_sub_dir: str
            Directory name where model files are located, e.g. 'model'
        author: str
            Author name
        gauge_short_code: str
            Short code of gauge
        station_elevation: str
            Elevation of gauge
        catchment_gauge_id: str
            Gauge id
        params: dict
            Default model parameters
        param_or_name: str
            Should parameter values or their names be used?
        template_type: str
            Either 'Ostrich' or 'Raven'
        data_dir: Path
            Data directory path

    """
    import shutil
    if template_type == "Raven":
        param_or_name = "init"
        file_name: str = f"{catchment_ch_id}_{model_type}.rvt"
        logger.debug(f"file_name = {file_name}")
        file_path: Path = Path(model_dir, file_name)
        logger.debug(f"file_path = {file_path}")
    if template_type == "Ostrich":
        param_or_name = "names"
        file_name: str = f"{catchment_ch_id}_{model_type}.rvt.tpl"
        file_path: Path = Path(model_dir, file_name)
        logger.debug(f"file_path = {file_path}")

    if not math.isnan(hru_info['GlaArea']):
        glacier: bool = True
    else:
        glacier: bool = False
    gauge_header = f":Gauge {gauge_short_code}\n"
    gauge_end = f":EndGauge{newline}{newline}"
    gauge_info = [
        f"  :Latitude    {gauge_lat}\n",
        f"  :Longitude {gauge_lon}\n",
        f"  :Elevation  {station_elevation}{newline}{newline}",
    ]

    flow_observation = [
        "# observed streamflow\n",
        f":RedirectToFile data_obs/{gauge_short_code}_Q_{catchment_gauge_id}_daily.rvt"
    ]

    if model_type == "HBV":
        gauge_correction = [
            f"  :RainCorrection    {params['HBV'][param_or_name]['HBV_Param_20']}{newline}",
            f"  :SnowCorrection    {params['HBV'][param_or_name]['HBV_Param_21']}{newline}{newline}"
        ]
        pet_monthly_ave, temp_monthly_ave = raven_tools.processing.raven_preprocess.pet_temp_monthly_ave(
            pet_filepath=Path(project_dir, data_dir, "forcings", "order_103168_PAY_ets150m0_1_data.txt"),
            temp_filepath=Path(project_dir, data_dir, "forcings", "order_103168_PAY_tre200h0_1_data.txt"))
        monthly_averages = [
            # The following line expands the list into a string with spaces between the values, omitting any brackets
            f"  :MonthlyAveEvaporation {' '.join(str(x) for x in pet_monthly_ave.to_list())}{newline}",
            f"  :MonthlyAveTemperature {' '.join(str(x) for x in temp_monthly_ave.to_list())}{newline}{newline}"
        ]
        gauge = [
            gauge_header,
            *gauge_info,
            *gauge_correction,
            *monthly_averages,
            gauge_end
        ]

    else:
        gauge = [
            gauge_header,
            *gauge_info,
            gauge_end
        ]

    with open(file_path, 'w') as ff:
        ff.writelines(f"{line}{newline}" for line in
                      create_header(author=author, catchment_ch_id=catchment_ch_id, model=model_type, rvx_type="rvt"))
        ff.write(f"# meteorological forcings\n")
        for f in forcing_block(start_year, end_year, catchment_ch_id=catchment_ch_id, model_type=model_type).values():
            for t in f:
                ff.write(f"{t}\n")

        ff.writelines(gauge)
        ff.writelines(flow_observation)
    if template_type == "Raven":
        dst_path: Path = Path(model_dir, model_sub_dir, file_name)
        shutil.copy(file_path, dst_path)


def generate_template_rvx(catchment_ch_id: str, hru_info: dict, csv_file=None, model_type=model_type,
                          params=default_params,
                          param_or_name="names",
                          start_year: int = start_year, end_year: int = end_year,
                          cali_end_year: str = cali_end_year,
                          glacier_module: bool = False,
                          data_dir=conf['ModelName']) -> dict:
    """Generates template text which can be written to .rvX file.

        Args:
            catchment_ch_id: str
                Catchment id
            hru_info: dict
                HRU information
            csv_file : str
                File path of the csv file with catchment information
            model_type : str
                Name of the model type
            params : dict
                Default model parameters
            param_or_name : str
                Should parameters values or their names be used?
            start_year: str
                Simulation start year
            end_year: str
                Simulation end year
            cali_end_year: str
                Calibration end year
            glacier_module: bool
                Set True if catchment contains glaciers
            data_dir: str
                Data directory
        Returns:
            rvx_params : dict
                Dictionary containing the parameters for the .rvX file
    """

    assert model_type in config.variables.supported_models, f"model_type expected GR4J, HYMOD, HMETS, HBV or MOHYSE, got {model_type} instead "
    end_date = f"{end_year + 1}-01-01 00:00:00"
    non_glaciated_area: float = float(hru_info['NonGlaArea'])
    non_glacier_altitude: float = float(hru_info['NonGlaAlti'])
    non_glacier_lat: float = float(hru_info['NonGlaLat'])
    non_glacier_lon: float = float(hru_info['NonGlaLon'])
    non_glacier_aspect: float = float(hru_info['NonGlaAspect'])
    non_glacier_slope: float = float(hru_info['NonGlaSlope'])
    if not math.isnan(hru_info['GlaArea']):
        glaciated_area: float = float(hru_info['GlaArea'])
        glacier_altitude: float = float(hru_info['GlaAlti'])
        glacier_lat: float = float(hru_info['GlaLat'])
        glacier_lon: float = float(hru_info['GlaLon'])
        glacier_aspect: float = float(hru_info['GlaAspect'])
        glacier_slope: float = float(hru_info['GlaSlope'])
        new = pd.read_csv(Path(data_dir, "Catchment", "HBV", f"hrus_{catchment_ch_id}.txt"), sep=",", header=None)
        if model_type == "HBV":
            with open(Path(data_dir, "Catchment", "HBV", f"hrus_{catchment_ch_id}.txt"), "r") as f:
                lines = f.readlines()
                lines = [line.strip("\n") for line in lines]
                hru_non_gla_list = []
                hru_list = \
                    [
                        ":HRUs",
                        "  :Attributes,  AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID,LAND_USE_CLASS, VEG_CLASS, SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT",
                        "  :Units     ,   km2,         m,      deg,       deg,     none,          none,      none,         none,            none,          none,   deg,    deg",
                        f"            1, {glaciated_area}, {glacier_altitude}, {glacier_lat}, {glacier_lon}, 1, GLACIER, GLACIER, GLACIER, [NONE], [NONE], {glacier_slope}, {glacier_aspect}",
                        *lines,
                        ":EndHRUs"
                    ]
        else:
            hru_list = \
                [
                    ":HRUs",
                    "  :Attributes,  AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID,LAND_USE_CLASS, VEG_CLASS, SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT",
                    "  :Units     ,   km2,         m,      deg,       deg,     none,          none,      none,         none,            none,          none,   deg,    deg",
                    f"            1, {non_glaciated_area},     {non_glacier_altitude},   {non_glacier_lat},     {non_glacier_lon},        1,        LU_ALL,   VEG_ALL,    DEFAULT_P,          [NONE],        [NONE],   {non_glacier_slope},  {non_glacier_aspect}",
                    f"            2, {glaciated_area}, {glacier_altitude}, {glacier_lat},     {glacier_lon}, 1,        GLACIER,   GLACIER,    GLACIER,          [NONE],        [NONE],   {glacier_slope},  {glacier_aspect}",
                    ":EndHRUs"
                ]
        land_use_classes = \
            [
                ":LandUseClasses",
                "   :Attributes, IMPERM, FOREST_COV",
                "   :Units, frac, frac",
                f"   LU_ALL, {int(csv_file.loc['a0425_clc18_5_aurbv1_0']['values']) / 100}, {int(csv_file.loc['a0418_clc18_5_afrtv1_0']['values']) / 100}",
                "    GLACIER, 0.0, 0.0",
                ":EndLandUseClasses"
            ]
        vegetation_classes = \
            [
                ":VegetationClasses",
                "   :Attributes, MAX_HT, MAX_LAI, MAX_LEAF_COND",
                "   :Units, m, none, mm_per_s",
                "   VEG_ALL, 0.0, 0.0, 0.0",
                "   GLACIER, 0.0, 0.0, 0.0",
                ":EndVegetationClasses"
            ]
    else:
        hru_list = \
            [
                ":HRUs",
                "  :Attributes,  AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID,LAND_USE_CLASS, VEG_CLASS, SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT",
                "  :Units     ,   km2,         m,      deg,       deg,     none,          none,      none,         none,            none,          none,   deg,    deg",
                f"            1, {non_glaciated_area},     {non_glacier_altitude},   {non_glacier_lat},     {non_glacier_lon},        1,        LU_ALL,   VEG_ALL,    DEFAULT_P,          [NONE],        [NONE],   {non_glacier_slope},  {non_glacier_aspect}",
                ":EndHRUs"
            ]
        land_use_classes = \
            [
                ":LandUseClasses",
                "   :Attributes, IMPERM, FOREST_COV",
                "   :Units, frac, frac",
                f"   LU_ALL, {int(csv_file.loc['a0425_clc18_5_aurbv1_0']['values']) / 100}, {int(csv_file.loc['a0418_clc18_5_afrtv1_0']['values']) / 100}",
                ":EndLandUseClasses"
            ]
        vegetation_classes = \
            [
                ":VegetationClasses",
                "   :Attributes, MAX_HT, MAX_LAI, MAX_LEAF_COND",
                "   :Units, m, none, mm_per_s",
                "   VEG_ALL, 0.0, 0.0, 0.0",
                ":EndVegetationClasses"
            ]

    if not glacier_module:
        gr4j_hydro_proc = \
            [
                ":HydrologicProcesses",
                "   :Precipitation            PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE",
                "   :SnowTempEvolve           SNOTEMP_NEWTONS    SNOW_TEMP",
                "   :SnowBalance              SNOBAL_CEMA_NIEGE  SNOW            PONDED_WATER",
                "   :OpenWaterEvaporation     OPEN_WATER_EVAP    PONDED_WATER    ATMOSPHERE     	 # Pn",
                "   :Infiltration             INF_GR4J           PONDED_WATER    MULTIPLE       	 # Ps-",
                "   :SoilEvaporation          SOILEVAP_GR4J      PRODUCT_STORE   ATMOSPHERE     	 # Es",
                "   :Percolation              PERC_GR4J          PRODUCT_STORE   TEMP_STORE     	 # Perc",
                "   :Flush                    RAVEN_DEFAULT      SURFACE_WATER   TEMP_STORE     	 # Pn-Ps",
                "   :Split                    RAVEN_DEFAULT      TEMP_STORE      CONVOLUTION[0] CONVOLUTION[1] 0.9  # Split Pr",
                "   :Convolve                 CONVOL_GR4J_1      CONVOLUTION[0]  ROUTING_STORE  	 # Q9",
                "   :Convolve                 CONVOL_GR4J_2      CONVOLUTION[1]  TEMP_STORE     	 # Q1",
                "   :Percolation              PERC_GR4JEXCH      ROUTING_STORE   GW_STORE       	 # F(x1)",
                "   :Percolation              PERC_GR4JEXCH2     TEMP_STORE      GW_STORE       	 # F(x1)",
                "   :Flush                    RAVEN_DEFAULT      TEMP_STORE      SURFACE_WATER  	 # Qd",
                "   :Baseflow                 BASE_GR4J          ROUTING_STORE   SURFACE_WATER  	 # Qr",
                ":EndHydrologicProcesses"
            ]
    else:
        gr4j_hydro_proc = \
            [
                ":HydrologicProcesses",
                "   :Precipitation            PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE",
                "   :SnowTempEvolve           SNOTEMP_NEWTONS    SNOW_TEMP",
                "   :SnowBalance              SNOBAL_CEMA_NIEGE  SNOW            PONDED_WATER",
                "   :OpenWaterEvaporation     OPEN_WATER_EVAP    PONDED_WATER    ATMOSPHERE     	 # Pn",
                "   :Infiltration             INF_GR4J           PONDED_WATER    MULTIPLE       	 # Ps-",
                "   :Flush                    RAVEN_DEFAULT      PONDED_WATER    GLACIER",
                "     :-->Conditional HRU_TYPE IS GLACIER",
                "   :GlacierMelt              GMELT_SIMPLE_MELT  GLACIER_ICE     GLACIER",
                "   :GlacierRelease           GRELEASE_LINEAR_STORAGE  GLACIER   SURFACE_WATER",
                "   :SoilEvaporation          SOILEVAP_GR4J      PRODUCT_STORE   ATMOSPHERE     	 # Es",
                "   :Percolation              PERC_GR4J          PRODUCT_STORE   TEMP_STORE     	 # Perc",
                "   :Flush                    RAVEN_DEFAULT      SURFACE_WATER   TEMP_STORE     	 # Pn-Ps",
                "     :-->Conditional HRU_TYPE IS_NOT GLACIER",
                "   :Split                    RAVEN_DEFAULT      TEMP_STORE      CONVOLUTION[0] CONVOLUTION[1] 0.9  # Split Pr",
                "   :Convolve                 CONVOL_GR4J_1      CONVOLUTION[0]  ROUTING_STORE  	 # Q9",
                "   :Convolve                 CONVOL_GR4J_2      CONVOLUTION[1]  TEMP_STORE     	 # Q1",
                "   :Percolation              PERC_GR4JEXCH      ROUTING_STORE   GW_STORE       	 # F(x1)",
                "   :Percolation              PERC_GR4JEXCH2     TEMP_STORE      GW_STORE       	 # F(x1)",
                "   :Flush                    RAVEN_DEFAULT      TEMP_STORE      SURFACE_WATER  	 # Qd",
                "   :Baseflow                 BASE_GR4J          ROUTING_STORE   SURFACE_WATER  	 # Qr",
                ":EndHydrologicProcesses"
            ]

    gr4j = {
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
                        f"   GLACIER, 0",
                        f"   DEFAULT_P, 4, SOIL_PROD, {params['GR4J'][param_or_name]['GR4J_X1']}, SOIL_ROUT, 0.300, SOIL_TEMP, 1.000, SOIL_GW, 1.000,",
                        ":EndSoilProfiles"
                    ],
                "Vegetation Classes":
                    vegetation_classes,
                "Land Use Classes":
                    land_use_classes,
                "Global Parameters":
                    [
                        ":GlobalParameter RAINSNOW_TEMP       0.0",
                        ":GlobalParameter RAINSNOW_DELTA      1.0",
                        f":GlobalParameter AIRSNOW_COEFF     {params['GR4J'][param_or_name]['Airsnow_Coeff']} # [1/d] = 1.0 - CEMANEIGE_X2 = 1.0 - x6",
                        f"#CN_X2 = {params['GR4J'][param_or_name]['GR4J_Cemaneige_X2']}",
                        f":GlobalParameter AVG_ANNUAL_SNOW    {params['GR4J'][param_or_name]['Cemaneige_X1']} # [mm]  =       CEMANEIGE_X1 =       x5",
                        "#:GlobalParameter PRECIP_LAPSE     0.0004 I assume not necessary for gridded data",
                        "#:GlobalParameter ADIABATIC_LAPSE  0.0065 not necessary for gridded data"
                    ],
                "Soil Parameters":
                    [
                        ":SoilParameterList",
                        "   :Parameters, POROSITY, GR4J_X3, GR4J_X2",
                        "   :Units, none, mm, mm / d",
                        f"   [DEFAULT], 1.0, {params['GR4J'][param_or_name]['GR4J_X3']}, {params['GR4J'][param_or_name]['GR4J_X2']}",
                        ":EndSoilParameterList"
                    ],
                "Land Use Parameters":
                    [
                        ":LandUseParameterList",
                        "   :Parameters, GR4J_X4, MELT_FACTOR",
                        "   :Units, d, mm / d / C",
                        f"   [DEFAULT], {params['GR4J'][param_or_name]['GR4J_X4']}, 3.5",
                        ":EndLandUseParameterList"
                    ]
            },
        "rvh":
            {"Subbasins":
                [
                    ":SubBasins",
                    "  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,       GAUGED",
                    "  :Units     ,          none,          none,   none,          km,         none",
                    f"            1,        {catchment_ch_id},            -1,   NONE,       _AUTO,     1",
                    ":EndSubBasins"
                ],
                "HRUs":
                    hru_list
            },
        "rvi":
            {"Model Organisation":
                [
                    f":StartDate             {start_year}-01-01 00:00:00",
                    f":EndDate               {end_date}",
                    ":TimeStep              1.0",
                    ":Method                ORDERED_SERIES",
                    f":RunName               {catchment_ch_id}_GR4J"
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
                        ":OroPrecipCorrect      OROCORR_SIMPLELAPSE",
                        f":EvaluationPeriod CALIBRATION {start_year}-01-01 {cali_end_year}-12-31",
                        f":EvaluationPeriod VALIDATION {int(cali_end_year) + 1}-01-01 {end_year}-12-31"
                    ],
                "Soil Layer Alias Definitions":
                    [
                        ":Alias PRODUCT_STORE      SOIL[0]",
                        ":Alias ROUTING_STORE      SOIL[1]",
                        ":Alias TEMP_STORE         SOIL[2]",
                        ":Alias GW_STORE           SOIL[3]"
                    ],
                "Hydrologic Process Order":
                    gr4j_hydro_proc,
                "Output Options":
                    [
                    ]
            },
        "rvc":
            {
                "Soil Profiles":
                    [
                        "# SOIL[0] = GR4J_X1 * 1000. / 2.0 (initialize to 1/2 full)",
                        "# SOIL[1] = 0.3m * 1000. / 2.0   (initialize to 1/2 full)"
                    ],
                "HRU States":
                    [
                        ":HRUStateVariableTable",
                        "   :Attributes SOIL[0] SOIL[1]",
                        "   :Units      mm      mm",
                        f"   1           {params['GR4J'][param_or_name]['GR4J_Soil_0']},   15.0",
                        ":EndHRUStateVariableTable"
                    ]
            },
        "rvt":
            {}
    }
    hbv = {
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
                        "    GLACIER, 0",
                        f"   DEFAULT_P,      3,    TOPSOIL,            {params['HBV'][param_or_name]['HBV_Param_17']},   FAST_RES,    100.0, SLOW_RES,    100.0",
                        ":EndSoilProfiles"
                    ],
                "Vegetation Classes":
                    vegetation_classes,
                "Vegetation Parameters":
                    [
                        ":VegetationParameterList",
                        "   :Parameters,  MAX_CAPACITY, MAX_SNOW_CAPACITY,  TFRAIN,  TFSNOW,",
                        "   :Units,                 mm,                mm,    frac,    frac,",
                        "       VEG_ALL,             10000,             10000,    0.88,    0.88,",
                        ":EndVegetationParameterList"
                    ],
                "Land Use Classes":
                    land_use_classes,
                "Global Parameters":
                    [
                        f":GlobalParameter RAINSNOW_TEMP       {params['HBV'][param_or_name]['HBV_Param_01']}",
                        ":GlobalParameter RAINSNOW_DELTA      1.0 #constant",
                        f":GlobalParameter PRECIP_LAPSE     {params['HBV'][param_or_name]['HBV_Param_12']} # I assume not necessary for gridded data, HBV_PARA_12=PCALT",
                        f":GlobalParameter ADIABATIC_LAPSE  {params['HBV'][param_or_name]['HBV_Param_13']} # not necessary for gridded data, HBV_PARA_13=TCALT",
                        f":GlobalParameter SNOW_SWI  {params['HBV'][param_or_name]['HBV_Param_04']} #HBV_PARA_04"
                    ],
                "Land Use Parameters":
                    [
                        ":LandUseParameterList",
                        "  :Parameters,   MELT_FACTOR, MIN_MELT_FACTOR,   HBV_MELT_FOR_CORR, REFREEZE_FACTOR, HBV_MELT_ASP_CORR",
                        "  :Units     ,        mm/d/K,          mm/d/K,                none,          mm/d/K,              none",
                        "  #              HBV_PARA_02,        CONSTANT,         HBV_PARA_18,     HBV_PARA_03,          CONSTANT",
                        f"    [DEFAULT],  {params['HBV'][param_or_name]['HBV_Param_02']},             2.2,        {params['HBV'][param_or_name]['HBV_Param_18']},    {params['HBV'][param_or_name]['HBV_Param_03']},              0.48",
                        ":EndLandUseParameterList",
                        "",
                        ":LandUseParameterList",
                        " :Parameters, HBV_MELT_GLACIER_CORR,   HBV_GLACIER_KMIN, GLAC_STORAGE_COEFF, HBV_GLACIER_AG",
                        " :Units     ,                  none,                1/d,                1/d,           1/mm",
                        "   #                       CONSTANT,           CONSTANT,        HBV_PARA_19,       CONSTANT,",
                        f"   [DEFAULT],                  1.64,               0.05,       {params['HBV'][param_or_name]['HBV_Param_19']},           0.05",
                        ":EndLandUseParameterList"

                    ],
                "Soil Parameters":
                    [
                        f"#For Ostrich:HBV_Alpha= {params['HBV'][param_or_name]['HBV_Param_15']}",
                        ":SoilParameterList",
                        "  :Parameters,                POROSITY,FIELD_CAPACITY,     SAT_WILT,     HBV_BETA, MAX_CAP_RISE_RATE,  MAX_PERC_RATE,  BASEFLOW_COEFF,            BASEFLOW_N",
                        "  :Units     ,                    none,          none,         none,         none,              mm/d,           mm/d,             1/d,                  none",
                        "  #                        HBV_PARA_05,   HBV_PARA_06,  HBV_PARA_14,  HBV_PARA_07,       HBV_PARA_16,       CONSTANT,        CONSTANT,              CONSTANT,",
                        f"    [DEFAULT],            {params['HBV'][param_or_name]['HBV_Param_05']},  {params['HBV'][param_or_name]['HBV_Param_06']}, {params['HBV'][param_or_name]['HBV_Param_14']}, {params['HBV'][param_or_name]['HBV_Param_07']},      {params['HBV'][param_or_name]['HBV_Param_16']},            0.0,             0.0,                   0.0",
                        "  #                                                        CONSTANT,                                     HBV_PARA_08,     HBV_PARA_09, 1+HBV_PARA_15=1+ALPHA,",
                        f"     FAST_RES,                _DEFAULT,      _DEFAULT,          0.0,     _DEFAULT,          _DEFAULT,   {params['HBV'][param_or_name]['HBV_Param_08']},    {params['HBV'][param_or_name]['HBV_Param_09']},              {params['HBV'][param_or_name]['HBV_Param_15b']}",
                        "  #                                                        CONSTANT,                                                      HBV_PARA_10,              CONSTANT,",
                        f"     SLOW_RES,                _DEFAULT,      _DEFAULT,          0.0,     _DEFAULT,          _DEFAULT,       _DEFAULT,    {params['HBV'][param_or_name]['HBV_Param_10']},                   1.0",
                        ":EndSoilParameterList"
                    ]
            },
        "rvh":
            {"Subbasins":
                [
                    ":SubBasins",
                    "  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,       GAUGED",
                    "  :Units     ,          none,          none,   none,          km,         none",
                    f"            1,        {catchment_ch_id},            -1,   NONE,       _AUTO,     1",
                    ":EndSubBasins"
                ],
                "HRUs":
                    hru_list,
                "Subbasin Properties":
                    [
                        ":SubBasinProperties",
                        "#                       HBV_PARA_11, DERIVED FROM HBV_PARA_11,",
                        "#                            MAXBAS,                 MAXBAS/2,",
                        "   :Parameters,           TIME_CONC,             TIME_TO_PEAK,",
                        "   :Units,                        d,                        d,",
                        f"              1,          {params['HBV'][param_or_name]['HBV_Param_11']},                  {params['HBV'][param_or_name]['HBV_Param_11b']},",
                        ":EndSubBasinProperties"
                    ]
            },
        "rvi":
            {"Model Organisation":
                [
                    f":StartDate             {start_year}-01-01 00:00:00",
                    f":EndDate               {end_date}",
                    ":TimeStep              1.0",
                    f":RunName               {catchment_ch_id}_HBV"
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
                        ":SoilModel                  SOIL_MULTILAYER 3",
                        f":EvaluationPeriod   CALIBRATION   {start_year}-01-01   {cali_end_year}-12-31",
                        f":EvaluationPeriod   VALIDATION    {int(cali_end_year) + 1}-01-01   {end_year}-12-31"
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
                        "   :CanopySnowEvap    CANEVP_ALL         CANOPY_SNOW     ATMOSPHERE",
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
                        f"{params['HBV'][param_or_name]['HBV_Param_17b']}",
                        ":EndInitialConditions"
                    ]
            },
        "rvt":
            {}
    }
    hmets = {
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
                        "   GLACIER, 0",
                        "   # DEFAULT_P, 2, TOPSOIL,          x(20)/1000, PHREATIC,         x(21)/1000,",
                        f"  DEFAULT_P, 2, TOPSOIL,     {params['HMETS'][param_or_name]['HMETS_Param_20b']}, PHREATIC,     {params['HMETS'][param_or_name]['HMETS_Param_21b']},",
                        ":EndSoilProfiles"
                    ],
                "Vegetation Classes":
                    vegetation_classes,
                "Land Use Classes":
                    land_use_classes,
                "Global Parameters":
                    [
                        f":GlobalParameter  SNOW_SWI_MIN {params['HMETS'][param_or_name]['HMETS_Param_09a']} # x(9)",
                        f":GlobalParameter  SNOW_SWI_MAX {params['HMETS'][param_or_name]['HMETS_Param_09b']} # x(9)+x(10) = {params['HMETS'][param_or_name]['HMETS_Param_09a']} + {params['HMETS'][param_or_name]['HMETS_Param_10']}",
                        f":GlobalParameter  SWI_REDUCT_COEFF {params['HMETS'][param_or_name]['HMETS_Param_11']} # x(11)",
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
                        f"      [DEFAULT],  {params['HMETS'][param_or_name]['HMETS_Param_05a']}, {params['HMETS'][param_or_name]['HMETS_Param_05b']},  {params['HMETS'][param_or_name]['HMETS_Param_07']},  {params['HMETS'][param_or_name]['HMETS_Param_08']},  {params['HMETS'][param_or_name]['HMETS_Param_13']},  {params['HMETS'][param_or_name]['HMETS_Param_14']},   {params['HMETS'][param_or_name]['HMETS_Param_12']},     {params['HMETS'][param_or_name]['HMETS_Param_16']},",
                        f"#      x(5),       x(5)+x(6) = {params['HMETS'][param_or_name]['HMETS_Param_05a']} + {params['HMETS'][param_or_name]['HMETS_Param_06']},            x(7),            x(8),           x(13),           x(14),            x(12),              x(16),",
                        ":EndLandUseParameterList",
                        "",
                        ":LandUseParameterList",
                        "   :Parameters,     GAMMA_SHAPE,     GAMMA_SCALE,    GAMMA_SHAPE2,    GAMMA_SCALE2,",
                        "   :Units,               -,             1/d,               -,             1/d,",
                        f"      [DEFAULT],  {params['HMETS'][param_or_name]['HMETS_Param_01']},  {params['HMETS'][param_or_name]['HMETS_Param_02']},  {params['HMETS'][param_or_name]['HMETS_Param_03']},  {params['HMETS'][param_or_name]['HMETS_Param_04']},",
                        "#      x(1),            x(2),            x(3),            x(4),",
                        ":EndLandUseParameterList"
                    ],
                "Soil Parameters":
                    [
                        ":SoilParameterList",
                        "   :Parameters,        POROSITY,      PERC_COEFF,  PET_CORRECTION, BASEFLOW_COEFF",
                        "   :Units,               -,             1/d,               -,            1/d",
                        f"      TOPSOIL,             1.0,  {params['HMETS'][param_or_name]['HMETS_Param_17']},  {params['HMETS'][param_or_name]['HMETS_Param_15']}, {params['HMETS'][param_or_name]['HMETS_Param_18']}",
                        f"      PHREATIC,             1.0,             0.0,             0.0, {params['HMETS'][param_or_name]['HMETS_Param_19']}",
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
                    f"            1,        {catchment_ch_id},            -1,   NONE,       _AUTO,     1",
                    ":EndSubBasins"
                ],
                "HRUs":
                    hru_list
            },
        "rvi":
            {"Model Organisation":
                [
                    f":StartDate               {start_year}-01-01 00:00:00",
                    f":EndDate                {end_date}",
                    ":TimeStep                1.0",
                    ":Method                  ORDERED_SERIES",
                    f":RunName                 {catchment_ch_id}_HMETS"
                ],
                "Model Options":
                    [
                        ":PotentialMeltMethod     POTMELT_HMETS",
                        ":RainSnowFraction        RAINSNOW_DATA",
                        "#:Evaporation             PET_DATA",
                        ":Evaporation            PET_OUDIN",
                        ":CatchmentRoute          ROUTE_DUMP",
                        ":Routing                 ROUTE_NONE",
                        ":SoilModel               SOIL_TWO_LAYER",
                        f":EvaluationPeriod   CALIBRATION   {start_year}-01-01   {cali_end_year}-12-31",
                        f":EvaluationPeriod   VALIDATION    {int(cali_end_year) + 1}-01-01   {end_year}-12-31"
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
                    ]
            },
        "rvc":
            {"Initial Storage":
                [
                    "# initialize to 1/2 full",
                    "# x(20b)/2",
                    f"#:UniformInitialConditions SOIL[0] {params['HMETS'][param_or_name]['HMETS_Param_20b']}",
                    "# x(21b)/2",
                    f"#:UniformInitialConditions SOIL[1] {params['HMETS'][param_or_name]['HMETS_Param_21b']}"
                ],
                "HRUs":
                    [
                        ":HRUStateVariableTable # (according to rchlumsk-BMSC-cf9a83c, modelname.rvc.tpl: formerly :InitialConditionsTable)",
                        ":Attributes SOIL[0] SOIL[1]",
                        ":Units mm mm",
                        f"1 {params['HMETS'][param_or_name]['HMETS_Param_20a']} {params['HMETS'][param_or_name]['HMETS_Param_21a']}",
                        ":EndHRUStateVariableTable",
                    ]
            },
        "rvt":
            {}
    }
    hymod = {
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
                        "   GLACIER, 0,",
                        "   # DEFAULT_P,      2, TOPSOIL,  HYMOD_PARA_2, GWSOIL, 10.0",
                        f"   DEFAULT_P, 2, TOPSOIL, {params['HYMOD'][param_or_name]['HYMOD_Param_02']}, GWSOIL, 10.0",
                        ":EndSoilProfiles"
                    ],
                "Land Use Classes":
                    land_use_classes,
                "Vegetation Classes":
                    vegetation_classes,
                "Global Parameters":
                    [
                        f":GlobalParameter RAINSNOW_TEMP {params['HYMOD'][param_or_name]['HYMOD_Param_03']}",
                        "   #:GlobalParameter      RAINSNOW_TEMP    HYMOD_PARA_3"
                    ],
                "Soil Parameters":
                    [
                        ":SoilParameterList",
                        "   :Parameters, POROSITY, PET_CORRECTION, BASEFLOW_COEFF,",
                        "   :Units, -, -, 1 / d,",
                        "       # TOPSOIL,            1.0 ,    HYMOD_PARA_8,               0.0,",
                        "       #  GWSOIL,            1.0 ,             1.0,   HYMOD_PARA_4=Ks,",
                        f"       TOPSOIL, 1.0, {params['HYMOD'][param_or_name]['HYMOD_Param_08']}, 0.0,",
                        f"       GWSOIL, 1.0, 1.0, {params['HYMOD'][param_or_name]['HYMOD_Param_04']},",
                        ":EndSoilParameterList"
                    ],
                "Land Use Parameters":
                    [
                        ":LandUseParameterList",
                        "   :Parameters, MELT_FACTOR, DD_MELT_TEMP, PDM_B,",
                        "   :Units, mm / d / K, degC, -,",
                        "       # [DEFAULT],    HYMOD_PARA_5,    HYMOD_PARA_6,  HYMOD_PARA_7=Bexp,",
                        f"       [DEFAULT], {params['HYMOD'][param_or_name]['HYMOD_Param_05']}, {params['HYMOD'][param_or_name]['HYMOD_Param_06']}, {params['HYMOD'][param_or_name]['HYMOD_Param_07']},",
                        ":EndLandUseParameterList"
                    ]
            },
        "rvh":
            {"Subbasins":
                [
                    ":SubBasins",
                    "  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,       GAUGED",
                    "  :Units     ,          none,          none,   none,          km,         none",
                    f"            1,        {catchment_ch_id},            -1,   NONE,       _AUTO,     1",
                    ":EndSubBasins"
                ],
                "HRUs":
                    hru_list,
                "Subbasin Properties":
                    [
                        ":SubBasinProperties",
                        "#                         HYMOD_PARA_1,                  3,",
                        "   :Parameters,           RES_CONSTANT,     NUM_RESERVOIRS,",
                        "   :Units,                         1/d,                  -,",
                        f"              1,          {params['HYMOD'][param_or_name]['HYMOD_Param_01']},                  3,",
                        ":EndSubBasinProperties"
                    ]
            },
        "rvi":
            {"Model Organisation":
                [
                    f":StartDate          {start_year}-01-01 00:00:00",
                    f":EndDate            {end_date}",
                    ":TimeStep           1.0",
                    ":Method             ORDERED_SERIES",
                    f":RunName            {catchment_ch_id}_HYMOD"
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
                        f":EvaluationPeriod   CALIBRATION   {start_year}-01-01   {cali_end_year}-12-31",
                        f":EvaluationPeriod   VALIDATION    {int(cali_end_year) + 1}-01-01   {end_year}-12-31"
                    ],
                "Hydrologic Process Order":
                    [
                        ":HydrologicProcesses",
                        "   :Precipitation     PRECIP_RAVEN       ATMOS_PRECIP    MULTIPLE",
                        "   :SnowBalance       SNOBAL_SIMPLE_MELT SNOW            PONDED_WATER",
                        "   :Infiltration      INF_PDM            PONDED_WATER    MULTIPLE",
                        "#  :Flush            RAVEN_DEFAULT      SURFACE_WATER   SOIL[1]   HYMOD_PARAM_9=ALPHA",
                        f"  :Flush             RAVEN_DEFAULT      SURFACE_WATER   SOIL[1]          {params['HYMOD'][param_or_name]['HYMOD_Param_09']}",
                        "   :SoilEvaporation   SOILEVAP_PDM       SOIL[0]         ATMOSPHERE",
                        "   :Baseflow          BASE_LINEAR        SOIL[1]         SURFACE_WATER",
                        ":EndHydrologicProcesses"
                    ],
                "Output Options":
                    [
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
    mohyse = {
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
                        "   GLACIER, 0",
                        "#  DEFAULT_P,      2, TOPSOIL, MOHYSE_PARA_5, GWSOIL, 10.0",
                        f"   DEFAULT_P,      2, TOPSOIL,     {params['MOHYSE'][param_or_name]['MOHYSE_Param_05']}, GWSOIL, 10.0",
                        ":EndSoilProfiles"
                    ],
                "Vegetation Classes":
                    vegetation_classes,
                "Vegetation Parameters":
                    [
                        ":VegetationParameterList",
                        "   :Parameters,    SAI_HT_RATIO,  RAIN_ICEPT_PCT,  SNOW_ICEPT_PCT,",
                        "   :Units,               -,               -,               -, ",
                        "       [DEFAULT],             0.0,             0.0,             0.0,   ",
                        ":EndVegetationParameterList"
                    ],
                "Land Use Classes":
                    land_use_classes,
                "Global Parameters":
                    [
                        "#:GlobalParameter      RAINSNOW_TEMP              -2.0",
                        ":GlobalParameter       TOC_MULTIPLIER              1.0",
                        f"# :GlobalParameter     MOHYSE_PET_COEFF  MOHYSE_PARA_01",
                        f":GlobalParameter       MOHYSE_PET_COEFF         {params['MOHYSE'][param_or_name]['MOHYSE_Param_01']}"
                    ],
                "Land Use Parameters":
                    [
                        ":LandUseParameterList",
                        "   :Parameters,     MELT_FACTOR,       AET_COEFF, FOREST_SPARSENESS, DD_MELT_TEMP,",
                        "   :Units,          mm/d/K,            mm/d,                 -,         degC,",
                        "#      [DEFAULT],   MOHYSE_PARA_3,   MOHYSE_PARA_2,               0.0,MOHYSE_PARA_4, ",
                        f"      [DEFAULT],          {params['MOHYSE'][param_or_name]['MOHYSE_Param_03']},          {params['MOHYSE'][param_or_name]['MOHYSE_Param_02']},               0.0,       {params['MOHYSE'][param_or_name]['MOHYSE_Param_04']},",
                        ":EndLandUseParameterList"
                    ],
                "Soil Parameters":
                    [
                        ":SoilParameterList",
                        "   :Parameters,        POROSITY,  PET_CORRECTION,        HBV_BETA,  BASEFLOW_COEFF,      PERC_COEFF, ",
                        "   :Units,               -,               -,               -,             1/d,             1/d, ",
                        "#      TOPSOIL,            1.0 ,             1.0,             1.0,   MOHYSE_PARA_7,   MOHYSE_PARA_6,",
                        "#      GWSOIL,            1.0 ,             1.0,             1.0,   MOHYSE_PARA_8,             0.0,",
                        f"      TOPSOIL,            1.0 ,             1.0,             1.0,          {params['MOHYSE'][param_or_name]['MOHYSE_Param_07']},          {params['MOHYSE'][param_or_name]['MOHYSE_Param_06']},",
                        f"      GWSOIL,            1.0 ,             1.0,             1.0,          {params['MOHYSE'][param_or_name]['MOHYSE_Param_08']},             0.0,",
                        ":EndSoilParameterList"
                    ]
            },
        "rvh":
            {"Subbasins":
                [
                    ":SubBasins",
                    "  :Attributes,          NAME, DOWNSTREAM_ID,PROFILE,REACH_LENGTH,       GAUGED",
                    "  :Units     ,          none,          none,   none,          km,         none",
                    f"            1,        {catchment_ch_id},            -1,   NONE,       _AUTO,     1",
                    ":EndSubBasins"
                ],
                "HRUs":
                    hru_list,
                "Subbasin Properties":
                    [
                        ":SubBasinProperties",
                        "#          1.0 / MOHYSE_PARA_10,   MOHYSE_PARA_9",
                        "   :Parameters,     GAMMA_SCALE,     GAMMA_SHAPE,",
                        "   :Units,                  1/d,               -",
                        f"              1,          {params['MOHYSE'][param_or_name]['MOHYSE_Param_10']},                  {params['MOHYSE'][param_or_name]['MOHYSE_Param_09']},",
                        ":EndSubBasinProperties"
                    ]
            },
        "rvi":
            {"Model Organisation":
                [
                    f":StartDate               {start_year}-01-01 00:00:00",
                    f":EndDate                {end_date}",
                    ":TimeStep                1.0",
                    ":Method                  ORDERED_SERIES",
                    f":RunName                 {catchment_ch_id}_MOHYSE"
                ],
                "Model Options":
                    [
                        ":SoilModel             SOIL_TWO_LAYER",
                        ":PotentialMeltMethod   POTMELT_DEGREE_DAY",
                        ":Routing               ROUTE_NONE",
                        ":CatchmentRoute        ROUTE_GAMMA_CONVOLUTION",
                        ":Evaporation           PET_MOHYSE",
                        ":DirectEvaporation",
                        ":RainSnowFraction      RAINSNOW_DATA",
                        f":EvaluationPeriod   CALIBRATION   {start_year}-01-01   {cali_end_year}-12-31",
                        f":EvaluationPeriod   VALIDATION    {int(cali_end_year) + 1}-01-01   {end_year}-12-31"
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
    rvx_params = {
        "GR4J": gr4j,
        "HBV": hbv,
        "HMETS": hmets,
        "HYMOD": hymod,
        "MOHYSE": mohyse
    }
    logger.debug("rvx_params dictionary created.")
    return rvx_params[model_type]


def generate_template_ostrich(catchment_ch_id: str,
                              model_type: str = model_type,
                              params: dict = default_params,
                              catchment_name: str = catchment_name,
                              author: str = author,
                              generation_date: str = generation_date,
                              max_iterations: int = 20) -> dict:
    """
    Generates template text which can be written to .rvp file

    Args:
        catchment_ch_id: str
            Catchment id
        model_type : str
            Name of model type, e.g. 'GR4J'
        params : dict
            Dictionary containing parameters values
        catchment_name : str
            Catchment name
        author: str
            Author name
        generation_date: str
            File generation date
        max_iterations: int
            Maximum number of Ostrich iteration runs.

    Return:
        ost_params[model_type] : dict
            Dictionary with Ostrich parameters
    """

    file_name = f"{catchment_ch_id}_{model_type}"
    assert model_type in config.variables.supported_models, f"model_type expected GR4J, HYMOD, HMETS, HBV or MOHYSE, got {model_type} instead "
    module_root_dir: Path = Path().resolve()
    response_variables = [
        f"# Reads the Nash-Sutcliffe value from a csv file. Semicolon is a filename separator",
        f"BeginResponseVars",
        f"#name	  filename			        keyword		line	col	token                               augmented?",
        f"KGE_NP_CALI      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	2	',' yes",
        f"PBIAS_CALI      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	3	',' yes",
        f"RMSE_CALI      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	4	',' yes",
        f"VE_CALI           ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	5	',' yes",
        f"KGE_NP_Cost      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	6	',' no",
        f"PBIAS_Cost      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_CALIBRATION	0	7	',' no",
        f"KGE_NP_VALI      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_VALIDATION	0	2	',' yes",
        f"PBIAS_VALI      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_VALIDATION	0	3	',' yes",
        f"RMSE_VALI      ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_VALIDATION	0	4	',' yes",
        f"VE_VALI           ./model/output/{file_name}_Diagnostics.csv;	HYDROGRAPH_VALIDATION	0	5	',' yes",
        f"EndResponseVars"
    ]
    tied_response_variable_mohyse = "#"
    if model_type == "MOHYSE":
        # tied_response_variable_mohyse = f"#{params['MOHYSE']['names']['MOHYSE_Param_08b']}   3 {params['MOHYSE']['names']['MOHYSE_Param_06']} {params['MOHYSE']['names']['MOHYSE_Param_07']} {params['MOHYSE']['names']['MOHYSE_Param_08']} wsum 1.0 1.0 1.0"
        pass
    else:
        tied_response_variable_mohyse = "#"
    tied_response_variables = [
        f"#Negative Non-Parametric Kling-Gupta efficiency",
        f"#BeginTiedRespVars",
        f"#KGE_NP_Cost 1 KGE_NP linear -1.00 -1.00",
        tied_response_variable_mohyse,
        f"#EndTiedRespVars",
    ]

    gcop_options = [
        f"BeginGCOP",
        f"CostFunction KGE_NP_Cost",
        f"PenaltyFunction APM",
        f"EndGCOP",
    ]
    algorithm_settings = [
        f"BeginParallelDDSAlg",
        f"PerturbationValue 0.20",
        f"MaxIterations {max_iterations}",
        f"UseRandomParamValues",
        f"# UseInitialParamValues",
        f"EndParallelDDSAlg"
    ]
    random_seed = [
        f"# Randomsed control added",
        f"#RandomSeed 3333",
    ]
    general_options = [
        f"ProgramType  	    ParallelDDS",
        f"ObjectiveFunction   GCOP",
        f"ModelExecutable     ./Ost-RAVEN.sh",
        f"PreserveBestModel   ./save_best.sh",
        f"",
        f"ModelSubdir processor_",
        f"",
        f"# OstrichWarmStart yes",
    ]

    #    ost_raven = [
    #        f"cp ./{file_name}.rvc model/{file_name}.rvc",
    #        f"cp ./{file_name}.rvh model/{file_name}.rvh",
    #        f"cp ./{file_name}.rvi model/{file_name}.rvi",
    #        f"cp ./{file_name}.rvp model/{file_name}.rvp",
    #        f"cp ./{file_name}.rvt model/{file_name}.rvt{newline}",
    #        f"## cd into the model folder",
    #        f"cd model{newline}",
    #        f"# Run Raven.exe",
    #        f"./Raven.exe {file_name} -o output/",
    #        f"cd output",
    #        f"DIAG_FILE=$(pwd)/{file_name}_Diagnostics.csv",
    #        f"HYDROGRAPH_FILE=$(pwd)/{file_name}_Hydrographs.csv",
    #        f"source {poetry_location}",
    #        f"python ./raven_diag.py \"$HYDROGRAPH_FILE\" \"$DIAG_FILE\"{newline}",
    #        f"exit 0",
    #    ]

    ost_raven_script_cp_lines = []
    if model_type in ["GR4J", "HMETS"]:
        ost_raven_script_cp_lines = [
            f"cp ./{file_name}.rvc model/{file_name}.rvc",
            f"cp ./{file_name}.rvp model/{file_name}.rvp",
        ]
    elif model_type == "HBV":
        ost_raven_script_cp_lines = [
            f"cp ./{file_name}.rvc model/{file_name}.rvc",
            f"cp ./{file_name}.rvh model/{file_name}.rvh",
            f"cp ./{file_name}.rvp model/{file_name}.rvp",
            f"cp ./{file_name}.rvt model/{file_name}.rvt",
        ]
    elif model_type == "MOHYSE":
        ost_raven_script_cp_lines = [
            f"cp ./{file_name}.rvh model/{file_name}.rvh",
            f"cp ./{file_name}.rvp model/{file_name}.rvp",
        ]

    ost_raven_header = [
        f"#!/bin/bash{newline}",
        f"set -e{newline}",
        f"# Get the latest version of the diagnostics script and copy it to the model folder",
        f"cp /storage/homefs/pz09y074/raven_master_files/raven_tools/raven_tools/processing/raven_diag.py model/output/raven_diag.py{newline}",
        f"# Copy the latest model files to the model folder",
    ]

    ost_raven_footer = [
        f"## cd into the model folder",
        f"cd model{newline}",
        f"# Run Raven.exe",
        f"./Raven.exe {file_name} -o output/",
        f"cd output",
        f"DIAG_FILE=$(pwd)/{file_name}_Diagnostics.csv",
        f"HYDROGRAPH_FILE=$(pwd)/{file_name}_Hydrographs.csv",
        f"source {poetry_location}",
        f"python ./raven_diag.py \"$HYDROGRAPH_FILE\" \"$DIAG_FILE\"{newline}",
        f"exit 0",
    ]

    ost_raven = ost_raven_header + ost_raven_script_cp_lines + ost_raven_footer

    ost_mpi_script_cp_lines = []
    if model_type in ["GR4J", "HMETS"]:
        ost_mpi_script_cp_lines = [
            f"cp ./{file_name}.rvp model/{file_name}.rvp",
            f"cp ./{file_name}.rvc model/{file_name}.rvc"
        ]
    elif model_type == "HBV":
        ost_mpi_script_cp_lines = [
            f"cp ./{file_name}.rvc model/{file_name}.rvc",
            f"cp ./{file_name}.rvh model/{file_name}.rvh",
            f"cp ./{file_name}.rvp model/{file_name}.rvp",
            f"cp ./{file_name}.rvt model/{file_name}.rvt"
        ]
    elif model_type == "MOHYSE":
        ost_mpi_script_cp_lines = [
            f"cp ./{file_name}.rvh model/{file_name}.rvh",
            f"cp ./{file_name}.rvp model/{file_name}.rvp"
        ]
    ost_mpi_header = [
        f"#!/bin/bash{newline}{newline}",
        f"# match assignment to location of OSTRICH installation{newline}",
        f"set -e{newline}"
    ]
    ost_mpi_footer = [
        f"OSTRICH_MPI=./OstrichMPI{newline}{newline}",
        f"mpirun $OSTRICH_MPI{newline}"
    ]
    ost_mpi_script = {
        "Ostrich MPI run":
            ost_mpi_header + ost_mpi_script_cp_lines + ost_mpi_footer
    }

    #    ost_mpi_script = {
    #        "Ostrich MPI run":
    #            [
    #                f"#!/bin/bash{newline}{newline}",
    #                f"# match assignment to location of OSTRICH installation{newline}",
    #                f"cp ./{file_name}.rvi model/{file_name}.rvi",
    #                f"cp ./{file_name}.rvh model/{file_name}.rvh",
    #                f"cp ./{file_name}.rvt model/{file_name}.rvt",
    #                f"cp ./{file_name}.rvp model/{file_name}.rvp",
    #                f"cp ./{file_name}.rvc model/{file_name}.rvc",
    #                f"OSTRICH_MPI=./OstrichMPI{newline}{newline}",
    #                f"mpirun $OSTRICH_MPI{newline}"
    #            ]
    #    }

    model_info = [
        f"# Model Type: {model_type}",
        f"# Catchment Name: {catchment_name}",
        f"# Catchment ID: {catchment_ch_id}",
        f"# Author: {author}",
        f"# Generation Date: {generation_date}"
    ]

    extra_dirs = [
        f"BeginExtraDirs",
        f"model",
        f"EndExtraDirs"
    ]

    constraints = [
        f"BeginConstraints",
        f"# not needed when no constraints, but PenaltyFunction statement above is required",
        f"# name     type     penalty    lwr   upr   resp.var",
        f"EndConstraints",
    ]
    # ost_save_output:
    #     "Save Ostrich Run Output":
    #         [
    #             "# !/bin/bash\n",
    #             f"set -e\n",
    #             f"echo \"saving Ostrich run...\n",
    #             f"if \[ ! -e "
    #             "echo - e "Processor: $1" >>.. / preserve_out.txt",
    #             "echo - e "Trial: $2" >>../ preserve_out.txt",
    # "echo - e "Counter: $3" >>../ preserve_out.txt",
    # "echo - e "Obj. Function Category: $4" >>../ preserve_out.txt",
    # "exit 0"
    #         ]
    hbv = {
        "ost_in":
            {
                "Model Info":
                    model_info,
                "General Options":
                    general_options,
                "Extra Directories":
                    extra_dirs,
                "File Pairs":
                    [
                        f"BeginFilePairs",
                        f"{file_name}.rvp.tpl;	{file_name}.rvp",
                        f"{file_name}.rvc.tpl;  {file_name}.rvc",
                        f"{file_name}.rvh.tpl;  {file_name}.rvh",
                        f"{file_name}.rvt.tpl;  {file_name}.rvt",
                        f"EndFilePairs"
                    ],
                "Parameter Specification":
                    [
                        f"#Parameter/DV Specification",
                        f"#name,initial value, lower bound, upper bound, input, output, internal transformations",
                        f"#name exactly as in *.tpl",
                        f"BeginParams",
                        f"#parameter	   init.	 low		high	tx_in  tx_ost tx_out",
                        f"{params['HBV']['names']['HBV_Param_01']}              random          {params['HBV']['lower']['HBV_Param_01']}                {params['HBV']['upper']['HBV_Param_01']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_02']}              random          {params['HBV']['lower']['HBV_Param_02']}                {params['HBV']['upper']['HBV_Param_02']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_03']}              random          {params['HBV']['lower']['HBV_Param_03']}                {params['HBV']['upper']['HBV_Param_03']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_04']}              random          {params['HBV']['lower']['HBV_Param_04']}                {params['HBV']['upper']['HBV_Param_04']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_05']}              random          {params['HBV']['lower']['HBV_Param_05']}                {params['HBV']['upper']['HBV_Param_05']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_06']}              random          {params['HBV']['lower']['HBV_Param_06']}                {params['HBV']['upper']['HBV_Param_06']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_07']}              random          {params['HBV']['lower']['HBV_Param_07']}                {params['HBV']['upper']['HBV_Param_07']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_08']}              random          {params['HBV']['lower']['HBV_Param_08']}                {params['HBV']['upper']['HBV_Param_08']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_09']}              random          {params['HBV']['lower']['HBV_Param_09']}                {params['HBV']['upper']['HBV_Param_09']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_10']}              random          {params['HBV']['lower']['HBV_Param_10']}                {params['HBV']['upper']['HBV_Param_10']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_11']}              random          {params['HBV']['lower']['HBV_Param_11']}                {params['HBV']['upper']['HBV_Param_11']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_12']}              random          {params['HBV']['lower']['HBV_Param_12']}                {params['HBV']['upper']['HBV_Param_12']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_13']}              random          {params['HBV']['lower']['HBV_Param_13']}                {params['HBV']['upper']['HBV_Param_13']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_14']}              random          {params['HBV']['lower']['HBV_Param_14']}                {params['HBV']['upper']['HBV_Param_14']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_15']}              random          {params['HBV']['lower']['HBV_Param_15']}                {params['HBV']['upper']['HBV_Param_15']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_16']}              random          {params['HBV']['lower']['HBV_Param_16']}                {params['HBV']['upper']['HBV_Param_16']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_17']}              random          {params['HBV']['lower']['HBV_Param_17']}                {params['HBV']['upper']['HBV_Param_17']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_18']}              random          {params['HBV']['lower']['HBV_Param_18']}                {params['HBV']['upper']['HBV_Param_18']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_19']}              random          {params['HBV']['lower']['HBV_Param_19']}                {params['HBV']['upper']['HBV_Param_19']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_20']}              random          {params['HBV']['lower']['HBV_Param_20']}                {params['HBV']['upper']['HBV_Param_20']}        none   none     none",
                        f"{params['HBV']['names']['HBV_Param_21']}              random          {params['HBV']['lower']['HBV_Param_21']}                {params['HBV']['upper']['HBV_Param_21']}        none   none     none",
                        f"EndParams"
                    ],
                "Tied Parameters":
                    [
                        f"BeginTiedParams",
                        f"# 1-parameter linear (TLIN = 2*XVAL) ",
                        f"{params['HBV']['names']['HBV_Param_11b']} 1 {params['HBV']['names']['HBV_Param_11']} linear 0.5 0.00 free",
                        f"{params['HBV']['names']['HBV_Param_17b']} 1 {params['HBV']['names']['HBV_Param_17']} linear 500 0.00 free",
                        f"{params['HBV']['names']['HBV_Param_15b']} 1 {params['HBV']['names']['HBV_Param_15']} linear 1.0 1.0 free",
                        f"EndTiedParams"
                    ],
                "Response Variables":
                    response_variables,
                "Tied Response Variables":
                    tied_response_variables,
                "GCOP Options":
                    gcop_options,
                "Constraints":
                    constraints,
                "Random Seed Control":
                    random_seed,
                "Algorithm Settings":
                    algorithm_settings
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
                        f"cp model/{file_name}.rvc                    model_best/{file_name}.rvc",
                        f"cp model/{file_name}.rvh                    model_best/{file_name}.rvh",
                        f"cp model/{file_name}.rvp                    model_best/{file_name}.rvp",
                        f"cp model/{file_name}.rvt                    model_best/{file_name}.rvt",
                        f"cp model/output/{file_name}_Diagnostics.csv model_best/{file_name}_Diagnostics.csv",
                        f"cp model/output/{file_name}_Hydrographs.csv model_best/{file_name}_Hydrographs.csv{newline}",
                        f"exit 0",
                    ]
            },
        "ost_raven":
            {
                "Ost-Raven":
                    ost_raven
            },
        "ost_mpi_script":
            ost_mpi_script
    }
    gr4j = {
        "ost_in":
            {
                "Model Info":
                    model_info,
                "General Options":
                    general_options,
                "Extra Directories":
                    extra_dirs,
                "File Pairs":
                    [
                        f"BeginFilePairs",
                        f"{file_name}.rvp.tpl;	{file_name}.rvp",
                        f"{file_name}.rvc.tpl;  {file_name}.rvc",
                        f"EndFilePairs"
                    ],
                "Parameter Specification":
                    [
                        f"#Parameter/DV Specification",
                        f"#name,initial value, lower bound, upper bound, input, output, internal transformations",
                        f"#name exactly as in *.tpl",
                        f"BeginParams",
                        f"#parameter	   init.	 low		high	tx_in  tx_ost tx_out",
                        f"{params['GR4J']['names']['GR4J_X1']}		random		{params['GR4J']['lower']['GR4J_X1']}		{params['GR4J']['upper']['GR4J_X1']}	none   none 	none",
                        f"{params['GR4J']['names']['GR4J_X2']}  	random	 	{params['GR4J']['lower']['GR4J_X2']}		{params['GR4J']['upper']['GR4J_X2']}	none   none 	none",
                        f"{params['GR4J']['names']['GR4J_X3']}  	random		{params['GR4J']['lower']['GR4J_X3']}		{params['GR4J']['upper']['GR4J_X3']}	none   none 	none",
                        f"{params['GR4J']['names']['GR4J_X4']}  	random		{params['GR4J']['lower']['GR4J_X4']}		{params['GR4J']['upper']['GR4J_X4']} 	none   none	none",
                        f"{params['GR4J']['names']['Cemaneige_X1']}  	random		{params['GR4J']['lower']['Cemaneige_X1']}		{params['GR4J']['upper']['Cemaneige_X1']}	none   none	none",
                        f"{params['GR4J']['names']['GR4J_Cemaneige_X2']}  	random		{params['GR4J']['lower']['GR4J_Cemaneige_X2']}		{params['GR4J']['upper']['GR4J_Cemaneige_X2']}	none   none 	none",
                        f"EndParams"
                    ],
                "Tied Parameters":
                    [
                        f"BeginTiedParams",
                        f"# 1-parameter linear (TLIN = 2*XVAL) ",
                        f"{params['GR4J']['names']['GR4J_Soil_0']}   1 {params['GR4J']['names']['GR4J_X1']} linear 500 0.00 free # SOIL[0] ",
                        f"{params['GR4J']['names']['Airsnow_Coeff']}   1 {params['GR4J']['names']['GR4J_Cemaneige_X2']} linear -1.00 0.00 free #Airsnow_Coeff"
                        f"EndTiedParams"
                    ],
                "Response Variables":
                    response_variables,
                "Tied Response Variables":
                    tied_response_variables,
                "GCOP Options":
                    gcop_options,
                "Constraints":
                    constraints,
                "Random Seed Control":
                    random_seed,
                "Algorithm Settings":
                    algorithm_settings
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
                        f"cp model/{file_name}.rvc                    model_best/{file_name}.rvc",
                        f"cp model/output/{file_name}_Diagnostics.csv model_best/{file_name}_Diagnostics.csv",
                        f"cp model/output/{file_name}_Hydrographs.csv model_best/{file_name}_Hydrographs.csv{newline}",
                        f"exit 0",
                    ]
            },
        "ost_raven":
            {
                "Ost-Raven":
                    ost_raven
            },
        "ost_mpi_script":
            ost_mpi_script,
        # "ost_save_output":
    }
    hmets = {
        "ost_in":
            {
                "Model Info":
                    model_info,
                "General Options":
                    general_options,
                "Extra Directories":
                    extra_dirs,
                "File Pairs":
                    [
                        f"BeginFilePairs",
                        f"{file_name}.rvp.tpl;	{file_name}.rvp",
                        f"{file_name}.rvc.tpl;  {file_name}.rvc",
                        f"EndFilePairs"
                    ],
                "Parameter Specification":
                    [
                        f"#Parameter/DV Specification",
                        f"#name,initial value, lower bound, upper bound, input, output, internal transformations",
                        f"#name exactly as in *.tpl",
                        f"BeginParams",
                        f"#parameter	   init.	 low		high	tx_in  tx_ost tx_out",
                        f"{params['HMETS']['names']['HMETS_Param_01']}          random          {params['HMETS']['lower']['HMETS_Param_01']}            {params['HMETS']['upper']['HMETS_Param_01']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_02']}          random          {params['HMETS']['lower']['HMETS_Param_02']}            {params['HMETS']['upper']['HMETS_Param_02']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_03']}          random          {params['HMETS']['lower']['HMETS_Param_03']}            {params['HMETS']['upper']['HMETS_Param_03']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_04']}          random          {params['HMETS']['lower']['HMETS_Param_04']}            {params['HMETS']['upper']['HMETS_Param_04']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_05a']}          random          {params['HMETS']['lower']['HMETS_Param_05a']}            {params['HMETS']['upper']['HMETS_Param_05a']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_06']}          random          {params['HMETS']['lower']['HMETS_Param_06']}            {params['HMETS']['upper']['HMETS_Param_06']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_07']}          random          {params['HMETS']['lower']['HMETS_Param_07']}            {params['HMETS']['upper']['HMETS_Param_07']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_08']}          random          {params['HMETS']['lower']['HMETS_Param_08']}            {params['HMETS']['upper']['HMETS_Param_08']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_09a']}          random          {params['HMETS']['lower']['HMETS_Param_09a']}            {params['HMETS']['upper']['HMETS_Param_09a']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_10']}          random          {params['HMETS']['lower']['HMETS_Param_10']}            {params['HMETS']['upper']['HMETS_Param_10']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_11']}          random          {params['HMETS']['lower']['HMETS_Param_11']}            {params['HMETS']['upper']['HMETS_Param_11']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_12']}          random          {params['HMETS']['lower']['HMETS_Param_12']}            {params['HMETS']['upper']['HMETS_Param_12']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_13']}          random          {params['HMETS']['lower']['HMETS_Param_13']}            {params['HMETS']['upper']['HMETS_Param_13']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_14']}          random          {params['HMETS']['lower']['HMETS_Param_14']}            {params['HMETS']['upper']['HMETS_Param_14']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_15']}          random          {params['HMETS']['lower']['HMETS_Param_15']}            {params['HMETS']['upper']['HMETS_Param_15']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_16']}          random          {params['HMETS']['lower']['HMETS_Param_16']}            {params['HMETS']['upper']['HMETS_Param_16']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_17']}          random          {params['HMETS']['lower']['HMETS_Param_17']}            {params['HMETS']['upper']['HMETS_Param_17']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_18']}          random          {params['HMETS']['lower']['HMETS_Param_18']}            {params['HMETS']['upper']['HMETS_Param_18']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_19']}          random          {params['HMETS']['lower']['HMETS_Param_19']}            {params['HMETS']['upper']['HMETS_Param_19']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_20b']}          random          {params['HMETS']['lower']['HMETS_Param_20b']}            {params['HMETS']['upper']['HMETS_Param_20b']}    none   none     none",
                        f"{params['HMETS']['names']['HMETS_Param_21b']}          random          {params['HMETS']['lower']['HMETS_Param_21b']}            {params['HMETS']['upper']['HMETS_Param_21b']}    none   none     none",
                        f"EndParams"
                    ],
                "Tied Parameters":
                    [
                        f"BeginTiedParams",
                        f"# 1-parameter linear (TLIN = 2*XVAL) ",
                        f"{params['HMETS']['names']['HMETS_Param_05b']}  2 {params['HMETS']['names']['HMETS_Param_05a']} {params['HMETS']['names']['HMETS_Param_06']} linear 0.00 1.00 1.00 0.00 free",
                        f"{params['HMETS']['names']['HMETS_Param_09b']}  2 {params['HMETS']['names']['HMETS_Param_09a']} {params['HMETS']['names']['HMETS_Param_10']} linear 0.00 1.00 1.00 0.00 free",
                        f"{params['HMETS']['names']['HMETS_Param_20a']}  1 {params['HMETS']['names']['HMETS_Param_20b']} linear 500 0.00 free",
                        f"{params['HMETS']['names']['HMETS_Param_21a']}  1 {params['HMETS']['names']['HMETS_Param_21b']} linear 500 0.00 free",
                        f"EndTiedParams"
                    ],
                "Response Variables":
                    response_variables,
                "Tied Response Variables":
                    tied_response_variables,
                "GCOP Options":
                    gcop_options,
                "Constraints":
                    constraints,
                "Random Seed Control":
                    random_seed,
                "Algorithm Settings":
                    algorithm_settings
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
                    ost_raven
            },
        "ost_mpi_script":
            ost_mpi_script
    }
    hymod = {
        "ost_in":
            {
                "Model Info":
                    model_info,
                "General Options":
                    general_options,
                "Extra Directories":
                    extra_dirs,
                "File Pairs":
                    [
                        f"BeginFilePairs",
                        f"{file_name}.rvh.tpl;	{file_name}.rvh",
                        f"{file_name}.rvp.tpl;	{file_name}.rvp",
                        f"{file_name}.rvi.tpl;	{file_name}.rvi",
                        f"EndFilePairs"
                    ],
                "Parameter Specification":
                    [
                        f"#Parameter/DV Specification",
                        f"#name,initial value, lower bound, upper bound, input, output, internal transformations",
                        f"#name exactly as in *.tpl",
                        f"BeginParams",
                        f"#parameter	   init.	 low		high	tx_in  tx_ost tx_out",
                        f"{params['HYMOD']['names']['HYMOD_Param_01']}          random          {params['HYMOD']['lower']['HYMOD_Param_01']}            {params['HYMOD']['upper']['HYMOD_Param_01']}    none   none     none",
                        f"{params['HYMOD']['names']['HYMOD_Param_02']}          random          {params['HYMOD']['lower']['HYMOD_Param_02']}            {params['HYMOD']['upper']['HYMOD_Param_02']}    none   none     none",
                        f"{params['HYMOD']['names']['HYMOD_Param_03']}          random          {params['HYMOD']['lower']['HYMOD_Param_03']}            {params['HYMOD']['upper']['HYMOD_Param_03']}    none   none     none",
                        f"{params['HYMOD']['names']['HYMOD_Param_04']}          random          {params['HYMOD']['lower']['HYMOD_Param_04']}            {params['HYMOD']['upper']['HYMOD_Param_04']}    none   none     none",
                        f"{params['HYMOD']['names']['HYMOD_Param_05']}          random          {params['HYMOD']['lower']['HYMOD_Param_05']}            {params['HYMOD']['upper']['HYMOD_Param_05']}    none   none     none",
                        f"{params['HYMOD']['names']['HYMOD_Param_06']}          random          {params['HYMOD']['lower']['HYMOD_Param_06']}            {params['HYMOD']['upper']['HYMOD_Param_06']}    none   none     none",
                        f"{params['HYMOD']['names']['HYMOD_Param_07']}          random          {params['HYMOD']['lower']['HYMOD_Param_07']}            {params['HYMOD']['upper']['HYMOD_Param_07']}    none   none     none",
                        f"{params['HYMOD']['names']['HYMOD_Param_08']}          random          {params['HYMOD']['lower']['HYMOD_Param_08']}            {params['HYMOD']['upper']['HYMOD_Param_08']}    none   none     none",
                        f"{params['HYMOD']['names']['HYMOD_Param_09']}          random          {params['HYMOD']['lower']['HYMOD_Param_09']}            {params['HYMOD']['upper']['HYMOD_Param_09']}    none   none     none",
                        f"EndParams"
                    ],
                "Response Variables":
                    response_variables,
                "Tied Response Variables":
                    tied_response_variables,
                "GCOP Options":
                    gcop_options,
                "Constraints":
                    constraints,
                "Random Seed Control":
                    random_seed,
                "Algorithm Settings":
                    algorithm_settings
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
                    ost_raven
            },
        "ost_mpi_script":
            ost_mpi_script
    }
    mohyse = {
        "ost_in":
            {
                "Model Info":
                    model_info,
                "General Options":
                    general_options,
                "Extra Directories":
                    extra_dirs,
                "File Pairs":
                    [
                        f"BeginFilePairs",
                        f"{file_name}.rvp.tpl;	{file_name}.rvp",
                        f"{file_name}.rvh.tpl;  {file_name}.rvh",
                        f"EndFilePairs"
                    ],
                "Parameter Specification":
                    [
                        f"#Parameter/DV Specification",
                        f"#name,initial value, lower bound, upper bound, input, output, internal transformations",
                        f"#name exactly as in *.tpl",
                        f"BeginParams",
                        f"#parameter	   init.	 low		high	tx_in  tx_ost tx_out",
                        f"{params['MOHYSE']['names']['MOHYSE_Param_01']}                random          {params['MOHYSE']['lower']['MOHYSE_Param_01']}          {params['MOHYSE']['upper']['MOHYSE_Param_01']}  none   none     none",
                        f"{params['MOHYSE']['names']['MOHYSE_Param_02']}                random          {params['MOHYSE']['lower']['MOHYSE_Param_02']}          {params['MOHYSE']['upper']['MOHYSE_Param_02']}  none   none     none",
                        f"{params['MOHYSE']['names']['MOHYSE_Param_03']}                random          {params['MOHYSE']['lower']['MOHYSE_Param_03']}          {params['MOHYSE']['upper']['MOHYSE_Param_03']}  none   none     none",
                        f"{params['MOHYSE']['names']['MOHYSE_Param_04']}                random          {params['MOHYSE']['lower']['MOHYSE_Param_04']}          {params['MOHYSE']['upper']['MOHYSE_Param_04']}  none   none     none",
                        f"{params['MOHYSE']['names']['MOHYSE_Param_05']}                random          {params['MOHYSE']['lower']['MOHYSE_Param_05']}          {params['MOHYSE']['upper']['MOHYSE_Param_05']}  none   none     none",
                        f"{params['MOHYSE']['names']['MOHYSE_Param_06']}                random          {params['MOHYSE']['lower']['MOHYSE_Param_06']}          {params['MOHYSE']['upper']['MOHYSE_Param_06']}  none   none     none",
                        f"{params['MOHYSE']['names']['MOHYSE_Param_07']}                random          {params['MOHYSE']['lower']['MOHYSE_Param_07']}          {params['MOHYSE']['upper']['MOHYSE_Param_07']}  none   none     none",
                        f"{params['MOHYSE']['names']['MOHYSE_Param_08']}                random          {params['MOHYSE']['lower']['MOHYSE_Param_08']}          {params['MOHYSE']['upper']['MOHYSE_Param_08']}  none   none     none",
                        f"{params['MOHYSE']['names']['MOHYSE_Param_09']}                random          {params['MOHYSE']['lower']['MOHYSE_Param_09']}          {params['MOHYSE']['upper']['MOHYSE_Param_09']}  none   none     none",
                        f"{params['MOHYSE']['names']['MOHYSE_Param_10']}                random          {params['MOHYSE']['lower']['MOHYSE_Param_10']}          {params['MOHYSE']['upper']['MOHYSE_Param_10']}  none   none     none",
                        f"EndParams"
                    ],
                "Response Variables":
                    response_variables,
                "Tied Response Variables":
                    tied_response_variables,
                "GCOP Options":
                    gcop_options,
                "Constraints":
                    constraints,
                "Random Seed Control":
                    random_seed,
                "Algorithm Settings":
                    algorithm_settings
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
                    ost_raven
            },
        "ost_mpi_script":
            ost_mpi_script
    }
    ost_params = {
        "GR4J": gr4j,
        "HBV": hbv,
        "HMETS": hmets,
        "HYMOD": hymod,
        "MOHYSE": mohyse
    }

    return ost_params[model_type]


def subsection_header(title: str) -> list[str]:
    """Generates subsection header for a .rvX file from given title.
    Args:
        title: str
            Title to be used in subsection header.
    Returns:
        subsection_head: list[str]
            List of header line as strings.
    """

    subsection_head: list[str] = [
        subsection_header_line,
        f"# ----{title}--------------------------------------------",
        subsection_header_line
    ]
    logger.debug("Subsection header written, returning and leaving function subsection_header() afterwards.")
    return subsection_head


def write_rvx(catchment_ch_id: str,
              hru_info: dict,
              model_dir: str = model_dir,
              model_type: str = model_type,
              data_dir=conf['ModelName'],
              project_dir: Path = project_dir,
              model_sub_dir: str = model_sub_dir,
              params: dict = default_params,
              template_type: str = "Raven",
              attribute_csv_dir: str = "Catchment",
              rvx_type: str = "rvi",
              author=conf['Author'],
              start_year: int = start_year,
              end_year: int = end_year,
              glacier_module: bool = False):
    """Writes .rvX file(s), either as an Ostrich or Raven template.
    Args:
        catchment_ch_id: str
            Catchment id
        hru_info: dict
            HRU information
        model_dir: str
            The directory where the model files are stored. Default is "model_dir".
        model_type: str
            The type of model to use. Default is "model_type".
        data_dir: str
            The directory where the input data is stored. Default is "data_dir".
        project_dir: Path
            The directory of the project. Default is "project_dir".
        model_sub_dir: str
            The sub-directory where the model files are stored. Default is "model_sub_dir".
        params: dict
            A dictionary of model parameters. Default is "default_params".
        template_type: str
            The type of template to use. Default is "Raven".
        attribute_csv_name: str
            The name of the attribute CSV file. Default is "CH-0057_attributes.csv".
        attribute_csv_dir: str
            The directory where the attribute CSV file is stored. Default is "Hydromap Attributes".
        rvx_type: str
            The type of RVX file to create. Default is "rvi".
        author: str
            Author name
        start_year: str
            Start year of simulation.
        end_year: str
            End year of simulation.
        glacier_module: bool
            Set True if catchment contains glacier.

    """

    import shutil
    assert model_type in config.variables.supported_models, f"Got model type: {model_type}, which is not supported, check variable \"" \
                                                            f"supported_models."
    attribute_csv_name = f"{raven_tools.config.variables.catchments[catchment_ch_id]['catchment_id']}_attributes.csv"
    csv_file = pd.DataFrame
    try:
        logger.debug(
            f"Trying to read catchment attribute CSV file {Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name)}...")
        csv_file = pandas.read_csv(Path(project_dir, data_dir, attribute_csv_dir, attribute_csv_name), sep=",",
                                   skiprows=[8],
                                   index_col='attribute_names', usecols=[0, 1])
    except:
        logger.exception("Error reading csv file into pandas...")
    file_name: str = f"{catchment_ch_id}_{model_type}.{rvx_type}"
    logger.debug(f".{rvx_type} filename set to {file_name}.")
    file_path: Path = Path(project_dir, model_dir, catchment_ch_id, model_type, file_name)
    logger.debug(f".{rvx_type} file path set to {file_path}.")
    template_sections = {}
    if template_type == "Raven":
        logger.debug(f"template_type is {template_type}.")
        logger.debug(f"Trying to generate .{rvx_type} template sections with function generate_template()...")
        template_sections = generate_template_rvx(model_type=model_type, hru_info=hru_info, csv_file=csv_file,
                                                  params=params,
                                                  param_or_name="init", start_year=start_year, end_year=end_year,
                                                  catchment_ch_id=catchment_ch_id, glacier_module=glacier_module,
                                                  data_dir=data_dir)
        logger.debug(f"Wrote .{rvx_type} template sections generated by generate_template() to dict template_sections")
    if template_type == "Ostrich":
        file_path: Path = Path(project_dir, model_dir, catchment_ch_id, model_type, file_name)
        logger.debug(f"template_type is {template_type}.")
        file_path = Path((str(file_path) + ".tpl"))
        logger.debug(f"New file_path: {file_path}")
        logger.debug(f"Trying to generate .{rvx_type}.tpl template sections with function generate_template()...")
        template_sections = generate_template_rvx(model_type=model_type, hru_info=hru_info, csv_file=csv_file,
                                                  params=params,
                                                  param_or_name="names", start_year=start_year, end_year=end_year,
                                                  catchment_ch_id=catchment_ch_id, glacier_module=glacier_module,
                                                  data_dir=data_dir)
        logger.debug(
            f"Wrote .{rvx_type}.tpl template sections generated by generate_template() to dict template_sections")

    logger.debug(f"Trying to write to file {file_path}")
    with open(file_path, 'w') as ff:
        ff.writelines(f"{line}{newline}" for line in
                      create_header(author=author, catchment_ch_id=catchment_ch_id, model=model_type,
                                    rvx_type=rvx_type))
        logger.debug("Header lines written.")
        ff.write(newline)
        logger.debug(f"template_sections dictionary:\n {template_sections}")
        for section in template_sections[rvx_type]:
            logger.debug(f"Current section: {section}")
            ff.writelines(f"{line}\n" for line in subsection_header(section))
            logger.debug("Subsection header written.")
            ff.writelines(f"{lin}\n" for lin in template_sections[rvx_type][section])
            logger.debug("Template section written")
            ff.write(newline)
    if template_type == "Raven":
        dst_path: Path = Path(project_dir, model_dir, catchment_ch_id, model_type, model_sub_dir, file_name)
        shutil.copy(file_path, dst_path)
        logger.debug("template_sections for-loop finished.")


def write_ostrich(
        model_dir: str = model_dir,
        model_type: str = model_type,
        project_dir: Path = project_dir,
        catchment_name: str = catchment_name,
        catchment_ch_id: str = "CH-0010",
        params: dict = default_params,
        ost_in: bool = True,
        save_best: bool = True,
        ost_raven: bool = True,
        ost_mpi_script: bool = True,
        max_iterations: int = 500):
    """Writes Ostrich input files ostIn.txt, save_best.sh and Ost-RAVEN.sh

    Args:
        model_dir: str
            The directory where the model files are stored. Default is "model_dir".
        model_type: str
            The type of model to use. Default is "model_type".
        project_dir: Path
            The directory of the project. Default is "project_dir".
        catchment_name: str
            The name of the catchment for which the RVX file should be created. Default is "catchment".
        catchment_ch_id: str
            Catchment id
        params: dict
            A dictionary of model parameters. Default is "default_params".
        ost_in: bool
            Flag to indicate if Ostrich input file should be created. Default is True.
        save_best: bool
            Flag to indicate if save_best.sh file should be created. Default is True.
        ost_raven: bool
            Flag to indicate if Ost-RAVEN.sh file should be created. Default is True.
        ost_mpi_script: bool
            Flag to indicate if Ostrich_MPI.sh file should be created.
        max_iterations: int
            Maximum number of Ostrich iterations.
    """
    assert model_type in config.variables.supported_models, f"Got model type: {model_type}, which is not supported, check variable \"" \
                                                            f"supported_models."
    ost_in_file_name: str = f"ostIn.txt"
    save_best_file_name: str = f"save_best.sh"
    ost_raven_file_name: str = f"Ost-RAVEN.sh"
    ost_shell_file_name: str = f"Ostrich_MPI.sh"
    logger.debug(f"filename w/o suffix set to {ost_in_file_name}.")
    template_sections = generate_template_ostrich(model_type=model_type, params=params, catchment_ch_id=catchment_ch_id,
                                                  catchment_name=catchment_name, max_iterations=max_iterations)
    logger.debug(f"Dictionary template_sections created.")
    logger.debug(f"Variable ost_in evaluated to {ost_in}")
    if ost_in:
        logger.debug(
            f"Trying to write to file: {Path(project_dir, model_dir, catchment_ch_id, model_type, ost_in_file_name)}")
        try:
            with open(Path(project_dir, model_dir, catchment_ch_id, model_type, ost_in_file_name), 'w+') as ff:
                logger.info(f"template_sections dictionary: {newline} {template_sections}")
                for section in template_sections["ost_in"]:
                    logger.info(f"Current section: {section}")
                    ff.writelines(f"{line}\n" for line in subsection_header(section))
                    logger.info("Subsection header written.")
                    ff.writelines(f"{lin}\n" for lin in template_sections["ost_in"][section])
                    logger.info("Template section written.")
                    ff.write(newline)
        except:
            logger.exception("There has been an error writing ostIn.txt")
    logger.debug(f"Variable save_best evaluated to {save_best}")
    if save_best:
        file_path: Path = Path(project_dir, model_dir, catchment_ch_id, model_type, save_best_file_name)
        logger.debug(
            f"Trying to write to file: {file_path}")
        with open(file_path, 'w') as ff:
            logger.info(f"template_sections dictionary: {newline} {template_sections}")
            for section in template_sections["save_best"]:
                ff.writelines(f"{lin}\n" for lin in template_sections["save_best"][section])
                logger.debug("Template section written.")
                ff.write(newline)
        logger.debug(f"Variable ost_raven evaluated to {ost_raven}")
        os.chmod(file_path, 0o775)
    if ost_raven:
        file_path: Path = Path(project_dir, model_dir, catchment_ch_id, model_type, ost_raven_file_name)
        logger.debug(
            f"Trying to write to file: {file_path}")
        with open(file_path, 'w') as ff:
            logger.info(f"template_sections dictionary: {newline} {template_sections}")
            for section in template_sections["ost_raven"]:
                ff.writelines(f"{lin}\n" for lin in template_sections["ost_raven"][section])
                logger.debug("Template section written.")
                ff.write(newline)
        os.chmod(file_path, 0o775)

    if ost_mpi_script:
        file_path: Path = Path(project_dir, model_dir, catchment_ch_id, model_type, ost_shell_file_name)
        with open(file_path, "w") as ff:
            for section in template_sections["ost_mpi_script"]:
                ff.writelines(f"{lin}\n" for lin in template_sections["ost_mpi_script"][section])
                logger.debug("Template section written.")
                ff.write(newline)
        os.chmod(file_path, 0o775)


if __name__ == '__main__':
    # model_dir = Path(project_config['ModelDir'])
    # model_type = Path("raven_broye_gr4j.rvt")
    start_year = 1981
    end_year = 2000
    # write_rvt(1981, 1982)
