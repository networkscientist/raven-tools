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


def write_rvh(model_dir=config['ModelDir'], model_name=config['ModelName'], data_dir=config['DataDir'],
              project_dir=config['ProjectDir'], catchment=config['Catchment'], model_sub_dir=config['ModelSubDir'],
              attribute_csv="CH-0057_attributes.csv"):
    assert model_name in supported_models, f"Model type {model_name} not in list of supported models."
    csv_file = pandas.read_csv(Path(project_dir, data_dir, "Hydromap Attributes", attribute_csv), sep=",",
                               skiprows=[8],
                               index_col='attribute_names')
    file_name: str = f"{catchment}_{model_name}.rvh"
    file_path: Path = Path(project_dir, model_dir, catchment, model_name, model_sub_dir, file_name)
    if model_name == "GR4J":
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

    elif model_name == "HYMOD":
        pass
    elif model_name == "HMETS":
        pass
    with open(file_path, 'w') as ff:
        ff.writelines(f"{line}{newline}"for line in header)
        ff.write(newline)
        ff.writelines(f"{line}{newline}" for line in sub_section_header("Sub-Basins"))
        ff.writelines(f"{line}{newline}" for line in sub_basins)
        ff.write(newline)
        ff.writelines(f"{line}{newline}" for line in sub_section_header("HRUs"))
        ff.writelines(f"{line}{newline}" for line in hrus)


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


def write_rvp(model_dir=config['ModelDir'], model_name=config['ModelName'], data_dir=config['DataDir'],
              project_dir=config['ProjectDir'], catchment=config['Catchment'], model_sub_dir=config['ModelSubDir']):
    """

    """

    csv_file = pandas.read_csv(Path(project_dir, data_dir, "Hydromap Attributes", "CH-0057_attributes.csv"), sep=",",
                               skiprows=[8],
                               index_col='attribute_names')
    file_name: str = f"{catchment}_{model_name}.rvp"
    file_path: Path = Path(project_dir, model_dir, catchment, model_name, model_sub_dir, file_name)
    soil_classes = [
        ":SoilClasses",
        "   :Attributes",
        "   :Units",
        "       SOIL_PROD",
        "       SOIL_ROUT",
        "       SOIL_TEMP",
        "       SOIL_GW",
        ":EndSoilClasses"
    ]
    soil_profiles = [
        "#     name,#horizons,{soiltype,thickness}x{#horizons}",
        "#     GR4J_X1 is thickness of first layer (SOIL_PROD), here 0.529",
        ":SoilProfiles",
        "   DEFAULT_P, 4, SOIL_PROD, 1.881244E+00, SOIL_ROUT, 0.300, SOIL_TEMP, 1.000, SOIL_GW, 1.000,",
        ":EndSoilProfiles"
    ]
    vegetation_classes = [
        ":VegetationClasses",
        "   :Attributes, MAX_HT, MAX_LAI, MAX_LEAF_COND",
        "   :Units, m, none, mm_per_s",
        "   VEG_ALL, 0.0, 0.0, 0.0",
        ":EndVegetationClasses"
    ]
    land_use_classes = [
        ":LandUseClasses",
        "   :Attributes, IMPERM, FOREST_COV",
        "   :Units, frac, frac",
        f"   LU_ALL, {int(csv_file.loc['a0425_clc18_5_aurbv1_0']['values']) / 100}, {int(csv_file.loc['a0418_clc18_5_afrtv1_0']['values']) / 100}",
        ":EndLandUseClasses"
    ]

    global_parameters = [
        ":GlobalParameter RAINSNOW_TEMP       0.0",
        ":GlobalParameter RAINSNOW_DELTA      1.0",
        ":GlobalParameter AIRSNOW_COEFF     6.649411E-01 # [1/d] = 1.0 - CEMANEIGE_X2 = 1.0 - x6",
        ":GlobalParameter AVG_ANNUAL_SNOW    16.9 # [mm]  =       CEMANEIGE_X1 =       x5",
        "#:GlobalParameter PRECIP_LAPSE     0.0004 I assume not necessary for gridded data",
        "#:GlobalParameter ADIABATIC_LAPSE  0.0065 not necessary for gridded data"
    ]

    soil_parameters = [
        ":SoilParameterList",
        "   :Parameters, POROSITY, GR4J_X3, GR4J_X2",
        "   :Units, none, mm, mm / d",
        "   [DEFAULT], 1.0, 2.742079E+02, -1.342446E+01",
        ":EndSoilParameterList"
    ]
    land_use_parameters = [
        ":LandUseParameterList",
        "   :Parameters, GR4J_X4, MELT_FACTOR",
        "   :Units, d, mm / d / C",
        "   [DEFAULT], 7.852040E-01, 2.830836E+01",
        ":EndLandUseParameterList"
    ]

    with open(file_path, 'w') as ff:
        newline = "\n"
        ff.writelines(f"{line}{newline}" for line in header)
        ff.write(newline)
        ff.writelines(f"{line}\n" for line in sub_section_header("Soil Classes"))
        ff.writelines(f"{line}\n" for line in soil_classes)
        ff.write(newline)
        ff.writelines(f"{line}\n" for line in sub_section_header("Soil Profiles"))
        ff.writelines(f"{line}\n" for line in soil_profiles)
        ff.write(newline)
        ff.writelines(f"{line}\n" for line in sub_section_header("Vegetation Classes"))
        ff.writelines(f"{line}\n" for line in vegetation_classes)
        ff.write(newline)
        ff.writelines(f"{line}\n" for line in sub_section_header("Land Use Classes"))
        ff.writelines(f"{line}\n" for line in land_use_classes)
        ff.write(newline)
        ff.writelines(f"{line}\n" for line in sub_section_header("Global Parameters"))
        ff.writelines(f"{line}\n" for line in global_parameters)
        ff.write(newline)
        ff.writelines(f"{line}\n" for line in sub_section_header("Soil Parameters"))
        ff.writelines(f"{line}\n" for line in soil_parameters)
        ff.write(newline)
        ff.writelines(f"{line}\n" for line in sub_section_header("Land Use Parameters"))
        ff.writelines(f"{line}\n" for line in land_use_parameters)


def sub_section_header(title: str):
    """

    :param title:
    :return:
    """
    sub_section_head = [
        sub_section_header_line,
        f"# ----{title}--------------------------------------------",
        sub_section_header_line
    ]
    return sub_section_head


if __name__ == '__main__':
    model_dir = Path(config['ModelDir'])
    model_name = Path("raven_broye_gr4j.rvt")
    start_year = 1981
    end_year = 2000
    write_rvt(1981, 1982, model_dir, model_name)
