import os
from pathlib import Path

import yaml

os.path.abspath(Path(os.path.abspath(__file__), ".."))
config_filepath = Path(os.path.abspath(Path(os.path.abspath(__file__), "..")),
                       "new_model_config.yaml")
params_filepath = Path(os.path.abspath(Path(os.path.abspath(__file__), "..")), "default_params.yaml")

# Read project_config.yaml file
try:
    with open(config_filepath, "r") as f:
        project_config = yaml.load(f, Loader=yaml.FullLoader)
except FileNotFoundError:
    logger.exception("project_config.yaml file could not be found!")

# Read default_params.yaml file
try:
    with open(params_filepath, "r") as f:
        default_params = yaml.load(f, Loader=yaml.FullLoader)
except FileNotFoundError:
    logger.exception("default_params.yaml file could not be found!")

supported_models = [
    "GR4J",
    "HYMOD",
    "HMETS",
    "HBV",
    "MOHYSE"
]

raven_filetypes = [
    "rvi",
    "rvh",
    "rvp",
    "rvc",
    "rvt"
]

# Station coordinates are from https://www.hydrodaten.admin.ch/en/{STATION_ID}.html

catchments = {
    "Ticino":
        {"ID": 2020,
         "lat": 721245,
         "lon": 117025,
         "short_code": "TicBel",
         "station_elevation": "220",
         "catchment_id": "CH-0053"
         },
    "Broye":
        {"ID": 2034,
         "lat": 561660,
         "lon": 187320,
         "short_code": "BroPay",
         "station_elevation": "441",
         "catchment_id": "CH-0057"
         },
    "Thur":
        {"ID": 2044,
         "lat": 693510,
         "lon": 272500,
         "short_code": "ThuAnd",
         "station_elevation": "356",
         "catchment_id": "CH-0058"
         },
    "Massa":
        {"ID": 2161,
         "lat": 643700,
         "lon": 137290,
         "short_code": "MasBla",
         "station_elevation": "1446",
         "catchment_id": "CH-0105"
         },
    "Weisse_Luetschine":
        {"ID": 2200,
         "lat": 635310,
         "lon": 164550,
         "short_code": "WLuZwe",
         "station_elevation": "650",
         "catchment_id": "CH-0118"
         },
    "Dischmabach":
        {"ID": 2327,
         "lat": 786220,
         "lon": 183370,
         "short_code": "DisDav",
         "station_elevation": "1668",
         "catchment_id": "CH-0159"
         }
}

camels_column_names = {
    "date": "date",
    "discharge_vol(m3/s)": "discharge"
}

forcings_dirs = [
    "RhiresD_v2.0_swiss.lv95",
    "SrelD_v2.0_swiss.lv95",
    "TabsD_v2.0_swiss.lv95",
    "TmaxD_v2.0_swiss.lv95",
    "TminD_v2.0_swiss.lv95"
]

forcing_netcdf_names = [
    "RhiresD",
    "SrelD",
    "TabsD",
    "TmaxD",
    "TminD"
]
