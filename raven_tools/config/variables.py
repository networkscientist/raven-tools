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
         },
    "Broye":
        {"ID": 2034,
         "lat": 561660,
         "lon": 187320
         },
    "Thur":
        {"ID": 2044,
         "lat": 693510,
         "lon": 272500
         },
    "Massa":
        {"ID": 2161,
         "lat": 643700,
         "lon": 137290
         },
    "Weisse Lütschine":
        {"ID": 2200,
         "lat": 635310,
         "lon": 164550
         },
    "Dischmabach":
        {"ID": 2327,
         "lat": 786220,
         "lon": 183370
         }
}
