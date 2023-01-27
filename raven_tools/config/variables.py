import logging
import os
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)
logger.debug("Logging from variables.py started.")

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
