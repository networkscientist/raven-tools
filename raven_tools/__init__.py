"""
API reference documentation for the 'raven_tools' package.
"""

import os
import yaml
import logging
import sys
from pathlib import Path
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
frm = logging.Formatter("{levelname}: {message} ({filename}/{funcName}/{lineno})",
                        style="{")
handler.setFormatter(frm)
logger.addHandler(handler)

os.path.abspath(Path(os.path.abspath(__file__),".."))
config_filepath = Path(os.path.abspath(Path(os.path.abspath(__file__),"..")),"config","config.yaml")
params_filepath = Path(os.path.abspath(Path(os.path.abspath(__file__),"..")),"config","default_params.yaml")

# Read config.yaml file
try:
    with open(config_filepath, "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
except FileNotFoundError:
    logger.exception("config.yaml file could not be found!")

# Read default_params.yaml file
try:
    with open(params_filepath, "r") as f:
        params = yaml.load(f, Loader=yaml.FullLoader)
except FileNotFoundError:
    logger.exception("default_params.yaml file could not be found!")