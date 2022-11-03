import logging
import os
import sys
from pathlib import Path

import yaml

import raven_run as rr

handler = logging.StreamHandler(sys.stdout)
frm = logging.Formatter("{levelname}: {message} ({filename}/{funcName}/{lineno})",
                        style="{")
handler.setFormatter(frm)
logger_raven_model = logging.getLogger("raven_model")
logger_raven_model.addHandler(handler)
logger_raven_model.setLevel(logging.DEBUG)
logger_raven_model.debug("Logging to console started")

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


class RavenModel:
    """

    """

    def __init__(self, model_type="model_type", catchment="catchment"):
        """

        :param model_type:
        :param catchment:
        """
        logger_raven_model.debug(f"Starting __init__ of {__name__}...")
        assert isinstance(model_type, str), f"model_type expected a string, got {type(model_type)} instead"
        assert isinstance(catchment, str), f"catchment expected a string, got {type(model_type)} instead"
        assert model_type in supported_models, f"model_type expected GR4J, HYMOD, HMETS, HBV or MOHYSE, got {model_type} instead "
        logger_raven_model.debug("Trying to open config.yaml...")
        with open("raven_tools/new_model_config.yaml", "r") as f:
            self.config = yaml.load(f, Loader=yaml.FullLoader)
        logger_raven_model.debug("config.yaml loaded.")
        logger_raven_model.debug("Trying to set self.X variables...")
        logger_raven_model.debug("Setting self.model_type...")
        self.model_type = model_type
        logger_raven_model.debug("Setting self.root_dir...")
        self.root_dir = Path(os.getcwd(), "RAVEN")
        logger_raven_model.debug("Setting self.catchment...")
        self.catchment = catchment
        logger_raven_model.debug("Setting self.attribute_csv_name (file name with catchment attributes...")
        self.attribute_csv = "CH-0057_attributes.csv"
        logger_raven_model.debug("Setting self.model_dir...")
        self.model_dir = Path(self.root_dir, "models", self.catchment, self.model_type)
        logger_raven_model.debug("Setting self.dirs...")
        logger_raven_model.debug("Setting self.model_type...")
        self.dirs = [
            Path(self.model_dir, "model"),
            Path(self.model_dir, "model", "output")
        ]
        logger_raven_model.debug("Self.X variables set.")
        logger_raven_model.debug(f"__init__ of {__name__} finished...")

    @property
    def model_type(self) -> str:
        """Returns model type.

        Returns model type. Supported models are: GR4J, HYMOD, HMETS, MOHYSE and HBV-EC.

        :return: Model type
        """
        logger_raven_model.debug("Getting model type...")
        return self._model_type

    @model_type.setter
    def model_type(self, value: str):
        """Set the model type

        Set the model type. Supported models are: GR4J, HYMOD, HMETS, MOHYSE and HBV-EC

        :param str value: Model type
        """
        logger_raven_model.debug("Setting model type...")
        self._model_type = value

    @property
    def root_dir(self):
        """

        :return:
        """
        logger_raven_model.debug("Getting root dir...")
        return self._root_dir

    @root_dir.setter
    def root_dir(self, value):
        """

        :param value:
        """
        logger_raven_model.debug("Setting root dir...")
        self._root_dir = value

    @property
    def model_dir(self):
        """

        :return:
        """
        logger_raven_model.debug("Getting model dir...")
        return self._model_dir

    @model_dir.setter
    def model_dir(self, value):
        """

        :param value:
        """
        logger_raven_model.debug("Setting model dir...")
        self._model_dir = Path(self.root_dir, "models", self.catchment, self.model_type)

    @property
    def catchment(self):
        """

        :return:
        """
        logger_raven_model.debug("Getting catchment name...")
        return self._catchment

    @catchment.setter
    def catchment(self, value):
        """

        :param value:
        """
        logger_raven_model.debug("Setting catchment name...")
        self._catchment = value

    @property
    def config(self):
        """

        :return:
        """
        logger_raven_model.debug("Getting configuration")
        return self._config

    @config.setter
    def config(self, value):
        """

        :param value:
        """
        logger_raven_model.debug("Setting configuration...")
        self._config = value

    @property
    def attribute_csv(self):
        """

        :return:
        """
        logger_raven_model.debug("Getting name of catchment attribute CSV file...")
        return

    @attribute_csv.setter
    def attribute_csv(self, value):
        """

        :param value:
        """
        logger_raven_model.debug("Setting name of catchment attribute CSV file...")
        self._attribute_csv = value

    def create_dirs(self):
        """

        """
        logger_raven_model.debug("Trying to create model directories...")
        try:
            for f in self.dirs:
                logger_raven_model.debug(f"Creating directory: {f}")
                f.mkdir(parents=True)
                print(f"Directory created: {f}")
        except FileExistsError:
            logger_raven_model.debug(f"Directory {f} already exists, skipping...")
            print("Directory already exists...")

    def write_rvh(self):
        """

        """
        rr.write_rvh(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                     catchment=self.catchment)

    def write_rvp(self, params=default_params, ostrich_template=False, raven_template=True):
        """

        """
        logger_raven_model.debug("Starting if-tree for template type...")
        if raven_template is True:
            logger_raven_model.debug("Variable raven_template is True...")
            logger_raven_model.debug("Trying to call rr.write_rvp function to create .rvp for Raven...")
            rr.write_rvp(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                         catchment=self.catchment, params=params, template_type="Raven")
            logger_raven_model.debug(".rvp for Raven created by rr.write_rvp function")
        if ostrich_template is True:
            logger_raven_model.debug("Variable ostrich_template is True...")
            logger_raven_model.debug("Trying to call rr.write_rvp function to create .rvp.tpl for Ostrich...")
            rr.write_rvp(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                         catchment=self.catchment, params=params, template_type="Ostrich")
            logger_raven_model.debug(".rvp.tpl for Ostrich created by rr.write_rvp function")
        if ostrich_template is False and raven_template is False:
            logger_raven_model.debug("Variables ostrich_template and raven_template set to False.")
            print("You have not selected a template type...")
            logger_raven_model.debug("No template file needed to be written by function rr.write_rvp.")
            pass
