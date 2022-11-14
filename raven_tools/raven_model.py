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

raven_filetypes = [
    "rvi",
    "rvh",
    "rvp",
    "rvc",
    "rvt"
]


class RavenModel:
    """Class that allows various operations on a selected model.

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
        logger_raven_model.debug(f"CWD: {os.getcwd()}")
        logger_raven_model.debug("Trying to open config.yaml...")
        with open("raven_tools/config/new_model_config.yaml", "r") as f:
            self.config = yaml.load(f, Loader=yaml.FullLoader)
        logger_raven_model.debug("config.yaml loaded.")
        logger_raven_model.debug("Trying to set self.X variables...")
        logger_raven_model.debug("Setting self.model_type...")
        self.model_type = model_type
        logger_raven_model.debug("Setting self.root_dir...")
        self.root_dir = Path(os.path.join(os.getcwd(), Path("RAVEN")))
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
        with open("raven_tools/config/default_params.yaml", "r") as f:
            self.default_params = yaml.load(f, Loader=yaml.FullLoader)

    def __getitem__(self, item):
        print(type(item), item)

    @property
    def model_type(self) -> str:
        """Returns model type.

        Returns model type. Supported models are: GR4J, HYMOD, HMETS, MOHYSE and HBV-EC.

        :return: Model type
        :rtype self._model_type: str
        """
        logger_raven_model.debug("Getting model type...")
        assert isinstance(self._model_type, str), f"model type should be str, is type{self._model_type} instead."
        return self._model_type

    @model_type.setter
    def model_type(self, value: str):
        """Set the model type

        Set the model type. Supported models are: GR4J, HYMOD, HMETS, MOHYSE and HBV-EC

        :param str value: Model type
        """
        logger_raven_model.debug("Setting model type...")
        assert isinstance(value, str), f"model_type expected string, got type{value} instead."
        self._model_type = value

    @property
    def root_dir(self) -> Path:
        """Get root directory

        :return self._root_dir: Root directory that contains the 'RAVEN' folder.
        :rtype self._root_dir: Path
        """
        logger_raven_model.debug("Getting root dir...")
        assert isinstance(self._root_dir, Path), f"root_dir should be Path, got type{self._root_dir} instead."
        return self._root_dir

    @root_dir.setter
    def root_dir(self, value: Path):
        """Set root directory aka project directory

        :param Path value: Root directory that contains the 'RAVEN' folder
        """
        logger_raven_model.debug("Setting root dir...")
        assert isinstance(value, Path), f"root_dir expected Path, got {type(self._root_dir)} instead."
        self._root_dir = value

    @property
    def model_dir(self) -> Path:
        """Get model directory

        :return self._model_dir:
        :rtype: Path
        """
        logger_raven_model.debug("Getting model dir...")
        assert isinstance(self._model_dir, Path), f"model_dir should be Path, is {type(self._model_dir)} instead."
        return self._model_dir

    @model_dir.setter
    def model_dir(self, value: Path):
        """Set model directory

        :param Path value: Path to model directory
        """
        logger_raven_model.debug("Setting model dir...")
        assert isinstance(value, Path), f"model_dir expected a Path, got {type(value)} instead."
        self._model_dir = value

    @property
    def catchment(self) -> str:
        """Return catchment name

        :return self._catchment: Catchment name
        :rtype: str
        """
        logger_raven_model.debug("Getting catchment name...")
        return self._catchment

    @catchment.setter
    def catchment(self, value: str):
        """Set catchment name.

        :param value: Catchment name
        :type value: str

        """
        logger_raven_model.debug("Setting catchment name...")
        self._catchment = value

    @property
    def config(self):
        """Returns configuration

        :return: Configuration
        """
        logger_raven_model.debug("Getting configuration")
        return self._config

    @config.setter
    def config(self, value):
        """Set configuration

        Set the configuration, if changed after model creation. Expects a python object created from a YAML config\
        file with structure as the template file.

        :param value: Configuration loaded from yaml.load()
        """
        logger_raven_model.debug("Setting configuration...")
        self._config = value

    @property
    def attribute_csv(self) -> str:
        """Returns name of catchment attribute CSV file.

        :return: Attribute CSV file name
        :rtype: str
        """
        logger_raven_model.debug("Getting name of catchment attribute CSV file...")
        return self._attribute_csvcsv

    @attribute_csv.setter
    def attribute_csv(self, value: str):
        """Set name of catchment attribute CSV file to be used.

        :param str value:
        """
        logger_raven_model.debug("Setting name of catchment attribute CSV file...")
        self._attribute_csv = value

    @property
    def default_params(self):
        """

        :return:
        """
        return self._default_params

    @default_params.setter
    def default_params(self, value: dict):
        self._default_params: dict = value

    def create_dirs(self):
        """Create model (sub-)directories.

        Creates model directories in the project's root dir. Follows the convention of rootdir/catchment/modelname/\
        subdirs

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

    def write_rvx(self, ostrich_template: bool = False, raven_template: bool = True, rvx_type: str = "rvi"):
        """Write .rvX file for Raven and/or Ostrich

        :param rvx_type: Suffix of the Raven file to be written (rvi,rvp,rvc,rvt or rvh)
        :param ostrich_template: Set True if Ostrich template should be generated.
        :param raven_template: Set True if Raven template should be generated.
        :type raven_template: bool
        :type ostrich_template: bool

        """
        assert rvx_type in raven_filetypes, f"Raven suffix .{rvx_type} is not in list of accepted suffixes."
        # TODO: Implement discharge and forcings time series.
        logger_raven_model.debug("Starting if-tree for template type...")
        if raven_template is True:
            logger_raven_model.debug("Variable raven_template is True...")
            logger_raven_model.debug(f"Trying to call rr.write_rvx function to create .{rvx_type} for Raven...")
            rr.write_rvx(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                         catchment=self.catchment, params=self.default_params, template_type="Raven", rvx_type=rvx_type)
            logger_raven_model.debug(f".{rvx_type} for Raven created by rr.write_rvx function")
        if ostrich_template is True:
            logger_raven_model.debug("Variable ostrich_template is True...")
            logger_raven_model.debug(f"Trying to call rr.write_rvx function to create .{rvx_type}.tpl for Ostrich...")
            rr.write_rvx(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                         catchment=self.catchment, params=self.default_params, template_type="Ostrich",
                         rvx_type=rvx_type)
            logger_raven_model.debug(f".{rvx_type}.tpl for Ostrich created by rr.write_rvx function")
        if ostrich_template is False and raven_template is False:
            logger_raven_model.debug("Variables ostrich_template and raven_template set to False.")
            print("You have not selected a template type...")
            logger_raven_model.debug("No template file needed to be written by function rr.write_rvx.")
            pass

    def write_ost_in(self):
        """Write .rvX file for Raven and/or Ostrich

        :param rvx_type: Suffix of the Raven file to be written (rvi,rvp,rvc,rvt or rvh)
        :param ostrich_template: Set True if Ostrich template should be generated.
        :param raven_template: Set True if Raven template should be generated.
        :type raven_template: bool
        :type ostrich_template: bool

        """

        logger_raven_model.debug("Variable raven_template is True...")
        logger_raven_model.debug(f"Trying to call rr.write_ostrich_files() function to create ostrich input files\
            for Raven...")
        rr.write_ost_in(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                        catchment=self.catchment, params=self.default_params)
        logger_raven_model.debug(f"ost_in.txt for Raven created by rr.write_rvx function")

    def write_save_best(self):
        logger_raven_model.debug("Variable raven_template is True...")
        logger_raven_model.debug(f"Trying to call rr.write_ostrich_files() function to create ostrich input files\
            for Raven...")
        rr.write_save_best(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                        catchment=self.catchment)
        logger_raven_model.debug(f"ost_in.txt for Raven created by rr.write_rvx function")
