import os
from pathlib import Path

import yaml

import raven_run as rr

logger = rt.logger
logger.debug("Logging from raven_model to console started")

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
    """Class that allows various operations on a selected model. There are methods to create the suggested directory
    structure, to create .rvX RAVEN configuration files and to write Ostrich configuration files.

    Args:

    Attributes:

    """

    def __init__(self, model_type: str = "GR4J", catchment: str = "default"):
        """
        Args:
            model_type (str): Name of model_type (GR4J, HYMOD, HMETS, HBV or MOHYSE)
            catchment (str): Name of catchment
        """
        logger.debug(f"Starting __init__ of {__name__}...")
        assert isinstance(model_type, str), f"model_type expected a string, got {type(model_type)} instead"
        assert isinstance(catchment, str), f"catchment expected a string, got {type(model_type)} instead"
        assert model_type in supported_models, f"model_type expected GR4J, HYMOD, HMETS, HBV or MOHYSE, got {model_type} instead "
        with open("raven_tools/config/new_model_config.yaml", "r") as f:
        logger.debug(f"CWD: {os.getcwd()}")
        logger.debug("Trying to open config.yaml...")
            self.config = yaml.load(f, Loader=yaml.FullLoader)
        logger.debug("config.yaml loaded.")
        logger.debug("Trying to set self.X variables...")
        logger.debug("Setting self.model_type...")
        self.model_type = model_type
        logger.debug("Setting self.root_dir...")
        self.root_dir = Path(os.path.join(os.getcwd(), Path("RAVEN")))
        logger.debug("Setting self.catchment...")
        self.catchment = catchment
        logger.debug("Setting self.attribute_csv_name (file name with catchment attributes...")
        self.attribute_csv = "CH-0057_attributes.csv"
        logger.debug("Setting self.model_dir...")
        self.model_dir = Path(self.root_dir, "models", self.catchment, self.model_type)
        logger.debug("Setting self.dirs...")
        logger.debug("Setting self.model_type...")
        self.dirs = [
            Path(self.model_dir, "model"),
            Path(self.model_dir, "model", "output"),
            Path(self.model_dir, "model", "data_obs")
        ]
        logger.debug("Setting self.data_dir...")
        self.data_dir = Path(self.root_dir, self.config['DataDir'])
        logger.debug("Self.X variables set.")
        logger.debug(f"__init__ of {__name__} finished...")
        with open("raven_tools/config/default_params.yaml", "r") as f:
            self.default_params = yaml.load(f, Loader=yaml.FullLoader)

    def __getitem__(self, item):
        print(type(item), item)

    @property
    def data_dir(self) -> Path:
        logger.debug("Getting data dir...")
        return self._data_dir

    @data_dir.setter
    def data_dir(self, value: Path):
        logger.debug("Setting data dir...")
        assert isinstance(value, Path), f"data_dir should be of type Path, is type {self._data_dir} instead."
        self._data_dir = value

    @property
    def model_type(self) -> str:
        """Returns model type. Supported models are GR4J, HYMOD, HMETS, MOHYSE and HBV-EC."""

        logger.debug("Getting model type...")
        assert isinstance(self._model_type, str), f"model type should be str, is type{self._model_type} instead."
        return self._model_type

    @model_type.setter
    def model_type(self, value: str):
        logger.debug("Setting model type...")
        assert isinstance(value, str), f"model_type expected string, got type{value} instead."
        self._model_type = value

    @property
    def root_dir(self) -> Path:
        """Get root directory that contains the 'RAVEN' folder."""

        logger.debug("Getting root dir...")
        assert isinstance(self._root_dir, Path), f"root_dir should be Path, got type{self._root_dir} instead."
        return self._root_dir

    @root_dir.setter
    def root_dir(self, value: Path):
        logger.debug("Setting root dir...")
        assert isinstance(value, Path), f"root_dir expected Path, got {type(self._root_dir)} instead."
        self._root_dir = value

    @property
    def model_dir(self) -> Path:
        """Get model directory"""

        logger.debug("Getting model dir...")
        assert isinstance(self._model_dir, Path), f"model_dir should be Path, is {type(self._model_dir)} instead."
        return self._model_dir

    @model_dir.setter
    def model_dir(self, value: Path):
        logger.debug("Setting model dir...")
        assert isinstance(value, Path), f"model_dir expected a Path, got {type(value)} instead."
        self._model_dir = value

    @property
    def catchment(self) -> str:
        """Return catchment name

        :return self._catchment: Catchment name
        :rtype: str
        """
        logger.debug("Getting catchment name...")
        return self._catchment

    @catchment.setter
    def catchment(self, value: str):
        """Set catchment name.

        :param value: Catchment name
        :type value: str

        """
        logger.debug("Setting catchment name...")
        self._catchment = value

    @property
    def config(self):
        """Returns configuration

        :return: Configuration
        """
        logger.debug("Getting configuration")
        return self._config

    @config.setter
    def config(self, value):
        """Set configuration

        Set the configuration, if changed after model creation. Expects a python object created from a YAML config\
        file with structure as the template file.

        :param value: Configuration loaded from yaml.load()
        """
        logger.debug("Setting configuration...")
        self._config = value

    @property
    def attribute_csv(self) -> str:
        """Returns name of catchment attribute CSV file.

        :return: Attribute CSV file name
        :rtype: str
        """
        logger.debug("Getting name of catchment attribute CSV file...")
        return self._attribute_csvcsv

    @attribute_csv.setter
    def attribute_csv(self, value: str):
        """Set name of catchment attribute CSV file to be used.

        :param str value:
        """
        logger.debug("Setting name of catchment attribute CSV file...")
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

        logger.debug("Trying to create model directories...")
        for f in self.dirs:
            try:
                logger.debug(f"Creating directory: {f}")
                f.mkdir(parents=True)
                print(f"Directory created: {f}")
            except FileExistsError:
                logger.debug(f"Directory {f} already exists, skipping...")
                print("Directory already exists...")
                pass

    def create_symlinks(self):
        logger.debug("Trying to create data symlinks...")
        try:
            src = ["RhiresD_v2.0_swiss.lv95",
                   "SrelD_v2.0_swiss.lv95",
                   "TabsD_v2.0_swiss.lv95",
                   "TmaxD_v2.0_swiss.lv95",
                   "TminD_v2.0_swiss.lv95"]
            for s in src:
                dst = Path(self.model_dir, "model", "data_obs", s)
                logger.debug(f"Symlink src: MeteoSwiss_gridded_products/{s}")
                logger.debug(f"Symlink dst: {dst}")
                os.symlink(Path(self.data_dir, "MeteoSwiss_gridded_products", s), dst)
                logger.debug(f"Symlink created")
                print(f"Symlink created:\n"
                      f"Source: MeteoSwiss_gridded_products/{s}\n"
                      f"Destination: {dst}")
            src = Path(self.data_dir, "Discharge", "BroPay_Q_2034_daily.rvt")
            logger.debug(f"Symlink src: {src}")
            logger.debug(f"Symlink dst: {dst}")
            os.symlink(src, dst)
            logger.debug(f"Symlink created")
            print(f"Symlink created:\n"
                  f"Source: {src}\n"
                  f"Destination: {dst}")
        except:
            logger.debug(f"Error creating symlinks...")
            logger.exception("Error creating symlinks...")
            print("There has been an error creating symlinks...")

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
        logger.debug("Starting if-tree for template type...")
        if raven_template is True:
            logger.debug("Variable raven_template is True...")
            logger.debug(f"Trying to call rr.write_rvx function to create .{rvx_type} for Raven...")
            rr.write_rvx(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                         catchment=self.catchment, params=self.default_params, template_type="Raven", rvx_type=rvx_type)
            logger.debug(f".{rvx_type} for Raven created by rr.write_rvx function")
        if ostrich_template is True:
            logger.debug("Variable ostrich_template is True...")
            logger.debug(f"Trying to call rr.write_rvx function to create .{rvx_type}.tpl for Ostrich...")
            rr.write_rvx(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                         catchment=self.catchment, params=self.default_params, template_type="Ostrich",
                         rvx_type=rvx_type)
            logger.debug(f".{rvx_type}.tpl for Ostrich created by rr.write_rvx function")
        if ostrich_template is False and raven_template is False:
            logger.debug("Variables ostrich_template and raven_template set to False.")
            print("You have not selected a template type...")
            logger.debug("No template file needed to be written by function rr.write_rvx.")
            pass

    def write_ost(self,
                  ost_in=True,
                  save_best=True,
                  ost_raven=True):
        """Write Ostrich files ostIn.txt, save_best.sh and ost-raven.sh

        :param ost_raven:
        :param ost_in:
        :param save_best:
        """

        logger.debug("Variable raven_template is True...")
        logger.debug(f"Trying to call rr.write_ostrich_files() function to create ostrich input files\
            for Raven...")
        rr.write_ostrich(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                         catchment=self.catchment, params=self.default_params, ost_in=ost_in, save_best=save_best,
                         ost_raven=ost_raven)
        logger.debug(f"ostIn.txt for Raven created by rr.write_rvx function")
