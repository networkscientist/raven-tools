"""
Work with a Raven class.
"""

import logging
import os
from pathlib import Path

# import raven_tools as rt

logger = logging.getLogger(__name__)
logger.debug("Logging from raven_model to console started")
from raven_tools import config
from raven_tools.processing import raven_run as rr
from raven_tools.processing import raven_preprocess as rpe


class RavenModel:
    """Class that allows various operations on a selected model. There are methods to create the suggested directory
    structure, to create .rvX RAVEN configuration files and to write Ostrich configuration files.

    Args:

    Attributes:

    """

    def __init__(self, model_type: str = "GR4J", catchment: str = "default", catchment_id: str = "CH-0057"):
        """
        Args:
            model_type (str): Name of model_type (GR4J, HYMOD, HMETS, HBV or MOHYSE)
            catchment (str): Name of catchment
        """
        logger.debug(f"Starting __init__ of {__name__}...")
        assert isinstance(model_type, str), f"model_type expected a string, got {type(model_type)} instead"
        assert isinstance(catchment, str), f"catchment expected a string, got {type(catchment)} instead"
        assert isinstance(catchment_id, str), f"catchment_id expected a string, got {type(catchment_id)} instead"
        self.supported_models = config.variables.supported_models
        self.raven_filetypes = config.variables.raven_filetypes
        assert model_type in self.supported_models, f"model_type expected GR4J, HYMOD, HMETS, HBV or MOHYSE, got {model_type} instead "
        logger.debug(f"CWD: {os.getcwd()}")
        logger.debug("Trying to open project_config.yaml...")
        try:
            self.conf = config.variables.project_config
        except:
            logger.debug("Error getting project_config file from __init__.py")
        logger.debug("project_config.yaml loaded.")
        logger.debug("Trying to set self.X variables...")
        logger.debug("Setting self.model_type...")
        self.model_type = model_type
        logger.debug("Setting self.root_dir...")
        self.root_dir = Path(os.path.join(os.getcwd(), Path("RAVEN")))
        logger.debug("Setting self.catchment...")
        self.catchment = catchment
        logger.debug("Setting self.catchment_id...")
        self.catchment_id = catchment_id
        logger.debug("Setting self.attribute_csv_name (file name with catchment attributes...")
        self.attribute_csv = f"{self.catchment_id}_attributes.csv"
        logger.debug("Setting self.model_dir...")
        self.model_dir = Path(self.root_dir, "models", self.catchment, self.model_type)
        logger.debug("Setting self.dirs...")
        logger.debug("Setting self.model_type...")
        self.dirs: list[Path] = [
            Path(self.model_dir, ""),
            Path(self.model_dir, "", "output"),
            Path(self.model_dir, "", "data_obs")
        ]
        logger.debug("Setting self.data_dir...")
        self.data_dir: Path = Path(self.root_dir, self.conf['DataDir'])
        logger.debug("Setting meteo folder...")
        self.meteo_dir = self.conf['MeteoSubDir']
        logger.debug("Setting start year...")
        self.start_year = self.conf['StartYear']
        logger.debug("Setting end year...")
        self.end_year = self.conf['EndYear']
        logger.debug("Self.X variables set.")
        logger.debug(f"__init__ of {__name__} finished...")
        try:
            self.default_params = config.variables.default_params
        except:
            logger("Error getting default_params file from __init__.py")

    def __getitem__(self, item):
        print(type(item), item)

    @property
    def data_dir(self) -> Path:
        """

        Returns:

        """
        logger.debug("Getting data dir...")
        return self._data_dir

    @data_dir.setter
    def data_dir(self, value: Path):
        logger.debug("Setting data dir...")
        assert isinstance(value, Path), f"data_dir should be of type Path, is type {type(value)} instead."
        self._data_dir = value

    @property
    def model_type(self) -> str:
        """Returns model type. Supported models are GR4J, HYMOD, HMETS, MOHYSE and HBV-EC."""

        logger.debug("Getting model type...")
        assert isinstance(self._model_type, str), f"model type should be str, is type{type(self._model_type)} instead."
        return self._model_type

    @model_type.setter
    def model_type(self, value: str):
        logger.debug("Setting model type...")
        assert isinstance(value, str), f"model_type expected string, got type{type(value)} instead."
        self._model_type = value

    @property
    def root_dir(self) -> Path:
        """Get root directory that contains the 'RAVEN' folder."""

        logger.debug("Getting root dir...")
        assert isinstance(self._root_dir, Path), f"root_dir should be Path, got type{type(self._root_dir)} instead."
        return self._root_dir

    @root_dir.setter
    def root_dir(self, value: Path):
        logger.debug("Setting root dir...")
        assert isinstance(value, Path), f"root_dir expected Path, got {type(value)} instead."
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
    def conf(self):
        """Returns configuration

        :return: Configuration
        """
        logger.debug("Getting configuration")
        return self._config

    @conf.setter
    def conf(self, value):
        """Set configuration

        Set the configuration, if changed after model creation. Expects a python object created from a YAML project_config\
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
        return self._attribute_csv

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

    @property
    def meteo_dir(self) -> str:
        """Return directory with meteorological forcings files.

        Returns:


        """
        return self._meteo_dir

    @meteo_dir.setter
    def meteo_dir(self, value: str):
        self._meteo_dir: str = value

    @property
    def start_year(self) -> int:
        """Returns start year."""

        logger.debug("Getting start year...")
        assert isinstance(self._start_year, int), f"Start year should be int, is type{type(self._start_year)} instead."
        return self._start_year

    @start_year.setter
    def start_year(self, value: int):
        logger.debug("Setting start year ...")
        assert isinstance(value, int), f"start_year expected int, got {type(value)} instead."
        self._start_year = value

    @property
    def end_year(self) -> int:
        """Returns end year."""

        logger.debug("Getting end year...")
        assert isinstance(self._end_year, int), f"end_year should be int, is type{type(self._end_year)} instead."
        return self._end_year

    @end_year.setter
    def end_year(self, value: int):
        logger.debug("Setting end year ...")
        assert isinstance(value, int), f"end_year expected int, got {type(value)} instead."
        self._end_year = value

    @property
    def supported_models(self) -> list[str]:
        """Returns supported_models."""
        assert isinstance(self._supported_models,
                          list), f"supported_models should be list[str], is type{type(self._supported_models)} instead."
        return self._supported_models

    @supported_models.setter
    def supported_models(self, value: list[str]):
        assert isinstance(value,
                          list), f"supported_models should be list[str], is type{type(value)} instead."
        self._supported_models = value

    @property
    def raven_filetypes(self) -> list[str]:
        """Returns raven_filetypes."""
        assert isinstance(self._raven_filetypes,
                          list), f"raven_filetypes should be list[str], is type{type(self._raven_filetypes)} instead."
        return self._raven_filetypes

    @raven_filetypes.setter
    def raven_filetypes(self, value: list[str]):
        assert isinstance(value,
                          list), f"raven_filetypes should be list[str], is type{type(value)} instead."
        self._raven_filetypes = value

    @property
    def catchment_id(self) -> str:
        """Returns catchment_id."""
        assert isinstance(self._catchment_id,
                          str), f"catchment_id should be str, is type {type(self._catchment_id)} instead."
        return self._catchment_id

    @catchment_id.setter
    def catchment_id(self, value: str):
        assert isinstance(value, str), f"catchment_id should be str, is type {type(value)} instead."
        self._catchment_id = value

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
        src = ["RhiresD_v2.0_swiss.lv95",
               "SrelD_v2.0_swiss.lv95",
               "TabsD_v2.0_swiss.lv95",
               "TmaxD_v2.0_swiss.lv95",
               "TminD_v2.0_swiss.lv95"]
        logger.debug("List with source folders created.")
        for s in src:
            dst = Path(self.model_dir, "model", "data_obs")
            logger.debug(f"Symlink src: MeteoSwiss_gridded_products/{s}")
            logger.debug(f"Symlink dst: {dst}")
            try:
                os.symlink(Path(self.data_dir, "MeteoSwiss_gridded_products", s), dst)
                logger.debug(f"Symlink created")
                print(f"Symlink created:\n"
                      f"Source: MeteoSwiss_gridded_products/{s}\n"
                      f"Destination: {dst}")
            except FileExistsError:
                logger.exception("Error creating symlink: File already exists")
        src = Path(self.data_dir, "Discharge", "BroPay_Q_2034_daily.rvt")
        logger.debug("Source Path created.")
        logger.debug(f"Symlink src: {src}")
        logger.debug(f"Symlink dst: {dst}")
        try:
            os.symlink(src, dst)
            logger.debug(f"Symlink created")
            print(f"Symlink created:\n"
                  f"Source: {src}\n"
                  f"Destination: {dst}")
        except FileExistsError:
            logger.exception("Error creating symlink: File already exists")

    def write_rvx(self, ostrich_template: bool = False, raven_template: bool = True, rvx_type: str = "rvi"):
        """Write .rvX file for Raven and/or Ostrich

        :param rvx_type: Suffix of the Raven file to be written (rvi,rvp,rvc,rvt or rvh)
        :param ostrich_template: Set True if Ostrich template should be generated.
        :param raven_template: Set True if Raven template should be generated.
        :type raven_template: bool
        :type ostrich_template: bool

        """
        assert rvx_type in config.variables.raven_filetypes, f"Raven suffix .{rvx_type} is not in list of accepted suffixes."
        # TODO: Implement discharge and forcings time series.
        logger.debug("Starting if-tree for template type...")
        if raven_template is True:
            logger.debug("Variable raven_template is True...")
            logger.debug(f"Trying to call rr.write_rvx function to create .{rvx_type} for Raven...")
            rr.write_rvx(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                         catchment=self.catchment, params=self.default_params, template_type="Raven",
                         rvx_type=rvx_type)
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
                         catchment=self.catchment, params=self.default_params, ost_in=ost_in,
                         save_best=save_best,
                         ost_raven=ost_raven)
        logger.debug(f"ostIn.txt for Raven created by rr.write_rvx function")

    def create_netcdf(self, merge=True):
        """Create the netCDF files for the chosen catchment

        """
        rpe.netcdf_clipper_multi()
        rpe.netcdf_to_dataset()
        if merge:
            pass
