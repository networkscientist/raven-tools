"""
Work with a Raven class.
"""
import functools
import logging
import os
import shutil
from pathlib import Path

import dask_geopandas
import geopandas as gpd
import pandas as pd

logger = logging.getLogger(__name__)
logger.debug("Logging from raven_model to console started")
from raven_tools import config
from raven_tools.processing import raven_run as rr
from raven_tools.processing import raven_preprocess as rpe
from pyproj import Transformer
import csv
import subprocess
import rioxarray as rxr
import math
from typing import List


class RavenModel:
    """Class that allows various operations on a selected model. There are methods to create the suggested directory
    structure, to create .rvX RAVEN configuration files and to write Ostrich configuration files.

    Attributes:
        bbox_filepath: Path to the bounding box file.
        catchment_ch_id: A string with the catchment id in the form of 'CH-XXXX'.
        supported_models:
        raven_filetypes:
        conf:
        ctm_info:
        gauge_lat:
        gauge_lon:
        stream_name:
        model_type:
        start_year:
        end_year:
        root_dir:
        gauge_id:
        gauge_short_code:
        station.elevation:
        attribute_csv:
        model_dir:
        model_sub_dir:
        dirs:
        data_dir:
        meteo_dir:
        start_year:
        end_year:
        glaciation_ratio:
        glacier_alti:
        default_params:
        raven_exe_path:
        ost_exe_path:
    """

    def __init__(self, model_type: str = "GR4J", catchment_ch_id: str = "CH-0057",
                 start_year: str = 1981, end_year: str = 2020):
        """Initializes the instance.

        Args:
            model_type: Name of model_type (GR4J, HYMOD, HMETS, HBV or MOHYSE)
            catchment_ch_id: Catchment id in the form of 'CH-XXXX'.
            start_year: First year to be used in calculations.
            end_year: Last year to be used in calculations.
        """
        self.bbox_filepath = Path(os.getcwd())
        logger.debug(f"Starting __init__ of {__name__}...")
        assert isinstance(model_type, str), f"model_type expected a string, got {type(model_type)} instead"
        assert isinstance(catchment_ch_id, str), f"catchment_id expected a string, got {type(catchment_ch_id)} instead"
        self.catchment_ch_id = catchment_ch_id
        logger.info(f"Self.catchment_ch_id set to {self.catchment_ch_id}.")
        self.supported_models = config.variables.supported_models
        self.raven_filetypes = config.variables.raven_filetypes
        self.forcing_type_names = config.variables.forcings_dirs
        assert model_type in self.supported_models, f"model_type expected GR4J, HYMOD, HMETS, HBV or MOHYSE, got {model_type} instead "
        logger.debug(f"CWD: {os.getcwd()}")
        logger.debug("Trying to open project_config.yaml...")
        try:
            self.conf = config.variables.project_config
        except:
            logger.debug("Error getting project_config file from __init__.py")

        try:
            self.ctm_info = config.variables.catchments[self.catchment_ch_id]['ID']
            self.gauge_lat, self.gauge_lon = ch1903_to_wgs84(config.variables.catchments[self.catchment_ch_id]["lat"],
                                                             config.variables.catchments[self.catchment_ch_id]["lon"])
        except:
            logger.exception("Error getting catchments info from __init__.py")
        self.stream_name = config.variables.catchments[self.catchment_ch_id]['stream_name']
        logger.debug("project_config.yaml loaded.")
        self.model_type = model_type
        self.start_year = start_year
        self.end_year = end_year
        self.root_dir = Path(os.path.join(os.getcwd(), Path("RAVEN")))
        self.gauge_id = config.variables.catchments[self.catchment_ch_id]['ID']
        self.gauge_short_code = config.variables.catchments[self.catchment_ch_id]['short_code']
        logger.info(f"Self.gauge_id set to {self.gauge_id}.")
        logger.info(f"Self.gauge_short_code set to {self.gauge_short_code}.")
        self.station_elevation = config.variables.catchments[self.catchment_ch_id]['station_elevation']
        logger.debug(f"Self.station_elevation set to {self.station_elevation}.")
        self.attribute_csv = f"{self.catchment_ch_id}_attributes.csv"
        self.model_dir = Path(self.root_dir, "models", self.catchment_ch_id, self.model_type)
        self.model_sub_dir = self.conf['ModelSubDir']
        self.dirs: list[Path] = [
            Path(self.model_dir, ""),
            Path(self.model_dir, self.model_sub_dir, "", "output"),
            Path(self.model_dir, self.model_sub_dir, "", "data_obs")
        ]
        self.data_dir: Path = Path(self.root_dir, self.conf['DataDir'])
        self.meteo_dir = self.conf['MeteoSubDir']
        self.start_year = self.conf['StartYear']
        self.end_year = self.conf['EndYear']
        logger.debug("Self.X variables set.")
        logger.debug(f"__init__ of {__name__} finished...")
        self.glaciation_ratio: float = 0.0
        self.glacier_alti: float = 0.0
        try:
            with open(Path(self.data_dir, "glaciers", f"glaciation_ratio_{self.catchment_ch_id}.txt"),
                      mode='r') as csv_file:
                csv_reader = csv.DictReader(csv_file, delimiter=';')
                for row in csv_reader:
                    self.glaciation_ratio = float(row['Gla_Ratio'])
                    self.glacier_alti = float(row['Gla_Alti'])
        except:
            logger.exception("Error reading glacier data")
        try:
            self.default_params = config.variables.default_params
        except:
            logger.exception("Error getting default_params file from __init__.py")
        self.raven_exe_path: Path = Path(self.conf['RavenExePath'])
        self.ost_exe_path: Path = Path(self.conf['OstrichExePath'])
        if self.model_type == "HBV":
            self.grid_weights_file_path = Path(self.data_dir, "MeteoSwiss_gridded_products", self.forcing_type_names[0],
                                               "out", f"grid_weights_{self.catchment_ch_id}_hbv.txt")
        else:
            self.grid_weights_file_path = Path(self.data_dir, "MeteoSwiss_gridded_products", self.forcing_type_names[0],
                                               "out", f"grid_weights_{self.catchment_ch_id}.txt")
        logger.debug(f"Raven exe path set: {self.raven_exe_path}")

    def __getitem__(self, item):
        print(type(item), item)

    @property
    def model_sub_dir(self) -> str:
        """Sub-directory in which the RAVEN model files are located."""
        assert isinstance(self._model_sub_dir,
                          str), f"model_sub_dir should be str, is type {type(self._model_sub_dir)} instead."
        return self._model_sub_dir

    @model_sub_dir.setter
    def model_sub_dir(self, value: str):
        assert isinstance(value, str), f"model_sub_dir should be str, is type {type(self._model_sub_dir)} instead."
        self._model_sub_dir = value

    @property
    def data_dir(self) -> Path:
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
    def gauge_id(self) -> int:
        """Returns gauge_id."""
        assert isinstance(self._gauge_id, int), f"gauge_id should be int, is type {type(self._gauge_id)} instead."
        return self._gauge_id

    @gauge_id.setter
    def gauge_id(self, value: int):
        assert isinstance(value, int), f"gauge_id should be int, is type {type(value)} instead."
        self._gauge_id = value

    @property
    def root_dir(self) -> Path:
        """Get root directory that contains the 'RAVEN' folder.

        Args:

        Returns:
            self._root_dir : Path
                Root directory that contains 'RAVEN' folder
        """
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
    def stream_name(self) -> str:
        """Return catchment name

        :return self._catchment: Catchment name
        :rtype: str
        """
        logger.debug("Getting catchment name...")
        return self._catchment

    @stream_name.setter
    def stream_name(self, value: str):
        """Set catchment name.

        :param value: Catchment name
        :type value: str

        """
        logger.debug("Setting catchment name...")
        self._catchment = value

    @property
    def ctm_info(self) -> int:
        """Returns ctm_info."""
        assert isinstance(self._ctm_info, int), f"ctm_info should be dict, is type {type(self._ctm_info)} instead."
        return self._ctm_info

    @ctm_info.setter
    def ctm_info(self, value: int):
        assert isinstance(value, int), f"ctm_info should be dict, is type {type(value)} instead."
        self._ctm_info = value

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
    def gauge_lat(self) -> float:
        """Returns gauge_lat."""
        assert isinstance(self._gauge_lat,
                          float), f"gauge_lat should be float, is type {type(self._gauge_lat)} instead."
        return self._gauge_lat

    @gauge_lat.setter
    def gauge_lat(self, value: float):
        assert isinstance(value, float), f"gauge_lat should be float, is type {type(value)} instead."
        self._gauge_lat = value

    @property
    def gauge_lon(self) -> float:
        """Returns gauge_lon."""
        assert isinstance(self._gauge_lon,
                          float), f"gauge_lon should be float, is type {type(self._gauge_lon)} instead."
        return self._gauge_lon

    @gauge_lon.setter
    def gauge_lon(self, value: float):
        assert isinstance(value, float), f"gauge_lon should be float, is type {type(value)} instead."
        self._gauge_lon = value

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
    def catchment_ch_id(self) -> str:
        """Returns catchment_ch_id."""
        assert isinstance(self._catchment_ch_id,
                          str), f"catchment_ch_id should be str, is type {type(self._catchment_ch_id)} instead."
        return self._catchment_ch_id

    @catchment_ch_id.setter
    def catchment_ch_id(self, value: str):
        assert isinstance(value, str), f"catchment_ch_id should be str, is type {type(value)} instead."
        self._catchment_ch_id = value

    @property
    def glaciation_ratio(self) -> float:
        """Returns glaciation_ratio."""
        assert isinstance(self._glaciation_ratio,
                          float), f"glaciation_ratio should be float, is type {type(self._glaciation_ratio)} instead."
        return self._glaciation_ratio

    @glaciation_ratio.setter
    def glaciation_ratio(self, value: float):
        assert isinstance(value, float), f"glaciation_ratio should be float, is type {type(value)} instead."
        self._glaciation_ratio = value

    @property
    def glacier_alti(self) -> float:
        """Returns glacier_alti."""
        assert isinstance(self._glacier_alti,
                          float), f"glacier_alti should be float, is type {type(self._glacier_alti)} instead."
        return self._glacier_alti

    @glacier_alti.setter
    def glacier_alti(self, value: float):
        assert isinstance(value, float), f"glacier_alti should be float, is type {type(value)} instead."
        self._glacier_alti = value

    @property
    def glacier_slope(self) -> float:
        """Returns glacier_slope."""
        assert isinstance(self._glacier_slope,
                          float), f"glacier_slope should be float, is type {type(self._glacier_slope)} instead."
        return self._glacier_slope

    @glacier_slope.setter
    def glacier_slope(self, value: float):
        assert isinstance(value, float), f"glacier_slope should be float, is type {type(value)} instead."
        self._glacier_slope = value

    @property
    def glacier_aspect(self) -> float:
        """Returns glacier_aspect."""
        assert isinstance(self._glacier_aspect,
                          float), f"glacier_aspect should be float, is type {type(self._glacier_aspect)} instead."
        return self._glacier_aspect

    @glacier_aspect.setter
    def glacier_aspect(self, value: float):
        assert isinstance(value, float), f"glacier_aspect should be float, is type {type(value)} instead."
        self._glacier_aspect = value

    @property
    def non_glacier_slope(self) -> float:
        """Returns non_glacier_slope."""
        assert isinstance(self._non_glacier_slope,
                          float), f"non_glacier_slope should be float, is type {type(self._non_glacier_slope)} instead."
        return self._non_glacier_slope

    @non_glacier_slope.setter
    def non_glacier_slope(self, value: float):
        assert isinstance(value, float), f"non_glacier_slope should be float, is type {type(value)} instead."
        self._non_glacier_slope = value

    @property
    def non_glacier_aspect(self) -> float:
        """Returns non_glacier_aspect."""
        assert isinstance(self._non_glacier_aspect,
                          float), f"non_glacier_aspect should be float, is type {type(self._non_glacier_aspect)} instead."
        return self._non_glacier_aspect

    @non_glacier_aspect.setter
    def non_glacier_aspect(self, value: float):
        assert isinstance(value, float), f"non_glacier_aspect should be float, is type {type(value)} instead."
        self._non_glacier_aspect = value

    @property
    def grid_weights_file_path(self) -> Path:
        """Returns grid_weights_file_path."""
        assert isinstance(self._grid_weights_file_path,
                          Path), f"grid_weights_file_path should be Path, is type {type(self._grid_weights_file_path)} instead."
        return self._grid_weights_file_path

    @grid_weights_file_path.setter
    def grid_weights_file_path(self, value: Path):
        assert isinstance(value, Path), f"grid_weights_file_path should be Path, is type {type(value)} instead."
        self._grid_weights_file_path = value

    @property
    def forcing_type_names(self) -> list[str]:
        """Returns forcing_type_names."""
        assert isinstance(self._forcing_type_names,
                          list), f"forcing_type_names should be List[str], is type {type(self._forcing_type_names)} instead."
        return self._forcing_type_names

    @forcing_type_names.setter
    def forcing_type_names(self, value: List[str]):
        assert isinstance(value, list), f"forcing_type_names should be List[str], is type {type(value)} instead."
        self._forcing_type_names = value

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

    @property
    def gauge_short_code(self) -> str:
        """Returns gauge_short_code."""
        assert isinstance(self._gauge_short_code,
                          str), f"gauge_short_code should be str, is type {type(self._gauge_short_code)} instead."
        return self._gauge_short_code

    @gauge_short_code.setter
    def gauge_short_code(self, value: str):
        assert isinstance(value, str), f"gauge_short_code should be str, is type {type(value)} instead."
        self._gauge_short_code = value

    @property
    def station_elevation(self) -> str:
        """Returns station_elevation."""
        assert isinstance(self._station_elevation,
                          str), f"station_elevation should be str, is type {type(self._station_elevation)} instead."
        return self._station_elevation

    @station_elevation.setter
    def station_elevation(self, value: str):
        assert isinstance(value, str), f"station_elevation should be str, is type {type(value)} instead."
        self._station_elevation = value

    @property
    def bbox_filepath(self) -> Path:
        """Returns bbox_filepath."""
        assert isinstance(self._bbox_filepath,
                          Path), f"bbox_filepath should be Path, is type {type(self._bbox_filepath)} instead."
        return self._bbox_filepath

    @bbox_filepath.setter
    def bbox_filepath(self, value: Path):
        assert isinstance(value, Path), f"bbox_filepath should be Path, is type {type(value)} instead."
        self._bbox_filepath = value

    @property
    def ost_exe_path(self) -> Path:
        """Returns ost_exe_path."""
        assert isinstance(self._ost_exe_path,
                          Path), f"ost_exe_path should be Path, is type {type(self._ost_exe_path)} instead."
        return self._ost_exe_path

    @ost_exe_path.setter
    def ost_exe_path(self, value: Path):
        assert isinstance(value, Path), f"ost_exe_path should be Path, is type {type(value)} instead."
        self._ost_exe_path = value

    def create_symlinks(self, forcings: bool = True, discharge: bool = True, raven_executable: bool = True,
                        ostrich_executable: bool = True, rvx_files: bool = True, raven_diag: bool = True, delete=False):
        from raven_tools.processing import raven_diag
        if forcings:
            src = ["RhiresD_v2.0_swiss.lv95",
                   "SrelD_v2.0_swiss.lv95",
                   "TabsD_v2.0_swiss.lv95",
                   "TmaxD_v2.0_swiss.lv95",
                   "TminD_v2.0_swiss.lv95"]
            logger.debug("List with source folders created.")
            for s in src:
                try:
                    src_path = Path("/storage/homefs/pz09y074/raven_master_files/RAVEN", "data",
                                    "MeteoSwiss_gridded_products", s)
                    dst = Path(self.model_dir, self.model_sub_dir, "data_obs", s)
                    logger.debug(f"Symlink src: {src_path}")
                    logger.debug(f"Symlink dst: {dst}")
                    os.symlink(src_path, dst, target_is_directory=True)
                    logger.debug(f"Symlink created")
                except FileExistsError:
                    logger.exception("Error creating symlink: File already exists")
        if discharge:
            discharge_filename = f"{self.gauge_short_code}_Q_{self.gauge_id}_daily.rvt"
            src = Path(self.data_dir, "Discharge", discharge_filename)
            dst = Path(self.model_dir, self.model_sub_dir, "data_obs", discharge_filename)
            logger.debug(f"Symlink src: {src}")
            logger.debug(f"Symlink dst: {dst}")
            if delete:
                try:
                    Path.unlink(dst)
                except:
                    logger.exception("There has been an exception unlinking...")
            else:
                try:
                    rcode = subprocess.call(['ln', "-sn", src, dst])
                    logger.debug(f"extent.sh executed with return code: {rcode}")
                    os.symlink(src, dst)
                    logger.debug(f"Symlink src: {src}")
                    logger.debug(f"Symlink dst: {dst}")
                except FileExistsError:
                    logger.exception("Error creating symlink: File already exists")

        if raven_executable:
            src = Path(self.raven_exe_path)
            dst = Path(self.model_dir, self.model_sub_dir, "Raven.exe")
            logger.debug(f"Symlink src: {src}")
            logger.debug(f"Symlink dst: {dst}")
            try:
                os.symlink(src, dst)
                logger.debug(f"Symlink created")
            except FileExistsError:
                logger.exception("Error creating symlink: File already exists")
                pass

        if ostrich_executable:
            src = Path(self.ost_exe_path)
            dst = Path(self.model_dir, "OstrichMPI")
            logger.debug(f"Symlink src: {src}")
            logger.debug(f"Symlink dst: {dst}")
            try:
                os.symlink(src, dst)
                logger.debug(f"Symlink created")
            except FileExistsError:
                logger.exception("Error creating symlink: File already exists")
                pass

        if rvx_files:
            for s in config.variables.raven_filetypes:
                src = Path(self.model_dir, f"{self.catchment_ch_id}_{self.model_type}.{s}")
                dst = Path(self.model_dir, self.model_sub_dir, f"{self.catchment_ch_id}_{self.model_type}.{s}")
                logger.debug(f"Symlink src: {src}")
                logger.debug(f"Symlink dst: {dst}")
                try:
                    os.symlink(src, dst)
                    logger.info(".rvX symlink created.")
                except FileExistsError:
                    logger.exception("Error creating symlink: File already exists.")
                    pass

        if raven_diag:
            shutil.copy(raven_diag.__file__, Path(self.model_dir, self.model_sub_dir, "output", "raven_diag.py"))

    def write_rvx(self, ostrich_template: bool = False, raven_template: bool = True, rvx_type: str = "rvi",
                  glacier_module: bool = False):
        """Write .rvX file for Raven and/or Ostrich

        :param rvx_type: Suffix of the Raven file to be written (rvi,rvp,rvc,rvt or rvh)
        :param ostrich_template: Set True if Ostrich template should be generated.
        :param raven_template: Set True if Raven template should be generated.
        :type raven_template: bool
        :type ostrich_template: bool

        """
        hru_info_dict = {}
        assert rvx_type in config.variables.raven_filetypes, f"Raven suffix .{rvx_type} is not in list of accepted suffixes."
        try:
            hru_info_df = pd.read_csv(Path(self.data_dir, "Catchment", "hru_info.csv"), na_values='-')
            hru_info_df = hru_info_df[hru_info_df['Ctm'] == self.catchment_ch_id]
            hru_info_dict = hru_info_df.to_dict('records')[0]
        except FileNotFoundError:
            logger.exception("Glacier properties file could not be found. Maybe you need to create it first.")
        if raven_template is True:
            logger.debug("Variable raven_template is True...")
            logger.debug(f"Trying to call rr.write_rvx function to create .{rvx_type} for Raven...")
            rr.write_rvx(catchment_ch_id=self.catchment_ch_id, hru_info=hru_info_dict, model_dir="models",
                         model_type=self.model_type,
                         project_dir=self.root_dir, params=self.default_params, template_type="Raven",
                         rvx_type=rvx_type, start_year=self.start_year, end_year=self.end_year,
                         glacier_module=glacier_module,
                         data_dir=self.data_dir)
            logger.debug(f".{rvx_type} for Raven created by rr.write_rvx function")
        if ostrich_template is True:
            logger.debug("Variable ostrich_template is True...")
            logger.debug(f"Trying to call rr.write_rvx function to create .{rvx_type}.tpl for Ostrich...")
            rr.write_rvx(catchment_ch_id=self.catchment_ch_id, hru_info=hru_info_dict, model_dir="models",
                         model_type=self.model_type,
                         project_dir=self.root_dir, params=self.default_params, template_type="Ostrich",
                         rvx_type=rvx_type, start_year=self.start_year, end_year=self.end_year,
                         glacier_module=glacier_module,
                         data_dir=self.data_dir)
            logger.debug(f".{rvx_type}.tpl for Ostrich created by rr.write_rvx function")
        if ostrich_template is False and raven_template is False:
            logger.debug("Variables ostrich_template and raven_template set to False.")
            print("You have not selected a template type...")
            logger.debug("No template file needed to be written by function rr.write_rvx.")
            pass

    def write_ost(self,
                  ost_in=True,
                  save_best=True,
                  ost_raven=True,
                  ost_mpi_script: bool = True,
                  max_iterations: int = 500):
        """Write Ostrich files ostIn.txt, save_best.sh and ost-raven.sh

        :param ost_raven:
        :param ost_in:
        :param save_best:
        """

        logger.debug("Variable raven_template is True...")
        logger.debug(f"Trying to call rr.write_ostrich_files() function to create ostrich input files\
            for Raven...")
        rr.write_ostrich(model_dir="models", model_type=self.model_type, project_dir=self.root_dir,
                         catchment_name=self.stream_name, catchment_ch_id=self.catchment_ch_id,
                         params=self.default_params, ost_in=ost_in,
                         save_best=save_best,
                         ost_raven=ost_raven,
                         ost_mpi_script=ost_mpi_script,
                         max_iterations=max_iterations)
        logger.debug(f"ostIn.txt for Raven created by rr.write_rvx function")

    def clip_netcdf(self, forcing_prefix="TabsD_v2.0_swiss.lv95"):
        """Create the netCDF files for the chosen catchment

        """
        extent_file_path = Path(self.data_dir, "Catchment", f"{self.catchment_ch_id}.shp")
        netcdf_file_path = Path(self.data_dir, "MeteoSwiss_gridded_products", forcing_prefix, "out",
                                f"{forcing_prefix}_{self.start_year}01010000_{self.end_year}12310000.nc")
        try:
            rpe.netcdf_clipper(netcdf_file_path=netcdf_file_path, extent_file_path=extent_file_path)
            logger.info(f"netcdf_file_path = {netcdf_file_path}")
        except:
            logger.exception(f"Error creating netCDF file {netcdf_file_path}")

    def merge_netcdf(self, forcing_prefix):
        forcing_dir: Path = Path(self.data_dir, "MeteoSwiss_gridded_products", forcing_prefix, "original_files")
        rpe.nc_merge(start_year=self.start_year, end_year=self.end_year, forcing_dir=forcing_dir,
                     forcing_prefix=forcing_prefix)

    def create_bbox(self):
        rpe.create_bbox(extent_shape_file_path=Path(self.data_dir, "Catchment", "reproject_2056",
                                                    f"{self.catchment_ch_id}.shp"),
                        bb_file_path=Path(self.data_dir, "Catchment", f"{self.catchment_ch_id}_bbox.shp"),
                        create_bbox_shp=True)

    def write_rvt(self, ostrich_template: bool = True, raven_template: bool = True):
        hru_info_dict = {}
        try:
            hru_info_df = pd.read_csv(Path(self.data_dir, "Catchment", "hru_info.csv"), na_values='-')
            hru_info_df = hru_info_df[hru_info_df['Ctm'] == self.catchment_ch_id]
            hru_info_dict = hru_info_df.to_dict('records')[0]
        except FileNotFoundError:
            logger.exception("Glacier properties file could not be found. Maybe you need to create it first.")
        write_rvt_partial = functools.partial(rr.write_rvt,
                                              start_year=self.start_year,
                                              end_year=self.end_year,
                                              model_dir=self.model_dir,
                                              model_type=self.model_type,
                                              project_dir=self.root_dir,
                                              catchment_ch_id=self.catchment_ch_id,
                                              gauge_lat=self.gauge_lat,
                                              gauge_lon=self.gauge_lon,
                                              model_sub_dir=self.model_sub_dir,
                                              gauge_short_code=self.gauge_short_code,
                                              station_elevation=self.station_elevation,
                                              catchment_gauge_id=str(self.ctm_info),
                                              params=self.default_params,
                                              hru_info=hru_info_dict,
                                              data_dir=self.data_dir)
        if ostrich_template:
            write_rvt_partial(template_type="Ostrich")
        if raven_template:
            write_rvt_partial(template_type="Raven")

    def camels_to_rvt(self):
        rpe.camels_to_rvt(data_dir=self.data_dir, gauge_id=self.gauge_id,
                          gauge_short_code=self.gauge_short_code, start_date=f"{self.start_year}-01-01",
                          end_date=f"{self.end_year + 1}-01-01")

    def create_grid_weights(self, read_from_file: False, forcing_name="RhiresD_v2.0_swiss.lv95"):
        global area_ratios
        hru_info = pd.read_csv(Path(self.data_dir, "Catchment/hru_info.csv"), sep=",", na_values='-')
        has_glacier = False
        if not math.isnan(hru_info[hru_info['Ctm'] == self.catchment_ch_id]['GlaArea'].iloc[0]):
            has_glacier = True
        # -------------------------------------------------------------------------------------------------------------
        # HBV Part
        # -------------------------------------------------------------------------------------------------------------
        grid_weights_file_path: Path = Path(self.data_dir, "MeteoSwiss_gridded_products", forcing_name, "out",
                                            f"grid_weights_{self.catchment_ch_id}.txt")
        netcdf_file_path = Path(self.data_dir, "MeteoSwiss_gridded_products", forcing_name, "out",
                                f"{forcing_name}_{self.start_year}01010000_{self.end_year}12310000_{self.catchment_ch_id}_clipped.nc")
        catchment_filepath = Path(self.data_dir, "Catchment", "reproject_2056",
                                  f"{config.variables.catchments[self.catchment_ch_id]['catchment_id']}.shp")
        non_gla_short_name, gla_short_name = "non_gla", "gla"
        grid_weights_file_path_hbv = Path(self.data_dir, "MeteoSwiss_gridded_products", forcing_name, "out",
                                          f"grid_weights_{self.catchment_ch_id}_hbv.txt")
        if read_from_file:
            if self.model_type == "HBV":
                weights_gla = pd.read_csv(grid_weights_file_path, sep="\s+", skiprows=7, skipfooter=1,
                                          names=["hru_id", "cell_id", "weight"])
                # weights = pd.concat(weights_glapd.read_csv(Path(self.data_dir, "Catchment", "HBV", f"weights_{self.catchment_ch_id}.txt"),
                #                       sep="\s+", header=None, names=["hru_id", "cell_id", "weight"])
            if not self.model_type == "HBV":
                pass

        if not read_from_file:
            if self.model_type == "HBV":
                grid_weights_dgdf = dask_geopandas.GeoDataFrame
                non_gla_out_path = Path(self.data_dir, "MeteoSwiss_gridded_products", forcing_name, "out",
                                        f"grid_weights_{self.catchment_ch_id}")
                basic_grid = rpe.create_grid(netcdf_filepath=netcdf_file_path,
                                             bounding_box_filename=self.bbox_filepath,
                                             out_path=non_gla_out_path,
                                             forcing_name=forcing_name, start_year=self.start_year).drop(
                    columns=['row', 'col', 'area'])
                if has_glacier:
                    glacier_shape_path = Path(self.data_dir, "glaciers", "SGI_2016_glaciers.shp")
                    grid_weights_dgdf = rpe.elevation_bands(ctm_ch_id=self.catchment_ch_id,
                                                            data_dir=self.data_dir,
                                                            catchment_filepath=catchment_filepath,
                                                            save_to_tif=False,
                                                            glacier_shape_path=glacier_shape_path,
                                                            basic_grid=basic_grid)
                if not has_glacier:
                    grid_weights_dgdf = rpe.elevation_bands(ctm_ch_id=self.catchment_ch_id,
                                                            data_dir=self.data_dir,
                                                            catchment_filepath=catchment_filepath,
                                                            save_to_tif=False, basic_grid=basic_grid)
                # rpe.dict_to_txt(dict=area_ratios, ctm_ch_id=self.catchment_ch_id, data_dir=self.data_dir)

                # band_list = rpe.create_elevation_band_tif_list_hbv(self.data_dir, 'non_glacier',
                #                                                    ctm_ch_id=self.catchment_ch_id)
                band_list = [col_name for col_name in grid_weights_dgdf.columns.to_list() if col_name != "polygons"]
                hru_id_list = list(range(1, len(band_list) + 1))
                # res, area_total = rpe.elevation_band_slope_aspect(data_dir=self.data_dir,
                #                                                   ctm_ch_id=self.catchment_ch_id,
                #                                                   band_list=band_list,
                #                                                   hru_id=hru_id, hru_info_df=hru_info, grid=basic_grid)
                # rpe.write_hru_info(self.data_dir, self.catchment_ch_id, area_total, res)
                weight_sum = 0
                rpe.write_hbv_weights_to_file(grid_weights_file_path_hbv, grid_weights_dgdf)
                slope_aspect_info_df, area_total, info = rpe.elevation_band_slope_aspect(data_dir=self.data_dir,
                                                                                         ctm_ch_id=self.catchment_ch_id,
                                                                                         band_list=band_list,
                                                                                         hru_id=hru_id_list,
                                                                                         hru_info_df=hru_info,
                                                                                         grid=basic_grid,
                                                                                         grid_weights_df=grid_weights_dgdf)

                logger.debug(f"grid weights written to file {grid_weights_file_path_hbv}")
            if not self.model_type == "HBV":
                pass


        # -------------------------------------------------------------------------------------------------------------
        # All other models
        # -------------------------------------------------------------------------------------------------------------
        else:
            if forcing_name == "RhiresD_v2.0_swiss.lv95":
                forcing_suffix = "ch01h.swiss.lv95"
            else:
                forcing_suffix = "ch01r.swiss.lv95"
            non_gla_out_path = Path(self.data_dir, "MeteoSwiss_gridded_products", forcing_name, "out",
                                    f"grid_weights_{self.catchment_ch_id}")
            basic_grid = rpe.create_grid(netcdf_filepath=netcdf_file_path, bounding_box_filename=self.bbox_filepath,
                                         out_path=non_gla_out_path,
                                         forcing_name=forcing_name, start_year=self.start_year).drop(
                columns=['row', 'col', 'area'])
            # ---------------------------------------------------------------------------------------------------------
            # Non-Glacier Part
            # ---------------------------------------------------------------------------------------------------------
            if math.isnan(hru_info[hru_info['Ctm'] == self.catchment_ch_id]['GlaArea'].iloc[0]):
                hru_short_name: str = "non_gla"
                non_gla_extent = gpd.read_file(catchment_filepath).drop(
                    columns=['ezgnr', 'ar_1903', 'watr_nm', 'place'])
                cols = ['id', 'geometry']
                combi_grid = rpe.calc_hru_weight(basic_grid=basic_grid, hru_extent_gdf=non_gla_extent,
                                                 hru_short_name=hru_short_name).drop(columns=cols)
            # ---------------------------------------------------------------------------------------------------------
            # Glacier part
            # ---------------------------------------------------------------------------------------------------------
            else:
                glacier_shape_path = Path(self.data_dir, "glaciers", "SGI_2016_glaciers.shp")
                gla_extent, non_gla_extent = rpe.hru_extent_from_shp(ctm_shp=catchment_filepath,
                                                                     hru_shp=glacier_shape_path)
                gla_extent = gla_extent.dissolve()
                cols = ['id', 'geometry']
                non_gla_grid = rpe.calc_hru_weight(basic_grid=basic_grid, hru_extent_gdf=non_gla_extent,
                                                   hru_short_name=non_gla_short_name).drop(columns=cols)
                cols = ['id', 'geometry', 'polygons']
                gla_grid = rpe.calc_hru_weight(basic_grid=basic_grid, hru_extent_gdf=gla_extent,
                                               hru_short_name=gla_short_name).drop(columns=cols)
                combi_grid: gpd.GeoDataFrame = non_gla_grid.join(other=gla_grid, rsuffix='_g', on='cell_id')
            rpe.write_grid_data_to_file(gdf=combi_grid, grid_weights_file_path=grid_weights_file_path)
            logger.debug(f"grid weights written to file {grid_weights_file_path}")

    def method_name(self, band_list, hru_id, res):
        for (bd, hru) in zip(band_list, hru_id):
            grid = gpd.read_file(Path(self.data_dir,
                                      f"MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/grid_weights_{self.catchment_ch_id}.shp"))
            dem = rxr.open_rasterio(filename=Path(self.data_dir,
                                                  f"DEM/hbv/non_glacier/dem_{self.catchment_ch_id}_{bd}_{str(int(bd) + 99)}.tif"),
                                    masked=True).squeeze()
            dem_df = dem.to_dataframe(name='alti').dropna()
            dem_gdf = gpd.GeoDataFrame(dem_df, crs='epsg:2056',
                                       geometry=gpd.points_from_xy(dem_df.reset_index().x, dem_df.reset_index().y))
            dem_gdf.reset_index(inplace=True)
            centroid = dem_gdf.dissolve().centroid
            ov = rpe.create_overlay(grd=grid, ctm_gdf=dem_gdf)
            cou = ov.value_counts('cell_id')
            cou = cou.to_frame().reset_index()
            dem_aspect = rxr.open_rasterio(
                filename=Path(self.data_dir,
                              f"DEM/hbv/aspects/dem_{self.catchment_ch_id}_{bd}_{str(int(bd) + 99)}_aspects.tif"),
                masked=True).squeeze()
            dem_slope = rxr.open_rasterio(
                filename=Path(self.data_dir,
                              f"DEM/hbv/slopes/dem_{self.catchment_ch_id}_{bd}_{str(int(bd) + 99)}_slopes.tif"),
                masked=True).squeeze()
            cou['Ctm'] = self.catchment_ch_id
            cou['Band'] = bd
            cou['hru_id'] = hru
            cou['aspect_mean'] = rpe.dem_mean_rasterio(dem_aspect)
            cou['slope_mean'] = rpe.dem_mean_rasterio(dem_slope)
            cou['alti_mean'] = rpe.dem_mean_rasterio(dem)
            cou['centroid_lat'], cou['centroid_lon'] = rpe.weighted_centroid_raster_points(
                raster_points_gdf=dem_gdf)
            lat, lon = rpe.crs_old_to_new(lat_old=dem_gdf.to_crs(epsg=32632).y.mean(),
                                          lon_old=dem_gdf.to_crs(epsg=32632).y.mean(), epsg_old=32632,
                                          epsg_new=4326)
            res = pd.concat([res, cou], ignore_index=True)

            # res = ov.value_counts('cell_id')
        return res

    def glacier_ratio(self, dem_tif_filenames):
        catchment_filepath: Path = Path(self.data_dir, "Catchment", "reproject_2056",
                                        f"{config.variables.catchments[self.catchment_ch_id]['catchment_id']}.shp")
        glacier_extent_filepath: Path = Path(self.data_dir, "glaciers", "SGI_2016_glaciers.shp")
        # glaciation_ratio_height = {
        #     "glaciation_ratio":
        #         [],
        #     "glaciation_height":
        #         []
        # }
        dem_filepaths = [Path(self.data_dir, "DEM", f) for f in dem_tif_filenames]
        try:
            glaciation_ratio, glacier_height, non_glaciation_ratio, non_gla_height, glaciated_centroid, non_glaciated_centroid = rpe.glaciation_ratio_height(
                catchment_filepath=catchment_filepath,
                glacier_shape_path=glacier_extent_filepath,
                dem_filepaths=dem_filepaths,
                dem_out_path=Path(self.data_dir, "DEM"),
                ctm_ch_id=self.catchment_ch_id)
            # glaciation_ratio_height["glaciation_ratio"].append(glaciation_ratio)
            # glaciation_ratio_height["glaciation_height"].append(glacier_height)
            return glaciation_ratio, glacier_height, non_glaciation_ratio, non_gla_height, glaciated_centroid, non_glaciated_centroid
        except ValueError:
            logger.exception("Error clipping DEM")
            pass

    def hru_dem_props(self, glacier: bool = False):
        columns = ["NonGlaArea",
                   "GlaArea",
                   "NonGlaAlti",
                   "GlaAlti",
                   "NonGlaLat",
                   "GlaLat",
                   "NonGlaLon",
                   "GlaLon",
                   "NonGlaAspect",
                   "GlaAspect",
                   "NonGlaSlope",
                   "GlaSlope"
                   ]
        catchment_filepath: Path = Path(self.data_dir, "Catchment", "reproject_2056",
                                        f"{config.variables.catchments[self.catchment_ch_id]['catchment_id']}.shp")
        glacier_extent_filepath: Path = Path(self.data_dir, "glaciers", "SGI_2016_glaciers.shp")
        aspect_slope = pd.DataFrame(columns=columns)
        aspect_slope.rename_axis("Ctm", inplace=True)
        ctm_glaciation_gdf, ctm_non_glaciation_gdf = rpe.hru_extent_from_shp(ctm_shp=catchment_filepath,
                                                                             hru_shp=glacier_extent_filepath)
        if glacier:
            filepath_aspect = Path(self.data_dir, "DEM", "dem_clipped", "aspects",
                                   f"dem_{self.catchment_ch_id}_glacier_aspects.tif")
            filepath_slope = Path(self.data_dir, "DEM", "dem_clipped", "slopes",
                                  f"dem_{self.catchment_ch_id}_glacier_slopes.tif")
            filepath_hru = Path(self.data_dir, "DEM", "dem_clipped", f"dem_{self.catchment_ch_id}_glacier.tif")
            aspect_slope.loc[self.catchment_ch_id, f"GlaAspect"] = rpe.dem_mean(filepath_aspect)
            aspect_slope.loc[self.catchment_ch_id, f"GlaSlope"] = rpe.dem_mean(filepath_slope)
            aspect_slope.loc[self.catchment_ch_id, f"GlaAlti"] = rpe.dem_mean(filepath_hru)
            ctm_total_area, non_glaciation_area, glaciation_area = rpe.area_ratio(
                catchment_filepath=catchment_filepath, glacier_shape_path=glacier_extent_filepath)
            aspect_slope.loc[self.catchment_ch_id, f"GlaArea"] = glaciation_area

            gla_lat, gla_lon = rpe.weighted_centroid_vector(ctm_glaciation_gdf)
            aspect_slope.loc[self.catchment_ch_id, f"GlaLat"] = gla_lat
            aspect_slope.loc[self.catchment_ch_id, f"GlaLon"] = gla_lon
        else:
            ctm_total_area, non_glaciation_area = rpe.area_ratio(catchment_filepath=catchment_filepath,
                                                                 glacier_shape_path=glacier_extent_filepath)

        aspect_slope.loc[self.catchment_ch_id, f"NonGlaArea"] = non_glaciation_area
        filepath_aspect = Path(self.data_dir, "DEM", "dem_clipped", "aspects",
                               f"dem_{self.catchment_ch_id}_aspects.tif")
        filepath_slope = Path(self.data_dir, "DEM", "dem_clipped", "slopes",
                              f"dem_{self.catchment_ch_id}_slopes.tif")
        filepath_hru = Path(self.data_dir, "DEM", "dem_clipped", f"dem_{self.catchment_ch_id}.tif")
        aspect_slope.loc[self.catchment_ch_id, f"NonGlaAspect"] = rpe.dem_mean(filepath_aspect)
        aspect_slope.loc[self.catchment_ch_id, f"NonGlaSlope"] = rpe.dem_mean(filepath_slope)

        aspect_slope.loc[self.catchment_ch_id, f"NonGlaAlti"] = rpe.dem_mean(filepath_hru)
        non_gla_lat, non_gla_lon = rpe.weighted_centroid_vector(ctm_non_glaciation_gdf)
        aspect_slope.loc[self.catchment_ch_id, f"NonGlaLat"] = non_gla_lat
        aspect_slope.loc[self.catchment_ch_id, f"NonGlaLon"] = non_gla_lon
        aspect_slope = aspect_slope[columns]
        return aspect_slope


def ch1903_to_wgs84(lat_1903, lon_1903):
    transformer = Transformer.from_crs("EPSG:21781", "EPSG:4326")
    lat_wgs84, lon_wgs84 = transformer.transform(lat_1903, lon_1903)
    return lat_wgs84, lon_wgs84
