import os
from pathlib import Path

import yaml


class RavenModel:
    """

    """

    def __init__(self, modeltype="modeltype", catchment="catchment"):
        logger.debug("Starting __init__...")
        assert isinstance(modeltype, str), f"modeltype expected a string, got {type(modeltype)} instead"
        assert isinstance(catchment, str), f"catchment expected a string, got {type(modeltype)} instead"
        logger.debug("Trying to open config.yaml...")
        with open("raven_tools/new_model_config.yaml", "r") as f:
            self.config = yaml.load(f, Loader=yaml.FullLoader)
        logger.debug("config.yaml loaded.")
        logger.debug("Trying to set self.X variables...")
        self.modeltype = modeltype
        self.rootdir = Path(os.getcwd(), "RAVEN")
        self.catchment = catchment
        self.modeldir = Path(self.rootdir, "models", self.catchment, self.modeltype)
        logger.debug("Self.X variables set.")
        self.dirs = [
            Path(self.modeldir,"model"),
            Path(self.modeldir,"model","output")
        ]

        self.
    @property
    def modeltype(self):
        print("Getting model type...")
        return (self._modeltype)

    @modeltype.setter
    def modeltype(self, value):
        print("Setting model type...")
        self._modeltype = value

    @property
    def rootdir(self):
        return (self._rootdir)

    @rootdir.setter
    def rootdir(self, value):
        self._rootdir = value

    @property
    def modeldir(self):
        return self._modeldir

    @modeldir.setter
    def modeldir(self, value):
        self._modeldir = Path(self.rootdir, "models", self.catchment, self.modeltype)

    @property
    def catchment(self):
        return self._catchment

    @catchment.setter
    def catchment(self, value):
        self._catchment = value

    def create_dirs(self):
        try:
            for f in self.dirs:
                f.mkdir(parents=True)
                print(f"Directory created: {f}")
        except FileExistsError:
            print("Directory already exists...")
