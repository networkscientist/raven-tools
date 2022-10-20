import os
from pathlib import Path

import yaml


class RavenModel:
    """

    """

    def __init__(self, modeltype = "modeltype", catchment = "catchment"):
        assert isinstance(modeltype, str), f"modeltype expected a string, got {type(modeltype)} instead"
        assert isinstance(catchment, str), f"catchment expected a string, got {type(modeltype)} instead"
        with open("raven_tools/config.yaml", "r") as f:
            config = yaml.load(f, Loader=yaml.FullLoader)
        self.modeltype = modeltype
        self.rootdir = Path(os.getcwd(), "RAVEN")
        self.catchment = catchment
        self.modeldir = Path(self.rootdir, "models", self.catchment, self.modeltype)
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
            Path(self.modeldir).mkdir(parents=True)
            print(f"Model directory created: {Path(self.modeldir)}")
        except FileExistsError:
            print("Directory already exists...")
