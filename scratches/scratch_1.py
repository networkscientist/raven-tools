import ravenpy.models
from ravenpy.models import GR4JCN_OST
import os
from pathlib import Path
import datetime as dt

model = GR4JCN_OST(workdir="/media/mainman/Data/RAVEN/testmodels/GR4J_tmp")
model = ravenpy.models.Raven(workdir="/media/mainman/Data/RAVEN/testmodels/GR4J_tmp")
original_dir = Path("/media/mainman/Data/RAVEN/testmodels/GR4J/model")
model_name = "raven_broye_gr4j"
file_types = [
    ".rvi",
    ".rvc",
    ".rvh",
    ".rvp",
    ".rvt"
]
file_names = []
for f in file_types:
    file_names.append(Path(original_dir,(model_name+f)))
    
forcings = [
    "data_obs/RhiresD_v2.0_swiss.lv95/merged/RhiresD_ch01h.swiss.lv95_198101010000_201912310000_clipped.nc",
    "data_obs/TabsD_v2.0_swiss.lv95/merged/TabsD_ch01r.swiss.lv95_198101010000_201912310000_clipped.nc",
    "data_obs/TmaxD_v2.0_swiss.lv95/merged/TmaxD_ch01r.swiss.lv95_198101010000_201912310000_clipped.nc",
    "data_obs/TminD_v2.0_swiss.lv95/merged/TminD_ch01r.swiss.lv95_198101010000_201912310000_clipped.nc",    
]

forcings = [Path(original_dir,f) for f in forcings]
model.configure()
model.configure([
    file_names[0],
    file_names[1],
    file_names[2],
    file_names[3],
    file_names[4],
])
model(forcings)
ravenpy.config.commands.RavenCommand.
model(
    ts=forcings,
    start_date=dt.datetime(1981,1,1),
    end_date = dt.datetime(2019,12,31)
)
model.config.commands.
