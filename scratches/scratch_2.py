# from ravenpy.models import Ostrich
# from ravenpy.models import Raven
# from ravenpy.models import GR4JCN
import ravenpy
from ravenpy.utilities.testdata import get_file
from pathlib import Path
import datetime as dt
original_dir = "/media/mainman/Data/RAVEN/tutorial/Salmon_GR4J"
model_name = "raven-gr4j-salmon"
file_types = [
    ".rvi",
    ".rvt",
    ".rvh",
    ".rvc",
    ".rvp"
]
rv_files = []
for f in file_types:
    rv_files.append(str(Path(original_dir,(model_name+f))))
# forcings = [
#     "data_obs/RhiresD_v2.0_swiss.lv95/merged/RhiresD_ch01h.swiss.lv95_198101010000_201912310000_clipped.nc",
#     "data_obs/TabsD_v2.0_swiss.lv95/merged/TabsD_ch01r.swiss.lv95_198101010000_201912310000_clipped.nc",
#     "data_obs/TmaxD_v2.0_swiss.lv95/merged/TmaxD_ch01r.swiss.lv95_198101010000_201912310000_clipped.nc",
#     "data_obs/TminD_v2.0_swiss.lv95/merged/TminD_ch01r.swiss.lv95_198101010000_201912310000_clipped.nc",
# ]
# forcings = [Path(original_dir,f) for f in forcings]
forcings = Path(original_dir, "data_obs")
workdir = "/media/mainman/Data/RAVEN/testmodels/Broye/tmp/testrun/"
identifier = "broye-gr4j"
description = "GR4J for the Broye catchment"
# raven_broye_gr4j = ravenpy.models.Raven(workdir=workdir,identifier=identifier,description=description)
#
# # config = dict(
# #     start_date=dt.datetime(2000,1,1)
# # )
# raven_broye_gr4j.configure(rv_files)
# raven_broye_gr4j(ts=forcings)
#
# raven_broye_gr4j.outputs["rv_config"]
#
# raven_broye_gr4j(forcings)
# raven_broye_gr4j.setup()
# ostrich_broye_gr4j = ravenpy.models.Ostrich(workdir=workdir)

# ---------------------------------------------------------------------------------------
model = ravenpy.models.GR4JCN(workdir=workdir,identifier=identifier,description=description)
model.configure(rv_files)
model.setup()
model.run(ts=forcings, start_date=dt.datetime(2000, 1, 1))
model(
    forcings,
    start_date=dt.datetime(2000, 1, 1),
    end_date=dt.datetime(2002, 1, 1),
    area=4250.6,
    elevation=843.0,
    latitude=54.4848,
    longitude=-123.3659,
    params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
)

model = ravenpy.models.GR4JCN_OST(workdir="/media/mainman/Data/RAVEN/testmodels/Broye/tmp/testrun/")
# model(ts='/media/mainman/Data/RAVEN/testmodels/Broye/GR4J_ravenpy_template/model/data_obs')
model.configure(config)
model.ost2raven()
model = ravenpy.models.Ostrich(workdir="/media/mainman/Data/RAVEN/testmodels/Broye/tmp/testrun/")
model.run(ts='/media/mainman/Data/RAVEN/testmodels/Broye/GR4J_ravenpy_template/model/data_obs')
model.configure(config)
model(ts='/media/mainman/Data/RAVEN/testmodels/Broye/GR4J_ravenpy_template/model/data_obs')
model(start_date=dt.datetime(1981,1,1))
model.setup()

model(forcings)

forcing = get_file(
    "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",github_url="https://github.com/Ouranosinc/raven-testdata"
)

display(tuple(config), forcing)

from ravenpy.models import Raven

model = Raven(workdir="/media/mainman/Data/RAVEN/testmodels/Broye/tmp/testrun")
model.configure(config)
model(forcing)
model.hydrograph
model.q_sim.plot()

import datetime as dt
from ravenpy.models import GR4JCN, HMETS, MOHYSE, HBVEC

model = GR4JCN()
model(
    forcing,
    start_date=dt.datetime(2000, 1, 1),
    end_date=dt.datetime(2002, 1, 1),
    area=4250.6,
    elevation=843.0,
    latitude=54.4848,
    longitude=-123.3659,
    params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
)


from ravenpy.utilities.nb_graphs import hydrographs
from ravenpy.utilities.graphs import hydrograph


