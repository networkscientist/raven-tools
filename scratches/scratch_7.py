from pathlib import Path

import raven_tools as rt

catchments = ["Ticino",
              "Broye",
              "Thur",
              "Massa",
              "Weisse_Luetschine",
              "Dischmabach"]
# catchments = ["Broye"]

for c in catchments:
    for f in rt.config.variables.forcings_dirs:
        model_instance = rt.model.raven_model.RavenModel(catchment=c, start_year=1981, end_year=2020)
        model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
        model_instance.clip_netcdf(forcing_prefix=f)
        model_instance.create_grid_weights(forcing_name=f)
print("Finished")
