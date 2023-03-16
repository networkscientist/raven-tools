from pathlib import Path

import raven_tools as rt

models_by_name = rt.config.variables.supported_models
# catchments_by_id = list(rt.config.variables.catchments.keys())
catchments_by_id = [key for key in rt.config.variables.catchments]

for c in catchments_by_id:
    for f in rt.config.variables.forcings_dirs:
        model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=c, start_year=1981, end_year=2020)
        model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
        model_instance.clip_netcdf(forcing_prefix=f)
        model_instance.create_grid_weights(forcing_name=f)
print("Finished")
