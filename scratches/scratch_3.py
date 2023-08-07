import os
from pathlib import Path

import raven_tools as rt

rt_path = Path(os.path.dirname(rt.__file__)).parent.absolute()

os.chdir(rt_path)

suffix = [
    "rvi", "rvh", "rvp", "rvc"
]
models_by_name = rt.config.variables.supported_models
catchments_by_id = [key for key in rt.config.variables.catchments]
models_by_name = ["HBV"]
catchments_by_id = ["CH-0053"]
start_year = 1981
end_year = 2020

# Do the following snippet for each forcing_dir
# ---------------------------------------------
# data_dir = Path("/home/sirian/Applications/Hydrology/RAVEN/data")

# Do the following snippet for each catchment and model type:
# -----------------------------------------------------------
rootdir = Path("/media/mainman/Work/RAVEN")

for c in catchments_by_id:
    for m in models_by_name:
        model_instance = rt.model.raven_model.RavenModel(model_type=m, catchment_ch_id=c, start_year=start_year,
                                                         end_year=end_year)
        # model_instance.model_dir = Path(rootdir, "models", c, m)
        #        model_instance.create_dirs()
        model_instance.root_dir = rootdir
        model_instance.data_dir = Path(rootdir, "data")
        # model_instance.camels_to_rvt()
        for s in suffix:
            model_instance.write_rvx(rvx_type=s, ostrich_template=True, raven_template=True)
        #        model_instance.write_rvt()
        model_instance.write_ost(max_iterations=50)
        # model_instance.create_symlinks(forcings=False, discharge=True, raven_executable=True, ostrich_executable=True,
        #                                rvx_files=False, raven_diag=False, delete=False)
