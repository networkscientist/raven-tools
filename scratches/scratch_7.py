import os
from pathlib import Path

import raven_tools as rt

models_by_name = rt.config.variables.supported_models
# catchments_by_id = list(rt.config.variables.catchments.keys())
catchments_by_id = [key for key in rt.config.variables.catchments]
# catchments_by_id = ["CH-0057"]
# catchments_by_id = ["CH-0058"]
dem_filepaths = []
# with open(f"/media/mainman/Work/RAVEN/data/DEM/linklist_{catchments_by_id[0]}.csv", "r") as f:
#     dem_filepaths = f.read().splitlines()
# for c in catchments_by_id:
#     for f in rt.config.variables.forcings_dirs:
#         model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=c, start_year=1981, end_year=2020)
#         model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
#         model_instance.clip_netcdf(forcing_prefix=f)
#         model_instance.create_grid_weights(forcing_name=f)
# print("Finished")
for c in catchments_by_id:
    with open(f"/media/mainman/Work/RAVEN/data/DEM/linklist_{c}.csv", "r") as f:
        dem_filepaths = f.read().splitlines()

    model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=c, start_year=1981, end_year=2020)
    model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
    for i, s in enumerate(dem_filepaths):
        dem_filepaths[i] = Path(model_instance.data_dir, "DEM", "out", os.path.basename(s))
    try:
        glaciation_ratio, glacier_height, non_glaciation_ratio, non_gla_height, glaciated_centroid, non_glaciated_centroid = model_instance.glacier_ratio(
            dem_tif_filenames=dem_filepaths)
        with open(
                f"/home/mainman/PycharmProjects/raven-tools/RAVEN/data/glaciers/glaciation_ratio_{c}.txt",
                'w') as f:
            f.write(f"Glac_Ratio;Alti;Non_Gla_Ratio;Non_Gla_Alti;Gla_Cent_X;Gla_Cent_Y;Non_Gla_Cent_X;Non_Gla_Cent_Y\n")
            f.write(
                f"{glaciation_ratio};{glacier_height};{non_glaciation_ratio};{non_gla_height};{glaciated_centroid.x};{glaciated_centroid.y};{non_glaciated_centroid.x};{non_glaciated_centroid.y}\n")
    except ValueError:
        print("There has been an error clipping DEM")
        pass

# with open(f"/home/mainman/PycharmProjects/raven-tools/RAVEN/data/glaciers/glaciation_ratio_{catchments_by_id[0]}.txt",
#           'w') as f:
#     for c in catchments_by_id:
#         model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=c, start_year=1981, end_year=2020)
#         model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
#         for i, s in enumerate(dem_filepaths):
#             dem_filepaths[i] = Path(model_instance.data_dir, "DEM", "out", os.path.basename(s))
#         glacier_area, glacier_height = model_instance.glaciation_ratio(dem_tif_filenames=dem_filepaths)
#         f.write(f"{c} - Area: {glacier_area}\n")
#         f.write(f"{c} - Height: {glacier_height}\n")

# model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=catchments_by_id[0], start_year=1981, end_year=2020)

# glaciation_ratio_height = model_instance.glaciation_ratio(dem_tif_filenames=dem_filepaths)
# print(glaciation_ratio_height)
