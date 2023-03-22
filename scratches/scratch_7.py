import raven_tools as rt

models_by_name = rt.config.variables.supported_models
# catchments_by_id = list(rt.config.variables.catchments.keys())
catchments_by_id = [key for key in rt.config.variables.catchments]
catchments_by_id = ["CH-0058"]
dem_filepaths = []
with open("/media/mainman/Work/RAVEN/data/DEM/dem_tif_filepaths.txt", "r") as f:
    dem_filepaths = f.read().splitlines()
# for c in catchments_by_id:
#     for f in rt.config.variables.forcings_dirs:
#         model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=c, start_year=1981, end_year=2020)
#         model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
#         model_instance.clip_netcdf(forcing_prefix=f)
#         model_instance.create_grid_weights(forcing_name=f)
# print("Finished")
with open("/home/mainman/PycharmProjects/raven-tools/RAVEN/data/glaciers/glaciation_ratio.txt", 'w') as f:
    for c in catchments_by_id:
        model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=c, start_year=1981, end_year=2020)
        glacier_area, glacier_height = model_instance.glaciation_ratio(dem_tif_filenames=dem_filepaths)
        f.write(f"{c} - Area: {glacier_area}\n")
        f.write(f"{c} - Height: {glacier_height}\n")

# model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=catchments_by_id[0], start_year=1981, end_year=2020)

# glaciation_ratio_height = model_instance.glaciation_ratio(dem_tif_filenames=dem_filepaths)
# print(glaciation_ratio_height)
