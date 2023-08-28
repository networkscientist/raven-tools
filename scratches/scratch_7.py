from pathlib import Path

import pandas as pd

import raven_tools as rt

# models_by_name = rt.config.variables.supported_models
models_by_name = ["HBV"]
# catchments_by_id = list(rt.config.variables.catchments.keys())
catchments_by_id = [key for key in rt.config.variables.catchments]
# catchments_by_id = ["CH-0140"]
# catchments_by_id = ["CH-0058"]
# dem_filepaths = []
# with open(f"/media/mainman/Work/RAVEN/data/DEM/linklist_{catchments_by_id[0]}.csv", "r") as f:
#     dem_filepaths = f.read().splitlines()
data_dir = Path("/media/mainman/Work/RAVEN/data")

hru_info_df_orig = pd.read_csv(Path(data_dir, "Catchment", "hru_info.csv"), na_values='-')

for c in catchments_by_id:
    # c = "CH-0053"
    hru_info_df = hru_info_df_orig[hru_info_df_orig['Ctm'] == c]
    hru_info_dict = hru_info_df.to_dict('records')[0]
    # for f in rt.config.variables.forcings_dirs:
    model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=c, start_year=1981, end_year=2020,
                                                     model_type=models_by_name[0])
    # model_instance.bbox_filepath = Path(model_instance.data_dir, "Catchment", f"{model_instance.catchment}_bbox.shp")
    model_instance.data_dir = data_dir
    # model_instance.create_bbox()
    # # model_instance.clip_netcdf(forcing_prefix=f)
    # model_instance.bbox_filepath = Path(model_instance.data_dir, "Catchment",
    #                                     f"{model_instance.catchment_ch_id}_bbox.shp")

    # if not math.isnan(hru_info_dict['GlaArea']):
    #     model_instance.create_grid_weights(forcing_name=f, glacier=True)
    # else:
    #     model_instance.create_grid_weights(forcing_name=f, glacier=False)
    model_instance.create_grid_weights(read_from_file=False)
print("Finished")

# This snippet gets the glaciation ratio for the non glaciated and for the glaciated HRU and writes the values with
# corresponding centroid coordinates into a file
# -----------------------------------------------------------------------------------------------------------------

# for c in catchments_by_id:
#     with open(f"/media/mainman/Work/RAVEN/data/DEM/linklist_{c}.csv", "r") as f:
#         dem_filepaths = f.read().splitlines()
#
#     model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=c, start_year=1981, end_year=2020)
#     model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
#     for i, s in enumerate(dem_filepaths):
#         dem_filepaths[i] = Path(model_instance.data_dir, "DEM", "out", os.path.basename(s))
#     try:
#         glaciation_ratio, glacier_height, non_glaciation_ratio, non_gla_height, glaciated_centroid, non_glaciated_centroid = model_instance.glacier_ratio(
#             dem_tif_filenames=dem_filepaths)
#         with open(
#                 f"/home/mainman/PycharmProjects/raven-tools/RAVEN/data/glaciers/glaciation_ratio_{c}.txt",
#                 'w') as f:
#             f.write(
#                 f"Gla_Ratio;Gla_Alti;Non_Gla_Ratio;Non_Gla_Alti;Gla_Cent_Lat;Gla_Cent_Lon;Non_Gla_Cent_Lat;Non_Gla_Cent_Lon\n")
#             f.write(
#                 f"{glaciation_ratio};{glacier_height};{non_glaciation_ratio};{non_gla_height};{glaciated_centroid.x};{glaciated_centroid.y};{non_glaciated_centroid.x};{non_glaciated_centroid.y}\n")
#     except ValueError:
#         print("There has been an error clipping DEM")
#         pass

# glacier_props = pd.DataFrame(columns=[
#     "NonGlaArea",
#     "GlaArea",
#     "NonGlaAlti",
#     "GlaAlti",
#     "NonGlaLat",
#     "GlaLat",
#     "NonGlaLon",
#     "GlaLon",
#     "NonGlaAspect",
#     "GlaAspect",
#     "NonGlaSlope",
#     "GlaSlope"
# ])
# glacier_props.rename_axis("Ctm", inplace=True)
# for c in catchments_by_id:
#     model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=c, start_year=1981, end_year=2020)
#     model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
#     glacier_info = pd.read_csv(
#         f"/home/mainman/PycharmProjects/raven-tools/RAVEN/data/glaciers/glaciation_ratio_{c}.txt", sep=';')
#     # model_instance.hru_aspect_slope(write_to_file=True)
#
#     if float(glacier_info['Gla_Ratio'][0]) != 0.0:
#         prop_df: pd.DataFrame = model_instance.hru_dem_props(glacier=True)
#
#     else:
#         prop_df: pd.DataFrame = model_instance.hru_dem_props(glacier=False)
#     glacier_props = pd.concat([glacier_props, prop_df])
# glacier_props.to_csv(Path("/media/mainman/Work/RAVEN/data", "Catchment", "hru_info.csv"), na_rep='-')

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
