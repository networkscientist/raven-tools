# import raven_model as rm
# import model.raven_model

import raven_tools as rt

# import raven_tools.model as rm
# import raven_run as rr
# gr4j_broye = rm.RavenModel(model_type="GR4J", catchment="Broye")
# gr4j_broye = rm.RavenModel(model_type="HYMOD", catchment="Broye")
# gr4j_broye = rm.RavenModel(model_type="GR4J", catchment="Broye")
# gr4j_broye.create_dirs()
suffix = [
    "rvi", "rvh", "rvp", "rvc"
]
# gr4j_broye.write_rvh()
# gr4j_broye.write_rvx()
models_by_name = rt.config.variables.supported_models
catchments_by_id = [key for key in rt.config.variables.catchments]
start_year = 1981
end_year = 2020
# m = "GR4J"
# model_instance = rt.model.raven_model.RavenModel(model_type=m, catchment=catchments[0])
# model_instance.start_year = 1981
# model_instance.end_year = 2020
# model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
# # model_instance.bbox_filepath = Path(model_instance.data_dir, "Catchment",
# #                                     f"{model_instance.catchment}_bbox.shp")
# forcing_dir = rt.config.variables.forcings_dirs[0]
# model_instance.create_netcdf(forcing_dir=forcing_dir, clip=True, merge=False)
# for n in rt.config.variables.forcings_dirs:
#     model_instance.create_netcdf(forcing_dir=n, clip=True, merge=True)
# for c in catchments:
#     model_instance.catchment = c
#     model_instance.create_symlinks(ostrich_executable=False, forcings=False, discharge=False)
# model_instance.create_grid_weights(forcing_name="RhiresD_v2.0_swiss.lv95")

# Do the following snippet for each forcing_dir
# ---------------------------------------------

# for c in catchments_by_id:
#     for f in rt.config.variables.forcings_dirs:
#         model_instance = rt.model.raven_model.RavenModel(catchment_ch_id=c, start_year=1981, end_year=2020)
#         model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
#         #         #         #         # model_instance.camels_to_rvt()
#         #         model_instance.create_symlinks()
#         model_instance.bbox_filepath = Path(model_instance.data_dir, "Catchment",
#                                             f"{model_instance.catchment}_bbox.shp")
# model_instance.create_grid_weights(forcing_name=f)
# model_instance.write_rvt(ostrich_template=True, raven_template=True)

# Do the following snippet for each catchment and model type:
# -----------------------------------------------------------

for c in catchments_by_id:
    for m in models_by_name:
        model_instance = rt.model.raven_model.RavenModel(model_type=m, catchment_ch_id=c, start_year=start_year,
                                                         end_year=end_year)
        # model_instance.create_dirs()
        # model_instance.camels_to_rvt()
        for s in suffix:
            model_instance.write_rvx(rvx_type=s, ostrich_template=True, raven_template=True)
        model_instance.write_rvt()
        model_instance.write_ost()
        # model_instance.create_symlinks(rvx_files=False)

# for f in rt.config.variables.forcings_dirs:
# model_instance.merge_netcdf(forcing_prefix=n)
# model_instance.clip_netcdf(forcing_prefix=f)
# model_instance.create_grid_weights(forcing_name=f)
#         for s in suffix:
#             model_instance.write_rvx(ostrich_template=True, rvx_type=s)
#             model_instance.write_ost()
# #         # model_instance.bbox_filepath = Path(model_instance.data_dir, "Catchment", f"{model_instance.catchment}_bbox.shp")
# #         # print(model_instance.bbox_filepath)
# # gr4j_broye.write_rvx(rvx_type="rvi")
# gr4j_broye.write_rvx(rvx_type="rvh")
# gr4j_broye.write_rvx(rvx_type="rvp")
# gr4j_broye.write_rvx(rvx_type="rvc")
# gr4j_broye.write_rvx(rvx_type="rvt")
# gr4j_broye.write_rvp()
# gr4j_broye.write_rvh(ostrich_template=True)
# gr4j_broye.write_rvi(ostrich_template=True)
# gr4j_broye.write_rvc(ostrich_template=True)


# rootdir = Path('/home/mainman/PycharmProjects/raven-tools/RAVEN')
# data_dir = Path('/mnt/data/RAVEN/data')
# generate_ostrich_template(model_type="GR4J", csv_file=(pandas.read_csv(Path(rootdir, data_dir, "Hydromap Attributes", "CH-0057_attributes.csv"), sep=",",
#                                skiprows=[8],
#                                index_col='attribute_names')),)
