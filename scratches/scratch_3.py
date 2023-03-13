# import raven_model as rm
# import model.raven_model
from pathlib import Path

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
models = [
    "GR4J",
    "HYMOD",
    "HMETS",
    "HBV",
    "MOHYSE"
]

# catchments = ["Ticino",
#               "Broye",
#               "Thur",
#               "Massa",
#               "Weisse_Luetschine",
#               "Dischmabach"]
c = "Dischmabach"
m = "GR4J"
model_instance = rt.model.raven_model.RavenModel(model_type=m, catchment=c)
model_instance.start_year = 1981
model_instance.end_year = 2020
model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
model_instance.bbox_filepath = Path(model_instance.data_dir, "Catchment",
                                    f"{model_instance.catchment}_bbox.shp")
forcing_dir = rt.config.variables.forcings_dirs[0]
model_instance.create_netcdf(forcing_dir=forcing_dir, clip=True, merge=False)
# for n in rt.config.variables.forcings_dirs:
#     model_instance.create_netcdf(forcing_dir=n, clip=True, merge=True)
# for c in catchments:
#     model_instance.catchment = c
#     model_instance.create_symlinks(ostrich_executable=False, forcings=False, discharge=False)
# model_instance.create_grid_weights(forcing_name="RhiresD_v2.0_swiss.lv95")

# for c in catchments:
#     for m in rt.config.variables.supported_models:
#         model_instance = rt.model.raven_model.RavenModel(model_type=m, catchment=c)
#         #         #         # print(model_instance.start_year)
#         #         #         # model_instance.create_dirs()
#         model_instance.start_year = 1981
#         model_instance.end_year = 2020
#         #         #         # model_instance.camels_to_rvt()
#         model_instance.create_symlinks()
#         # model_instance.bbox_filepath = Path(model_instance.data_dir, "Catchment",
#         #                                             f"{model_instance.catchment}_bbox.shp")
#         model_instance.data_dir = Path("/media/mainman/Work/RAVEN/data")
#         model_instance.write_rvt(ostrich_template=True, raven_template=True)
#                         for n in config.variables.forcings_dirs:
#                             model_instance.create_netcdf(forcing_dir=n, clip=True, merge=True)
#
#                 for n in config.variables.forcings_dirs:
#                     model_instance.create_grid_weights(forcing_name=n)
#         for s in suffix:
#             model_instance.write_rvx(ostrich_template=True, rvx_type=s)
#             model_instance.write_ost()
# #         # model_instance.create_netcdf(clip=False, merge=True)
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
