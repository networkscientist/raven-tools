# import raven_model as rm
# import model.raven_model
import config.variables
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

catchments = ["Ticino",
              "Broye",
              "Thur",
              "Massa",
              "Weisse Lütschine",
              "Dischmabach"]

for c in catchments:
    for m in rt.config.variables.supported_models:
        model_instance = rt.model.raven_model.RavenModel(model_type=m, catchment=c)
        print(model_instance.start_year)
        model_instance.create_dirs()
        model_instance.start_year = 1980
        model_instance.end_year = 1989
        model_instance.write_rvt()
        model_instance.camels_to_rvt()
        model_instance.create_symlinks()
        for n in config.variables.forcings_dirs:
            model_instance.create_netcdf(forcing_dir=n)
        for s in suffix:
            model_instance.write_rvx(ostrich_template=True, rvx_type=s)
            model_instance.write_ost()
# gr4j_broye.write_rvx(rvx_type="rvi")
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
