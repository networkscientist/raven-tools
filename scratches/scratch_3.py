import raven_model as rm

# import raven_run as rr
# gr4j_broye = rm.RavenModel(model_type="GR4J", catchment="Broye")
# gr4j_broye = rm.RavenModel(model_type="HYMOD", catchment="Broye")
# gr4j_broye = rm.RavenModel(model_type="GR4J", catchment="Broye")
# gr4j_broye.create_dirs()
suffix = [
    "rvi", "rvh", "rvp", "rvc", "rvt"
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
for m in models:
    model_instance = rm.RavenModel(model_type=m, catchment="Broye")
    model_instance.create_dirs()
    for s in suffix:
        model_instance.write_rvx(ostrich_template=True, rvx_type=s)
gr4j_broye.write_rvx(rvx_type="rvi")
gr4j_broye.write_rvx(rvx_type="rvh")
gr4j_broye.write_rvx(rvx_type="rvp")
gr4j_broye.write_rvx(rvx_type="rvc")
gr4j_broye.write_rvx(rvx_type="rvt")
# gr4j_broye.write_rvp()
# gr4j_broye.write_rvh(ostrich_template=True)
# gr4j_broye.write_rvi(ostrich_template=True)
# gr4j_broye.write_rvc(ostrich_template=True)


# rootdir = Path('/home/mainman/PycharmProjects/raven-tools/RAVEN')
# data_dir = Path('/mnt/data/RAVEN/data')
# generate_ostrich_template(model_type="GR4J", csv_file=(pandas.read_csv(Path(rootdir, data_dir, "Hydromap Attributes", "CH-0057_attributes.csv"), sep=",",
#                                skiprows=[8],
#                                index_col='attribute_names')),)
