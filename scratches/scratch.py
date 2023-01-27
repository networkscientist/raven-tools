import raven_run as rr

# with open("raven_tools/project_config.yaml", "r") as f:
#     project_config = yaml.load(f, Loader=yaml.FullLoader)
# model_dir = Path(project_config['ModelDir'])
# model_name = str(project_config['ModelName'])
rr.write_rvt(1981, 2010)
