import raven_run as rr
from pathlib import Path
import yaml

# with open("raven_tools/config.yaml", "r") as f:
#     config = yaml.load(f, Loader=yaml.FullLoader)
# model_dir = Path(config['ModelDir'])
# model_name = str(config['ModelName'])
rr.write_rvt(1981, 2010)
