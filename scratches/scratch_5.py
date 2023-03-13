# import datetime
from pathlib import Path

import raven_tools.processing.raven_preprocess as rpe

bbox_filename: Path = "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/tmp/bbox.shp"
rpe.create_bbox(Path("/media/mainman/Work/RAVEN/data/Catchment/CH-0159.shp"), bbox_filename, create_bbox_shp=True)
grid_path: Path = "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/tmp/CH-0159_grid.shp"
grid = rpe.create_grid(netcdf_filepath=Path(
    "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/tmp/CH-0159_clipped.nc"),
    bounding_box_filename=bbox_filename, out_path=grid_path, forcing_name="RhiresD", start_year=1981,
    export_shp=True)
