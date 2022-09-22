from typing import Any

import geopandas as gpd
from geopandas import GeoDataFrame
from pandas import DataFrame
from shapely.geometry import Polygon
import numpy as np
from pathlib import Path
import xarray as xr
import cfgrib
from xarray import Dataset

home_path = Path.home()
bounding_box_filename = Path("/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out"
                             "/bbox.shp")
extent_shape_file_path = Path("/media/mainman/Data/RAVEN/data/Catchment/Broye_Payerne.shp")  # The shape file that 
# defines the extent of the catchment 
netCDF_file_path = Path("/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out"
                        "/RhiresD_ch01h.swiss.lv95_196201010000_196212310000_clipped.nc")  # The netCDF file that will 
# be used to get the bounds/extent 

# Read in the bounding box shape file from the clipping folder as a GeoPandas DataFrame
catchment: DataFrame | GeoDataFrame | Any = gpd.read_file(bounding_box_filename)
# Read in the clipped netCDF file into a xarray Dataset
ds: Dataset = xr.open_dataset(netCDF_file_path, engine="netcdf4")
# Convert the xarray Dataset into a Pandas Dataframe
netcdf_df: DataFrame = ds.to_dataframe()

# TODO: Write a function that creates a shape file with the netCDF's grid
# TODO: Write a function that takes the Union of the grid shape file and the catchment shape file
# TODO: Write a function that calculates the area of each grid cell that is overlapping with the catchment area and return the ratio of the total catchment area
# TODO: Write a function that exports the grid weights into a RAVEN compatible file

# xmin, ymin, xmax, ymax = catchment.total_bounds
#
# length = 1000
# wide = 1000
#
# cols = list(np.arange(xmin, xmax + wide, wide))
# rows = list(np.arange(ymin, ymax + length, length))
#
# polygons = []
# for x in cols[:-1]:
#     for y in rows[:-1]:
#         polygons.append(Polygon([(x,y), (x+wide, y), (x+wide, y+length), (x, y+length)]))
#
# grid = gpd.GeoDataFrame({'geometry':polygons})
# grid.to_file("/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/grid.shp")
