"""Grid Weight Creator

This script creates a file with the grid weights as required by RAVEN.

It takes the catchment shape file and a netCDF gridded file.

Currently, the dependencies have to be installed through conda-forge (create a new environment for this), at least on
my computer.
"""

from pathlib import Path
from typing import Any, List

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from geopandas import GeoDataFrame
from pandas import DataFrame
from shapely.geometry import Polygon
from xarray import Dataset, DataArray

home_path = Path.home()
bounding_box_filename = Path(
    "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/bbox.shp")
extent_shape_file_path = Path("/media/mainman/Data/RAVEN/data/Catchment/Broye_Payerne.shp")  # The shape file that 
# defines the extent of the catchment 
netCDF_file_path = Path(
    "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/RhiresD_ch01h.swiss.lv95_196201010000_196212310000_clipped.nc")  # The netCDF file that will
# be used to get the bounds/extent 

# Read in the bounding box shape file from the clipping folder as a GeoPandas DataFrame
bbox: DataFrame | GeoDataFrame | Any = gpd.read_file(bounding_box_filename)
bbox.set_crs(21781, allow_override=True)
bbox = bbox.to_crs("EPSG:2056")

# Read in the clipped netCDF file into a xarray Dataset
ds: Dataset = xr.open_dataset(netCDF_file_path, engine="netcdf4")
# Since we're only interested in the grid (not the actual cell values), only use 1 day
# TODO: Check date range
#   Check that date range is actually for one single day (excluding the end date).
ds = ds.sel(time=slice('1962-01-01', '1962-01-02'))
# Select the actual data value column
# TODO: Check Dataset column - Necessary?
#   Since the value are never accessed, do we need to select the column?
xarr: DataArray | Dataset = ds['RhiresD']
# Convert Dataset to DataFrame and reset the index
df: DataFrame | None = xarr.to_dataframe().reset_index()
# Convert the DataFrame to a GeoDataFrame, using the northing and easting provided by the netCDF. Furthermore,
# change the projection to LV95
gdf: GeoDataFrame = gpd.GeoDataFrame(df.RhiresD, geometry=gpd.points_from_xy(df.E, df.N), crs="EPSG:2056")

# TODO: Write a function that creates a shape file with the netCDF's grid
# TODO: Write a function that takes the Union of the grid shape file and the catchment shape file
# TODO: Write a function that calculates the area of each grid cell that is overlapping with the catchment area and return the ratio of the total catchment area
# TODO: Write a function that exports the grid weights into a RAVEN compatible file

# From the GeoDataFrame that contains the netCDF grid, get the extent as points
xmin, ymin, xmax, ymax = gdf.total_bounds

# Since the netCDF has a cell size of 1000m, set this here
length: int = 1000
wide: int = 1000

# Create the northing and easting coordinates of the bounding vertex for each cell
cols: list[Any] = list(np.arange(xmin, xmax + wide, wide))
rows: list[Any] = list(np.arange(ymin, ymax + length, length))

# initialize the Polygon list
polygons: list[Polygon] = []

# Loop over each cell and create the corresponding Polygon, then append it to the Polygon list
for x in cols[:-1]:
    for y in rows[:-1]:
        polygons.append(Polygon([(x, y), (x + wide, y), (x + wide, y + length), (x, y + length)]))

# Create a GeoDataFrame from the Polygon list and set the projection accordingly
grid: GeoDataFrame = gpd.GeoDataFrame({'geometry': polygons})
grid = grid.set_crs(2056)

# Export the grid to a shape file
grid.to_file("/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/grid.shp")

# Read the catchment shape file into a GeoDataFrame and set the projection accordingly
catchment = gpd.read_file(extent_shape_file_path)
catchment.set_crs(2056)

# Combine the features of the catchment and grid layers into one new GeoDataFrame (so a union can take place) and set
# the projection accordingly.
gdf = gpd.GeoDataFrame(pd.concat([catchment, grid]))
gdf.set_crs(2056)

# Overlay the catchment with the grid and keep the features which intersect
res_union = catchment.overlay(grid, how='intersection')
# Set the CRS to WGS84 to preserve areas
res_union = res_union.to_crs(32632)

# For each intersected feature, compute the area and write into a new field
res_union["area"] = res_union['geometry'].area
# Convert the area to float
area_sum = float(res_union["area"].sum())

# Loop over each intersected feature and write the relative area (compared with the total catchment area) into a new
# field.
for index, row in res_union.iterrows():
    # To minimize numerical errors, multiply and then divide the result by 1000000
    res_union.at[index, "area_rel"] = ((res_union.loc[index]["area"] * 1000000) / area_sum) / 1000000
    res_union.at[index, "index"] = index

# Re-project back into the LV95 projection and export to a shape file
res_union = res_union.to_crs(2056)
res_union.to_file("/media/mainman/Data/RAVEN/data/union.shp")
