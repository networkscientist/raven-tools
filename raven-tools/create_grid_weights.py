"""Grid Weight Creator

This script creates a file with the grid weights as required by RAVEN.

It takes the catchment shape file and a netCDF gridded file.

Currently, the dependencies have to be installed through conda-forge (create a new environment for this), at least on
my computer.

Please note that the netCDF coordinates start with (x,y)=(1,1) bottom left.
"""

from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from geopandas import GeoDataFrame
from pandas import DataFrame
from shapely.geometry import Polygon
from xarray import Dataset, DataArray

home_path = Path.home()
output_file = Path(
    "../../../../../media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/GridWeights"
    ".txt")
bounding_box_filename = Path(
    "../../../../../media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/bbox.shp")
extent_shape_file_path = Path("/media/mainman/Data/RAVEN/data/Catchment/Broye_Payerne.shp")  # The shape file that 
# defines the extent of the catchment 
netCDF_file_path = Path(
    "../../../../../media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out"
    "/RhiresD_ch01h.swiss.lv95_196201010000_196212310000_clipped.nc")  # The netCDF file that will be used to get the bounds/extent


def create_overlay(grd: GeoDataFrame, ctm: GeoDataFrame):
    """Overlays a GeoDataFrame over another to create overlay Polygons

    :type ctm: GeoDataFrame
    :type grd: GeoDataFrame
    :rtype: GeoDataFrame
    :param grd:
    :param ctm:
    :return res_u:
    :return res_d:
    """

    res_u = ctm.overlay(grd, how='intersection')
    res_d = grd.overlay(ctm, how='difference')
    res_d = res_d.rename_geometry('geometry')
    res_u.set_index("cell_id")
    res_d.set_index("cell_id")
    return res_u, res_d


def calc_relative_area(gdf: GeoDataFrame):
    """Calculates the relative area of each polygon in a GeoDataFrame.

    Calculates the relative area of each polygon in a GeoDataFrame with EPSG=2056, writes it into a new column and returns the GeoDataFrame with EPSG=2056.
    :param GeoDataFrame gdf: GeoDataFrame in EPSG=2056
    :return GeoDataFrame gdf: GeoDataFrame in EPSG=2056 (LV95)
    """

    # Set the CRS to WGS84 to preserve areas
    gdf = gdf.to_crs(32632)
    # For each intersected feature, compute the area and write into a new field
    gdf["area"] = gdf["geometry"].area
    # Convert the area to float
    area_sum: float = float(gdf["area"].sum())
    for index, row in gdf.iterrows():
        # To minimize numerical errors, multiply and then divide the result by 1000000
        gdf.at[index, "area_rel"] = ((gdf.loc[index]["area"] * 1000000) / area_sum) / 1000000
        gdf.at[index, "index"] = int(index)
        # res_union.at[index, "cell_id"] = cell_id[index]
        # gdf.at[index, "col_id"] = int(col_id[index])
        # gdf.at[index, "row_id"] = int(row_id[index])
    # Re-project back into the LV95 projection and export to a shape file
    gdf = gdf.to_crs(2056)
    return gdf


def write_weights_to_file(grd, filename):
    """

    :type grd: object
    :type filename: object
    :param filename:
    :param grd:
    """
    # Write to GridWeights.txt
    data = write_grid_data(grd)
    with open(filename, 'w') as ff:
        ff.write(':GridWeights                     \n')
        ff.write('   #                                \n')
        ff.write('   # [# HRUs]                       \n')
        ff.write('   :NumberHRUs       {0}            \n'.format(1))
        ff.write('   :NumberGridCells  {0}            \n'.format(len(grid.index)))
        ff.write('   #                                \n')
        ff.write('   # [HRU ID] [Cell #] [w_kl]       \n')
        for idata in data:
            ff.write("   1   {0}   {1}\n".format(idata[1], idata[0]))
        ff.write(':EndGridWeights \n')


def write_grid_data(grd):
    """

    :type grd: object
    :param grd:
    :return:
    """
    # Loop over each intersected feature and write the relative area (compared with the total catchment area) into a new
    # field.
    data_to_write = []
    for index, row in grd.iterrows():
        data_to_write.append(
            [float(grd.loc[index]["area_rel"]), grd.loc[index]["cell_id"]])
    return data_to_write


def copy_rel_area_from_union_to_grid(uni, grd):
    """

    :type grd: object
    :param uni:
    :param grd:
    :return:
    """
    for index, row in uni.iterrows():
        cell_id_old_value = uni.at[index, "cell_id"]
        area_rel_old_value = uni.at[index, "area_rel"]
        ind = grd[grd['cell_id'] == cell_id_old_value].index.tolist()
        grd.at[ind[0], "area_rel"] = area_rel_old_value
    return grd


# Read in the bounding box shape file from the clipping folder as a GeoPandas DataFrame
bbox: GeoDataFrame = gpd.read_file(bounding_box_filename)
bbox.set_crs(21781, allow_override=True)
bbox = bbox.to_crs("EPSG:2056")

# Read in the clipped netCDF file into a xarray Dataset
ds: Dataset = xr.open_dataset(netCDF_file_path, engine="netcdf4")

# Since we're only interested in the grid (not the actual cell values), only use 1 day
ds = ds.sel(time=slice('1962-01-01', '1962-01-02'))
# Select the actual data value column
# TODO: Check Dataset column - Necessary?
#   Since the value are never accessed, do we need to select the column?
xarr: DataArray | Dataset = ds['RhiresD']
# Convert Dataset to DataFrame and reset the index
df: DataFrame = xarr.to_dataframe().reset_index()
# Convert the DataFrame to a GeoDataFrame, using the northing and easting provided by the netCDF. Furthermore,
# change the projection to LV95
combi_catchment_grid: GeoDataFrame = gpd.GeoDataFrame(df.RhiresD, geometry=gpd.points_from_xy(df.E, df.N),
                                                      crs="EPSG:2056")

# From the GeoDataFrame that contains the netCDF grid, get the extent as points
xmin, ymin, xmax, ymax = combi_catchment_grid.total_bounds

# Since the netCDF has a cell size of 1000m, set this here
length: int = 1000
wide: int = 1000

# Create the northing and easting coordinates of the bounding vertex for each cell
cols: list[Any] = list(np.arange(xmin, xmax + wide, wide))
rows: list[Any] = list(np.arange(ymin, ymax + length, length))

# initialize the Polygon list
polygons: list[Polygon] = []
cell_id: list[str] = []

# Create the GeoDataFrame with the grid, providing column names plus polygons for the geometry
grid_cols = ['row', 'col', 'cell_id', 'polygons', 'area', 'area_rel']
grid: GeoDataFrame = gpd.GeoDataFrame(columns=grid_cols, geometry='polygons')
# Loop over each cell and create the corresponding Polygon, then append it to the Polygon list
# Also write the cell id from cell_id=i_row*n_columns + i_column
for ix, x in enumerate(cols):
    for iy, y in enumerate(rows):
        polygons.append(Polygon([(x, y), (x + wide, y), (x + wide, y + length), (x, y + length)]))
        cid = int(iy * len(cols) + ix)
        cell_id.append(f"{str(cid)}")

# Use the polygon list in the GeoDataFrame and set the projection accordingly
grid["polygons"] = polygons
grid["cell_id"] = cell_id

# Every cell that is not within the catchment area will have area_rel set to zero
grid["area_rel"] = 0
grid = grid.set_crs(2056)

# Export the grid to a shape file
grid.to_file("/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/grid.shp")

# Read the catchment shape file into a GeoDataFrame and set the projection accordingly
catchment = gpd.read_file(extent_shape_file_path)
catchment.set_crs(2056)

# Combine the features of the catchment and grid layers into one new GeoDataFrame (so a union can take place) and set
# the projection accordingly.
combi_catchment_grid = gpd.GeoDataFrame(pd.concat([catchment, grid]))
combi_catchment_grid.set_crs(2056)

# Create union and difference overlay GeoDataFrames
res_union, res_diff = create_overlay(grid, catchment)

# Compute the relative area a.k.a grid weight and write to shape files
res_union = calc_relative_area(res_union)
res_union.to_file("/media/mainman/Data/RAVEN/data/union.shp")
res_diff = calc_relative_area(res_diff)
res_diff.to_file("/media/mainman/Data/RAVEN/data/difference.shp")
grid.set_index("cell_id")

# Write a RAVEN compatible grid weights file
copy_rel_area_from_union_to_grid(res_union, grid)
write_weights_to_file(grid, output_file)
