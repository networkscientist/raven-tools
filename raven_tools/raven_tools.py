import raven_model as rm

suffix = [
    "rvi", "rvh", "rvp", "rvc", "rvt"
]

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
    model_instance.create_symlinks()
    for s in suffix:
        model_instance.write_rvx(ostrich_template=True, rvx_type=s)
    model_instance.write_ost()

import os
src = '/home/mainman/PycharmProjects/raven-tools/raven_tools/preappend_dblquotes.py'
dst = '/home/mainman/PycharmProjects/raven-tools/raven_tools/preappend_dblquotes_symlink.py'
os.symlink(src, dst)

import raven_tools.raven_preprocess as rpe
from pathlib import Path
import geopandas as gpd
import pandas as pd

home_path = Path.home()
netCDF_file_path = Path(
        "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out"
        "/RhiresD_ch01h.swiss.lv95_196201010000_196212310000_clipped.nc")
bounding_box_filename = Path(
        "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/bbox.shp")
output_file = Path(
    "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/GridWeights.txt")
bounding_box_filename = Path(
    "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/bbox.shp")
# The shape file that defines the extent of the catchment
extent_shape_file_path = Path("/media/mainman/Data/RAVEN/data/Catchment/Broye_Payerne.shp")
# The netCDF file that will be used to get the bounds/extent
netCDF_file_path = Path(
    "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out"
    "/RhiresD_ch01h.swiss.lv95_196201010000_196212310000_clipped.nc")

grid = rpe.create_grid(netCDF_file_path, bounding_box_filename, export_shp=False)
# Read the catchment shape file into a GeoDataFrame and set the projection accordingly
catchment = gpd.read_file(extent_shape_file_path)
catchment.set_crs(2056)

# Combine the features of the catchment and grid layers into one new GeoDataFrame (so a union can take place) and set
# the projection accordingly.
combi_catchment_grid = gpd.GeoDataFrame(pd.concat([catchment, grid]))
combi_catchment_grid.set_crs(2056)

# Create union and difference overlay GeoDataFrames
res_union, res_diff = rpe.create_overlay(grid, catchment)

# Compute the relative area a.k.a grid weight and write to shape files
res_union = rpe.calc_relative_area(res_union)
res_union.to_file("/media/mainman/Data/RAVEN/data/union.shp")
res_diff = rpe.calc_relative_area(res_diff)
res_diff.to_file("/media/mainman/Data/RAVEN/data/difference.shp")
grid.set_index("cell_id")

# Write a RAVEN compatible grid weights file
grid: gpd.GeoDataFrame = rpe.copy_rel_area_from_union_to_grid(res_union, grid)
rpe.write_weights_to_file(grid, output_file)