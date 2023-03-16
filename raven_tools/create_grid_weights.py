"""
Grid Weights Creator

This script creates a file with the grid weights as required by RAVEN. It takes the catchment shape file and a netCDF gridded file.
Currently, the dependencies have to be installed through conda-forge (create a new environment for this), at least on
my computer.
Please note that the netCDF coordinates start with (x,y)=(1,1) bottom left.
"""
# TODO: Check wether deps still have to be installed through conda-forge...

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from geopandas import GeoDataFrame
from pandas import DataFrame
from shapely.geometry import Polygon
from xarray import Dataset


def create_overlay(grd: GeoDataFrame, ctm: GeoDataFrame):
    """Overlays a GeoDataFrame over another to create overlay Polygons

    Overlays two GeoDataFrame over each other and returns two new GeoDataFrames, one for the mode 'intersection',
    the second for the mode 'difference'

    Parameters
    ----------
    grd : GeoDataFrame
        Grid as given by the netCDF file
    ctm : GeoDataFrame
        Catchment as given by a shape file

    Returns
    -------
    res_u : GeoDataFrame
        Grid cells within the catchment area
    res_d : GeoDataFrame
        Grid cells outside the catchment area

    """

    res_u: GeoDataFrame = ctm.overlay(grd, how='intersection')
    print(f"res_u.crs = {res_u.crs}")
    res_d: GeoDataFrame = grd.overlay(ctm, how='difference')
    print(f"res_d.crs = {res_d.crs}")

    res_d = res_d.rename_geometry('geometry')
    res_u.set_index("cell_id")
    res_d.set_index("cell_id")
    return res_u, res_d


def create_overlay(grd: GeoDataFrame, ctm: GeoDataFrame):
    """Overlays a GeoDataFrame over another to create overlay Polygons

    Overlays two GeoDataFrame over each other and returns two new GeoDataFrames, one for the mode 'intersection',
    the second for the mode 'difference'

    Parameters
    ----------
    grd : GeoDataFrame
        Grid as given by the netCDF file
    ctm : GeoDataFrame
        Catchment as given by a shape file

    Returns
    -------
    res_u : GeoDataFrame
        Grid cells within the catchment area
    res_d : GeoDataFrame
        Grid cells outside the catchment area

    """

    res_u: GeoDataFrame = ctm.overlay(grd, how='intersection')
    print(f"res_u.crs = {res_u.crs}")

    res_u.set_index("cell_id")

    return res_u


def calc_relative_area(gdf2056: GeoDataFrame):
    """Calculates the relative area of each polygon in a GeoDataFrame.

    Calculates the relative area of each polygon in a GeoDataFrame with EPSG=2056, writes it into a new column and
    returns the GeoDataFrame with EPSG=2056.

    Parameters
    ----------
    gdf : GeoDataFrame
        GeoDataFrame in EPSG=2056

    Returns
    -------
    gdf : GeoDataFrame
        GeoDataFrame with relative areas of each polygon

    """

    # Set the CRS to WGS84 to preserve areas
    gdf2056.set_crs("EPSG:2056", allow_override=True)
    gdf32632 = gdf2056.to_crs(epsg='32632')
    # For each intersected feature, compute the area and write into a new field
    gdf32632["area"] = gdf32632["geometry"].area
    # Convert the area to float
    area_sum: float = float(gdf32632["area"].sum())
    for index, row in gdf32632.iterrows():
        # To minimize numerical errors, multiply and then divide the result by 1000000
        gdf32632.at[index, "area_rel"] = ((gdf32632.loc[index]["area"] * 1000000) / area_sum) / 1000000
        gdf32632.at[index, "index"] = int(index)
    # Re-project back into the LV95 projection and export to a shape file
    gdf2056: GeoDataFrame = gdf32632.to_crs(2056)
    return gdf2056


def write_weights_to_file(grd: GeoDataFrame, filename: Path):
    """Write grid weights to Raven compatible file

    Parameters
    ----------
    grd : GeoDataFrame
        Grid derived from the netCDF file
    filename : Path
        Path to the grid weight text file

    """
    # Write to GridWeights.txt
    data = write_grid_data(grd)
    # The following section has been adapted from Juliane Mai's derive_grid_weights.py script. It writes to a file
    # that is compatible with RAVEN
    with open(filename, 'w') as ff:
        ff.write(':GridWeights                     \n')
        ff.write('   #                                \n')
        ff.write('   # [# HRUs]                       \n')
        ff.write('   :NumberHRUs       {0}            \n'.format(1))  # Currently for GR4J, there's 1 HRU
        ff.write('   :NumberGridCells  {0}            \n'.format(len(grid.index)))
        ff.write('   #                                \n')
        ff.write('   # [HRU ID] [Cell #] [w_kl]       \n')
        for idata in data:
            ff.write("   1   {0}   {1}\n".format(idata[1], idata[0]))
        ff.write(':EndGridWeights \n')


def write_grid_data(grd: GeoDataFrame):
    """Loops over each grid cell and extracts the grid weights.

    Parameters
    ----------
    grd : GeoDataFrame
        Grid as derived from the netCDF file

    Returns
    -------
    data_to_write : list
        List with the relative areas/grid weights of each cell.

    """
    # Loop over each intersected feature and write the relative area (compared with the total catchment area) into a new
    # field.
    data_to_write: list[list[float]] = []
    for index, row in grd.iterrows():
        data_to_write.append(
            [float(grd.loc[index]["area_rel"]), grd.loc[index]["cell_id"]])
    return data_to_write


def copy_rel_area_from_union_to_grid(uni: GeoDataFrame, grd: GeoDataFrame):
    """Takes grid weights from a union GeoDataFrame and writes the to the grid GeoDataFrame.

    Parameters
    ----------
    uni : GeoDataFrame
        GeoDataFrame containing the grid cells within the catchment.
    grd : GeoDataFrame
        Grid GeoDataFrame as derived from netCDF file

    Returns
    -------
    grd : GeoDataFrame
        Grid GeoDataFrame with grid weights

    """

    # Loop over the union GeoDataFrame, take relative area/grid weight and write it to corresponding cell in grid
    # GeoDataFrame.
    for index, row in uni.iterrows():
        cell_id_old_value = uni.at[index, "cell_id"]
        area_rel_old_value = uni.at[index, "area_rel"]
        ind = grd[grd['cell_id'] == cell_id_old_value].index.tolist()
        grd.at[ind[0], "area_rel"] = area_rel_old_value
    return grd


if __name__ == '__main__':
    home_path = Path.home()
    output_file = Path(
        "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/tmp/script_GridWeights.txt")
    bounding_box_filename = Path(
        "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/tmp/bbox.shp")
    # The shape file that defines the extent of the catchment
    extent_shape_file_path = Path(
        "/media/mainman/Work/RAVEN/data/Catchment/reproject_2056/CH-0159.shp")
    # The netCDF file that will be used to get the bounds/extent
    netCDF_file_path = Path(
        "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/tmp/CH-0159_clipped.nc")

    # Read in the bounding box shape file from the clipping folder as a GeoPandas DataFrame
    # bbox: GeoDataFrame = gpd.read_file(bounding_box_filename)
    # bbox.set_crs(2056, allow_override=True)
    # bbox = bbox.to_crs("EPSG:2056")

    # Read in the clipped netCDF file into a xarray Dataset
    ds: Dataset = xr.open_dataset(netCDF_file_path, engine="netcdf4")

    # Since we're only interested in the grid (not the actual cell values), only use 1 day
    ds = ds.sel(time=slice('2000-01-01', '2000-01-02'))
    # Select the actual data value column
    xarr = ds['RhiresD']
    # Convert Dataset to DataFrame and reset the index
    df: DataFrame = xarr.to_dataframe().reset_index()
    # Convert the DataFrame to a GeoDataFrame, using the northing and easting provided by the netCDF. Furthermore,
    # change the projection to LV95
    combi_catchment_grid: GeoDataFrame = gpd.GeoDataFrame(df.RhiresD,
                                                          geometry=gpd.points_from_xy(df.E, df.N, crs='2056'),
                                                          crs="EPSG:2056")

    # From the GeoDataFrame that contains the netCDF grid, get the extent as points
    xmin, ymin, xmax, ymax = combi_catchment_grid.total_bounds

    # Since the netCDF has a cell size of 1000m, set this here
    length: int = 1000
    wide: int = 1000

    # Create the northing and easting coordinates of the bounding vertex for each cell
    cols: list = list(np.arange(xmin - 500, xmax + 500, wide))
    rows: list = list(np.arange(ymin - 500, ymax + 500, length))

    # initialize the Polygon list
    polygons: list[Polygon] = []
    cell_id: list[str] = []

    # Create the GeoDataFrame with the grid, providing column names plus polygons for the geometry
    grid_cols = ['row', 'col', 'cell_id', 'polygons', 'area', 'area_rel']
    grid: GeoDataFrame = gpd.GeoDataFrame(columns=grid_cols, geometry='polygons')
    grid.set_crs(epsg='2056')
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
    # grid = grid.set_crs(2056)

    # Export the grid to a shape file
    grid.to_file(
        "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/tmp/script_grid.shp")
    # _-----------------------

    # Read the catchment shape file into a GeoDataFrame and set the projection accordingly
    catchment = gpd.read_file(extent_shape_file_path)
    catchment.set_crs("EPSG:2056", allow_override=True)

    # Combine the features of the catchment and grid layers into one new GeoDataFrame (so a union can take place) and set
    # the projection accordingly.
    combi_catchment_grid = gpd.GeoDataFrame(pd.concat([catchment, grid]))
    combi_catchment_grid.set_crs("EPSG:2056", allow_override=True)

    # Create union and difference overlay GeoDataFrames
    catchment.set_crs(epsg='2056')
    res_union = create_overlay(grid, catchment)
    res_union.set_crs(2056)
    # res_diff.set_crs(2056)

    # Compute the relative area a.k.a grid weight and write to shape files
    res_union = calc_relative_area(res_union)
    res_union.set_crs(2056)
    # res_union.to_file(
    #     "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/tmp/script_union.shp")
    # res_diff.set_crs(2056, allow_override=True)
    # res_diff = calc_relative_area(res_diff)
    # res_diff.to_file(
    #     "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/tmp/script_difference.shp")
    grid.set_index("cell_id")

    # Write a RAVEN compatible grid weights file
    copy_rel_area_from_union_to_grid(res_union, grid)
    write_weights_to_file(grid, output_file)
