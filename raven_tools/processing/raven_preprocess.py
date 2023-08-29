"""
Tools to pre-process files for Raven models.

Functions:
    create_bbox_geometry(Path)
    create_bounding_shape(Path)
    netcdf_to_dataset(Path)
    dataset_to_netcdf(Dataset)
    netcdf_clipper(Path,Path,GeoDataFrame)
    netcdf_clipper_multi(Path,GeoDataFrame)
"""
import glob
import json
import logging
import math
import os
import re
import shutil
import subprocess
from datetime import datetime
from datetime import timedelta
from pathlib import Path
from typing import Union, Any

import dask.dataframe as dd
import dask_geopandas as dgpd
import geopandas as gpd
import netCDF4
import numpy as np
import pandas as pd
import rasterio.errors
import rasterio.mask
import rasterio.merge
import rioxarray as rxr
import xarray as xr
from geopandas import GeoDataFrame
from netCDF4 import Dataset
from numpy import exp
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union

from raven_tools import config

logger = logging.getLogger(__name__)
logger.debug("Entered raven_preprocess.py.")


def create_bbox_gdf(extent_shape_path: Path):
    """Create a bounding box Polygon for an input shape file.

    Args:
        extent_shape_path : Path
            Path to the shape file for which to create a bounding box.

    Returns:
        bbox_poly : Polygon
            A shapely polygon that defines the bounding box.
        ext_gdf : GeoDataFrame
            The GeoDataFrame with the data from the extent shape file. The extent has to be in EPSG:2056

    """
    logger.debug("Entered create_bbox_geometry function...")
    ext_gdf = gpd.read_file(extent_shape_path)
    ext_gdf.set_crs("EPSG:2056")
    # ext_gdf.to_crs("EPSG:2056", inplace=True)
    se = ext_gdf.geometry.total_bounds  # Get bounds and store in array
    bbox_poly = Polygon([[se[0], se[1]], [se[2], se[1]], [se[2], se[3]], [se[0], se[3]]])
    return bbox_poly, ext_gdf  # Return the Polygon and the GeoDataFrame


def create_bbox(extent_shape_file_path: Path, bb_file_path: Path, create_bbox_shp: bool = True):
    """Create a bounding box shape file.

    This function takes a shape file, creates a bounding box around the features and returns this box as a
    GeoDataFrame. Additionally, returns the original extent as a GeoDataFrame.

    Args:
        extent_shape_file_path : Path
            Path to the shape file for which to create a bounding box (in EPSG2056).
        bb_file_path : Path
            Full path to the bounding box shape file to be written (in EPSG2056).
        create_bbox_shp: bool
            True if bounding box shape file should be written to file.

    Returns:
        bbox_gdf : GeoDataFrame
            Bounding box as a GeoDataFrame (in EPSG2056).
        ext_gdf : GeoDataFrame
            The GeoDataFrame with the data from the extent shape file (in EPSG2056)
    """

    # Create the bounding box Polygon
    bbox_poly, ext_gdf = create_bbox_gdf(extent_shape_file_path)
    # Create the bounding shape in a GeoDataFrame
    bbox_gdf: GeoDataFrame = gpd.GeoDataFrame(pd.DataFrame(['p1'], columns=['geom']),
                                              crs={'init': 'epsg:2056'},
                                              geometry=[bbox_poly])
    if create_bbox_shp:
        # This writes the bounding box as a shape file
        bbox_gdf.to_file(str(bb_file_path))
    return bbox_gdf, ext_gdf


def netcdf_to_dataset(netcdf_file_path: Path) -> xr.Dataset:
    """Reads a netCDF file into an xarray dataset

    Args:
        netcdf_file_path : Path
            Path to the netCDF file. It has to be in CH1903+/LV95 (EPSG=2056).

    Returns:
        xds : xr.Dataset
            The netCDF data as a netCDF4 dataset in CH1903+/LV95 (EPSG=2056).

    """

    # Read in the netCDF file into an xarray Dataset
    xds: xr.Dataset = xr.open_dataset(netcdf_file_path)
    # Define the dimensions
    xds.rio.set_spatial_dims(x_dim="E", y_dim="N", inplace=True)
    # Set the CRS projection EPSG to 2056, since this is the one the netCDF files are in
    xds.rio.write_crs("EPSG:2056", inplace=True)
    # delete the attribute 'grid_mapping' to prevent an error
    vars_list = list(xds.data_vars)
    try:
        for var in vars_list:
            del xds[var].attrs['grid_mapping']
    except KeyError:
        logger.exception(f"KeyError {xds}")
    return xds


def dataset_to_netcdf(xds_to_write: xr.Dataset, netcdf_file_path: Path, catchment: str):
    """ Writes xarray dataset to netCDF

    Args:
        xds_to_write : xr.Dataset
            xarray Dataset to write.
        netcdf_file_path : Path
            netCDF file path to write to.
        catchment: Catchment name

    """
    try:
        os.mkdir(f"{netcdf_file_path.parent.parent}/out/{catchment}")
    except FileExistsError:
        pass
    # Write the netCDF file
    xds_to_write.to_netcdf(
        f"{netcdf_file_path.parent.parent}/out/{catchment}/{netcdf_file_path.stem}_{catchment}_clipped{netcdf_file_path.suffix}",
        mode="w")


def netcdf_clipper(netcdf_file_path: Path, extent_file_path: Path):
    """

    Args:
        netcdf_file_path: Path to netCDF file to clip.
        extent_file_path: Path to extent to which netCDF file should be clipped.
    """
    logger.debug("Trying to call extent.sh...")
    rcode = subprocess.call(['raven_tools/extent.sh', str(extent_file_path), str(netcdf_file_path)])
    logger.debug(f"extent.sh executed with return code: {rcode}")


def netcdf_pet_hamon(netcdf_file_path: Path, name_pattern: dict[str, str]):
    """Calculates PET with Hamon approach and adds values to netCDF file.
    %TODO: Check how the name_pattern works and document it properly
    Args:
        name_pattern : dict[str, str]
            Contains the pattern for the cdf_out_path variable.
        netcdf_file_path : Path
            Path to the netCDF file to calculate PET from

    """
    cdf_out_path = Path(str(netcdf_file_path).replace(list(name_pattern.keys())[0], list(name_pattern.values())[0]))
    # Copy the clipped file to a new file so as not to change the original file
    shutil.copyfile(netcdf_file_path, cdf_out_path)
    # Read the clipped netCDF file into a netCDF Dataset
    # cdf_dataset_in = Dataset(netcdf_file_path, format="NETCDF4")
    # Create the output file in append mode
    cdf_dataset_out: netCDF4.Dataset = Dataset(cdf_out_path, "r+", format="NETCDF4")
    # Create a new variable 'PET' in the output file, according to the 'TabsD' variable in the input file
    pet = cdf_dataset_out.createVariable('PET', np.float32, fill_value=-999.99, dimensions=('time', 'N', 'E'))
    pet.units = 'mm/d'
    pet.long_name = 'daily potential evapotranspiration (Hamon)'
    pet.setncatts({'grid_name': "ch01r.swiss.lv95",
                   'version': "v1.2",
                   'prod_date': "2022-10-01",
                   'coordinates': "lon lat",
                   'grid_mapping': u"swiss_lv95_coordinates"})
    # Get the values of 'TabsD' of the dataset as a nested array
    tabsd = cdf_dataset_out['TabsD'][:]
    latitude = cdf_dataset_out['lat'][:, 1]
    pet_array = cdf_dataset_out['PET'][:]
    day_length = np.empty((366, 60, 36))

    for dx, day in np.ndenumerate(tabsd):
        sol_dec = 0.409 * np.sin(2. * np.pi / 366. * (dx[0] + 1) - 1.39)
        l_rad = latitude[dx[1]] * np.pi / 180
        s_angle = np.arccos(-np.tan(sol_dec) * np.tan(l_rad))
        dl_h = 24 / np.pi * s_angle
        pet_value = (dl_h / 12) ** 2 * exp(day / 16)
        day_length[dx[0], dx[1], dx[2]] = pet_value
    cdf_dataset_out['PET'][:] = day_length

    # Close the dataset and write the file
    cdf_dataset_out.close()


def netcdf_elevation(netcdf_filepath: Path):
    """Adds elevation variable to netCDF file.

    Args:
        netcdf_filepath: Path to netCDF file for which elevation variable should be added.
    """

    infile = Dataset(filename=netcdf_filepath, mode="r+", format="NETCDF4")

    try:
        ele = infile.createVariable('Elevation', np.float32, fill_value=-999.99, dimensions=('time', 'N', 'E'))
        ele.units = 'meters'
        ele.long_name = 'elevation above sea level'
        ele.setncatts({'grid_name': "ch01r.swiss.lv95",
                       'version': "v1.2",
                       'prod_date': "2022-10-01",
                       'coordinates': "lon lat",
                       'grid_mapping': u"swiss_lv95_coordinates"})
    except RuntimeError:
        pass
    latitude = infile['lat'][:, 1]


def pet_temp_monthly_ave(pet_filepath: Path, temp_filepath: Path):
    """Calculates monthly averages for temperature and PET from CSV values

    Args:
        pet_filepath: Path
            Path to PET CSV file.
        temp_filepath: Path
            Path to temperature CSV file.

    Returns:
        pet_monthly: Monthly PET means.
        temp_monthly: Monthly temperature means.
    """

    pet_monthly_from_order = pd.read_csv(pet_filepath, sep=";")

    pet_monthly_from_order['time'] = pd.to_datetime(pet_monthly_from_order['time'], format='%Y%m')
    pet_monthly_from_order['month'] = pet_monthly_from_order['time'].dt.month
    pet_monthly = pet_monthly_from_order.groupby(pet_monthly_from_order.time.dt.month)['ets150m0'].mean()

    temp_hourly_from_order = pd.read_csv(temp_filepath, sep=";")
    temp_hourly_from_order['time'] = pd.to_datetime(temp_hourly_from_order['time'], format='%Y%m%d%H')
    temp_monthly = temp_hourly_from_order.groupby(temp_hourly_from_order.time.dt.month)['tre200h0'].mean()

    return pet_monthly, temp_monthly


def resample_netcdf_monthly():
    """Resamples netCDF file to monthly values

    """
    import xarray as xr
    ds = xr.open_dataset(
        "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/merged/RhiresD_ch01h.swiss.lv95_198101010000_202012310000_Dischmabach_clipped.nc")
    ds_resampled = ds.resample(time='m').mean()
    netcdf_file_path = "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/merged/resampled.nc"
    dataset_to_netcdf(ds_resampled, Path(netcdf_file_path), catchment="Dischmabach")

    monthly_data = ds.resample(freq='m', dim='time', how='mean')


def nc_merge(start_year: int, end_year: int, forcing_dir: Path, forcing_prefix: str, catchment_name: str = "Broye",
             legacy_mode: bool = False):
    """Merges multiple netCDF files into one.

    Args:
        legacy_mode: bool
            True if you want to use legacy mode
        catchment_name: str
            Name of catchment.
        forcing_prefix: str
            Prefix to meteo forcings.
        start_year : int
            Start year, e.g. first file to merge.
        end_year : int
            End Year, e.g. last file to merge.
        forcing_dir : Path
            Root directory where forcing files are located

    """
    if legacy_mode:
        logger.debug("Trying to call nc_combine.sh...")
        rcode = subprocess.call(
            ['raven_tools/nc_combine.sh', str(start_year), str(end_year), str(forcing_dir), catchment_name])
        logger.debug(f"nc_combine.sh executed with return code: {rcode}")
    else:
        logger.info("Trying to call nc_combine.sh...")
        rcode = subprocess.call(
            ['raven_tools/cdo_mergetime.sh', str(start_year), str(end_year), str(forcing_prefix), str(forcing_dir)])
        logger.debug(f"cdo_mergetime.sh executed with return code: {rcode}")


def create_grid(netcdf_filepath: Path, bounding_box_filename: Path, out_path: Path, forcing_name: str, start_year,
                export_shp: bool = True) -> GeoDataFrame:
    """Creates a grid GeoDataFrame and optionally exports to shape file

    Args:
        netcdf_filepath : Path
            Path to netCDF file whose grid should be used (in EPSG2056)
        bounding_box_filename : Path
            Path to bounding box shape file (in EPSG2056)
        out_path : Path
            If export_shp=True, write shape file to this Path
        forcing_name : str
            Used to select column in netCDF data
        start_year :
            Used to slice data so as to save computation time
        export_shp : bool
            Set to True if you want to export the grid into a shape file (with EPSG2056)

    Returns:
        grid : GeoDataFrame
            GeoDataFrame containing the grid in EPSG2056

    """
    # Read in the bounding box shape file from the clipping folder as a GeoPandas DataFrame
    logger.debug(f"bbox_filename: {bounding_box_filename}")

    # Read in the clipped netCDF file into a xarray Dataset
    ds: Dataset = xr.open_dataset(netcdf_filepath, engine="netcdf4")

    # Since we're only interested in the grid (not the actual cell values), only use 1 day
    # Try to use 0-1 instead of strings to avoid having to use concrete date values.
    # ds = ds.sel(time=slice('1962-01-01', '1962-01-02'))
    start_date = f"01-01-{start_year}"
    start_date = datetime.strptime(str(start_date), "%d-%m-%Y")
    ds = ds.sel(time=slice(start_date, (start_date + timedelta(days=1))))
    # Select the actual data value column
    column_name = re.findall("[A-Za-z]+", forcing_name)[0]
    xarr: xr.DataArray = ds[column_name]
    # Convert Dataset to DataFrame and reset the index
    df: pd.DataFrame = xarr.to_dataframe().reset_index()
    # Convert the DataFrame to a GeoDataFrame, using the northing and easting provided by the netCDF. Furthermore,
    # change the projection to LV95
    data_column = pd.Series(df[column_name])
    combi_catchment_grid: GeoDataFrame = gpd.GeoDataFrame(data_column,
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
    grid_cols = ['row', 'col', 'cell_id', 'polygons', 'area']
    grid: GeoDataFrame = gpd.GeoDataFrame(columns=grid_cols, geometry='polygons', crs="EPSG:2056")
    # Loop over each cell and create the corresponding Polygon, then append it to the Polygon list
    # Also write the cell id from cell_id=i_row*n_columns + i_column
    for ix, x in enumerate(cols):
        for iy, y in enumerate(rows):
            polygons.append(Polygon([(x, y), (x + wide, y), (x + wide, y + length), (x, y + length)]))
            cid = int(iy * len(cols) + ix)
            cell_id.append(f"{str(cid)}")

    # Use the polygon list in the GeoDataFrame and set the projection accordingly
    grid["polygons"], grid["cell_id"] = polygons, cell_id

    if export_shp:
        # Export the grid to a shape file
        grid.to_file(str(str(out_path) + ".shp"))
    logger.debug(f"grid weights shp file written: {out_path}")
    return grid


def create_overlay(grd: GeoDataFrame, ctm_gdf: GeoDataFrame) -> GeoDataFrame:
    """Overlays a GeoDataFrame over another to create overlay Polygons

    Overlays two GeoDataFrame over each other and returns a new GeoDataFrame for the mode 'intersection'.

    Args:
        grd : GeoDataFrame
            Grid as given by the netCDF file in EPSG2056
        ctm_gdf: GeoDataFrame
            Catchment GDF

    Returns:
        res_u : GeoDataFrame
            Grid cells within the catchment area in EPSG2056
    """

    ctm_gdf.set_crs(epsg="2056"), grd.set_crs(epsg='2056')
    ctm_gdf.to_crs(32632), grd.to_crs(32632)
    res_u: GeoDataFrame = ctm_gdf.overlay(grd, how='intersection')
    # res_u.set_index("cell_id", inplace=True)
    res_u.to_crs(2056)
    return res_u


def create_dask_overlay(underlay: gpd.GeoDataFrame, overlay: gpd.GeoDataFrame, keep_geom_type=True,
                        overlay_type='intersection'):
    """Overlays a GeoDataFrame over another to create overlay Polygons

    Overlays two GeoDataFrame over each other and returns a new GeoDataFrame for the mode 'intersection'.

    Args:
        underlay : GeoDataFrame
            GeoDataFrame to use as underlying, in EPSG:2056
        overlay: GeoDataFrame
            GeoDataFrame to use as overlying, in EPSG:2056
        keep_geom_type: bool
            Used as parameter to the gpd.overlay function
        overlay_type: str
            Type of overlay to be used, e.g. 'intersection' or 'difference'

    Returns:
        res_u : GeoDataFrame
            GeoDataFrame after overlay operation in EPSG:2056
    """
    overlay.set_crs(2056)
    underlay.set_crs(2056)
    overlay.to_crs(32632)
    underlay.to_crs(32632)
    # overlay = overlay.explode()
    # underlay = underlay.explode()
    # res_u: dgpd.GeoDataFrame = ctm_gdf.overlay(grd, how='intersection', keep_geom_type=keep_geom_type)
    res_u = gpd.overlay(underlay, overlay, overlay_type, keep_geom_type)
    # res_u.set_index("cell_id", inplace=True)
    res_u.to_crs(2056)
    return res_u


def calc_relative_area(hru_gdf: GeoDataFrame, hru_short_name: str) -> GeoDataFrame:
    """Calculates the relative area of each polygon in a GeoDataFrame.

    Calculates the relative area of each polygon in a GeoDataFrame with EPSG=2056, writes it into a new column and
    returns the GeoDataFrame with EPSG=2056.

    Args:
        hru_gdf : GeoDataFrame
            GeoDataFrame in EPSG=2056
        hru_short_name: str
            Short name of HRU to be used as column name prefix

    Returns:
        gdf : GeoDataFrame
            GeoDataFrame with relative areas of each polygon, in EPSG:2056

    """
    area_column_name: str = f"{hru_short_name}_area_rel"

    # Set the CRS to WGS84 to preserve areas
    hru_gdf = hru_gdf.to_crs(32632)

    # For each intersected feature, compute the area and write into a new field
    hru_gdf[area_column_name] = (hru_gdf["geometry"].area * 1000000) / np.sum(hru_gdf.geometry.area * 1000000)
    # Convert the area to float
    # area_sum: float = float(hru_gdf["area"].sum())
    # Re-project back into the LV95 projection and export to a shape file
    hru_gdf: GeoDataFrame = hru_gdf.to_crs(2056)
    return hru_gdf


def create_grid_data_dict(grd: GeoDataFrame, glacier: bool = False) -> dict:
    """Loops over each grid cell and extracts the grid weights.

    Args:
        grd : GeoDataFrame
            Grid as derived from the netCDF file
        glacier : bool
            True if glacier

    Returns:
        data_to_write : dict
            Dict with the relative areas/grid weights of each cell.

    """
    # Loop over each intersected feature and write the relative area (compared with the total catchment area) into a new
    # field.
    data_to_write = {}
    if glacier:
        column_names = ["non_gla_area_rel", "gla_area_rel"]
        hru_base_counter = 1
        for c in column_names:
            area_rel = grd[c][~grd[c].isna()]
            data_to_write[f'{hru_base_counter}'] = area_rel
            hru_base_counter += 1
    return data_to_write


def write_grid_data_to_file(gdf: gpd.GeoDataFrame, grid_weights_file_path: Path):
    """Writes grid weight data to rvt file.

    Args:
        gdf: GeoDataFrame
            GeoDataFrame containing the grid data to write.
        grid_weights_file_path: Path
            Path to grid weights file to be written to.
    """
    col_list = list(gdf.columns)
    col_list.remove('polygons')
    gdf.fillna(0, inplace=True)
    with open(grid_weights_file_path, 'w') as ff:
        ff.write(':GridWeights                     \n')
        ff.write('   #                                \n')
        ff.write('   # [# HRUs]                       \n')
        ff.write(f'   :NumberHRUs       {len(col_list)}            \n')  # Currently for GR4J, there's 1 HRU
        ff.write(f"   :NumberGridCells  {len(gdf)}            \n")
        ff.write('   #                                \n')
        ff.write('   # [HRU ID] [Cell #] [w_kl]       \n')
        for hru_num, col_name in enumerate(col_list):
            lst = [f"   {hru_num + 1}   {index}   {weight}\n" for (index, weight) in zip(gdf.index, gdf[col_name])]
            ff.writelines(lst)
        ff.write(':EndGridWeights \n')


def write_grid_data_to_file(weights: pd.DataFrame, grid_weights_file_path: Path):
    """Writes grid weight data to rvt file.

    Args:
        weights: DataFrame
            DataFrame containing the grid data to write.
        grid_weights_file_path: Path
            Path to grid weights file to be written to.
    """

    with open(grid_weights_file_path, 'w') as ff:
        ff.write(':GridWeights                     \n')
        ff.write('   #                                \n')
        ff.write('   # [# HRUs]                       \n')
        ff.write(
            f'   :NumberHRUs       {len(weights.hru_id.unique())}            \n')  # Currently for GR4J, there's 1 HRU
        ff.write(f"   :NumberGridCells  {len(weights)}            \n")
        ff.write('   #                                \n')
        ff.write('   # [HRU ID] [Cell #] [w_kl]       \n')
        lst = [f"   {hru_id}   {cell_id}   {weight}\n" for (hru_id, cell_id, weight) in
               zip(weights.hru_id, weights.cell_id, weights.weight)]
        ff.writelines(lst)
        ff.write(':EndGridWeights \n')


# def copy_rel_area_from_union_to_grid(res_union: GeoDataFrame, grid: GeoDataFrame, hru_short_name: str) -> GeoDataFrame:
#     """Takes grid weights from a union GeoDataFrame and writes the to the grid GeoDataFrame.
#
#     Args:
#         res_union : GeoDataFrame
#             GeoDataFrame containing the grid cells within the catchment.
#         grid : GeoDataFrame
#             Grid GeoDataFrame as derived from netCDF file
#
#     Returns:
#         grd : GeoDataFrame
#             Grid GeoDataFrame with grid weights
#
#     """
#
#     # Loop over the union GeoDataFrame, take relative area/grid weight and write it to corresponding cell in grid
#     # GeoDataFrame.
#     grid_column_name: str = f"{hru_short_name}_area_rel"
#     for index, row in res_union.iterrows():
#         cell_id_old_value = res_union.at[index, "cell_id"]
#         area_rel_old_value = res_union.at[index, grid_column_name]
#         ind = grid[grid['cell_id'] == cell_id_old_value].index.tolist()
#         grid.at[ind[0], grid_column_name] = area_rel_old_value
#     return grid


def camels_to_rvt(data_dir, gauge_id, gauge_short_code, start_date="2000-01-01", end_date="2000-12-31"):
    """Reads CAMELS CSV discharge file and creates RAVEN .rvt file.

    Reads a daily CAMELS discharge CSV file, with the file path read from a global parameter and converts it to a daily
    discharge .rvt file compatible with RAVEN. It takes start and end date as arguments to only use a selected date
    range.

    Args:
        data_dir:
        gauge_id:
        gauge_short_code:
        start_date: str
        end_date: str

    """
    logger.debug("Entered function camels_to_rvt.")
    # Read in the discharge data from .txt file
    camels_filename = f"CAMELS_CH_obs_based_{gauge_id}.txt"
    logger.debug(f"camels_filename = {camels_filename}")
    df_meteo: pd.DataFrame = pd.read_csv(Path(data_dir, "Discharge", camels_filename), sep=";")
    # Rename the column for easier keyboard typing
    df_meteo = df_meteo.rename(columns=config.variables.camels_column_names)

    # Convert time columns to datetime format
    df_meteo['date'] = pd.to_datetime(df_meteo['date'], format="%Y-%m-%d")
    logger.debug("Converted date column to datetime format")
    # Set the start time
    start_time = "0:00:00"
    # Subset according to start and end date
    df_meteo = subset_dataframe_time(df_meteo, start_date, end_date)
    # Uncomment the following line, if you want to get daily means. Otherwise, RAVEN will do it for you
    # df_meteo = df_meteo.resample('d', on='time').mean().dropna(how="all")
    # Invoke export_to_rvt_file to export
    out_filename = f"{gauge_short_code}_Q_{gauge_id}_daily.rvt"
    out_path = Path(data_dir, "Discharge", out_filename)
    export_to_rvt_file(start_date, start_time, df_meteo, out_path)


def subset_dataframe_time(dataframe: pd.DataFrame, start_date: str, end_date: str) -> pd.DataFrame:
    """Subsetting a dataframe using a time interval.

    Args:
        dataframe: pd.DataFrame
            DataFrame with data
        start_date: str
            Start date of subset
        end_date: str
            End date of subset

    Returns:
        subset_dataframe: pd.DataFrame
            DataFrame with new start and end date.

    """
    # Date to string conversion
    logger.debug("Entered function subset_dataframe_time.")
    start_date = datetime.strptime(start_date, "%Y-%m-%d")
    end_date = datetime.strptime(end_date, "%Y-%m-%d")
    logger.debug("Did strptime on start_date and end_date.")
    # Create data interval mask. Offsetting by 1 day to include end date.
    # TODO: Is this offsetting still necessary, since changing to the CAMEL files?
    mask = dataframe['date'].between(start_date, end_date + pd.DateOffset(days=1), inclusive="left")
    logger.debug("Created data interval mask.")
    # Apply the mask
    subset_dataframe = dataframe[mask]
    logger.debug("Applied data interval mask on dataframe -> return result.")
    return subset_dataframe


def export_to_rvt_file(start_date, start_time, df, out_path):
    """Writes RVT file from DataFrame

    Args:
        start_date: str
            Start date of simulation
        start_time: str
            End date of simulation
        df: pd.DataFrame
            DataFrame to write
        out_path:
            Path to RVT file to be written

    """
    with open(out_path, 'w') as f:
        # print(rvt_filename)
        f.write(f":ObservationData\tHYDROGRAPH\t1\tm3/s\n{start_date}\t{start_time}\t1\t{len(df)}\n")
        # For gauged precipitation data
        df_as_string = df.to_string(justify="right", header=False, index=False,
                                    columns=['discharge'], na_rep="-1.2345")
        f.write(df_as_string)
        # f.write(df_as  _string)
        f.write("\n:EndObservationData")


def glacier_extent(ctm_gdf: GeoDataFrame, glacier_gdf: GeoDataFrame) -> tuple[GeoDataFrame, GeoDataFrame]:
    """Creates glaciated and non-glaciated area in a catchment from catchment extent and glacier extent.

    Args:
        ctm_gdf: GeoDataFrame
            Catchment extent.
        glacier_gdf: GeoDataFrame
            Glacier extent.

    Returns:
        ctm_glaciation: GeoDataFrame
            Glaciated area.
        ctm_non_glaciation: GeoDataFrame
            Non-glaciated area.

    """
    ctm_glaciation: GeoDataFrame = ctm_gdf.overlay(glacier_gdf, how='intersection')
    ctm_non_glaciation: GeoDataFrame = ctm_gdf.overlay(glacier_gdf, how='difference')
    # ctm_glaciation.set_crs(epsg='2056')
    # ctm_non_glaciation.set_crs(epsg='2056')
    # ctm_glaciation.set_index("cell_id")
    # ctm_glaciation.set_crs(epsg='2056')
    return ctm_glaciation, ctm_non_glaciation


def hru_extent_from_shp(ctm_shp: Path, hru_shp: Path) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Creates a HRU and non_HRU area in a catchment from shape files.

    Args:
        ctm_shp: Path
            Path to catchment shape file
        hru_shp: Path
            Path to HRU shape file

    Returns:
        ctm_hru_intersection: GeoDataFrame
            Area of HRU within catchment
        ctm_hru_difference: GeoDataFrame
            Area within catchment that is not HRU

    """
    ctm_gdf, hru_gdf = gpd.read_file(ctm_shp), gpd.read_file(hru_shp)
    drop_col_match = ['ezgnr', 'ar_1903', 'watr_nm', 'place', 'gid', 'pk_glacier', 'sgi-id', 'name', 'rl_0', 'rl_1',
                      'rl_2', 'rl_3', 'i_code', 'year_acq', 'year_rel', 'area_km2', 'length_km', 'masl_min', 'masl_med',
                      'masl_mean', 'masl_max', 'slope_deg', 'aspect_deg']
    ctm_gdf, hru_gdf = (ctm_gdf.drop(columns=[col for col in drop_col_match if col in list(ctm_gdf.columns)]),
                        hru_gdf.drop(columns=[col for col in drop_col_match if col in list(hru_gdf.columns)]))
    hru_gdf_explode = hru_gdf.explode(ignore_index=True)
    ctm_hru_intersection: GeoDataFrame = ctm_gdf.overlay(hru_gdf, how='intersection')
    ctm_hru_difference: GeoDataFrame = ctm_gdf.overlay(hru_gdf_explode, how='difference', keep_geom_type=True)
    return ctm_hru_intersection, ctm_hru_difference


def area_from_ratio_dem_props(area_ratios: dict, ctm_ch_id: str, data_dir: Path) -> dict:
    """Return dictionary with area per elevation band.

    Args:
        area_ratios: dict
            Dictionary with lower bounds of elevation bands as keys.
        ctm_ch_id:
            Catchment id.
        data_dir:
            Path to project data directory.

    Returns:
        area_ratios: dict
            Dictionary with lower bounds of elevation bands as keys and respective areas as values.
    """
    props = pd.read_csv(Path(data_dir, f"Catchment/hru_info.csv"))
    non_gla_area = float(props[props['Ctm'] == ctm_ch_id]['NonGlaArea'])
    gla_area = float(props[props['Ctm'] == ctm_ch_id]['GlaArea'])
    total_area = non_gla_area + gla_area
    for lower_bound, ratio in area_ratios.items():
        area_ratios[lower_bound] = ratio * non_gla_area
    return area_ratios


def dict_to_txt(dict: dict, ctm_ch_id: str, data_dir: Path):
    """Writes elevation band areas from dictionary to text file.

    Args:
        dict: dict
            Elevation band areas to write
        ctm_ch_id:
            Catchment id.
        data_dir:
            Path to data directory.
    """

    with open(Path(data_dir, f"DEM/hbv/non_glacier/elevation_band_areas_{ctm_ch_id}.txt"), 'w') as f:
        for key, value in dict.items():
            f.write(f"{key}:{value}\n")


def area_ratio(catchment_filepath: Path, glacier_shape_path: Path):
    """Calculates the relative glaciated area of a catchment

    Args:
        catchment_filepath: Path
            Path to catchment extent shape file.
        glacier_shape_path: Path
            Path to glacier extent shape file.

    Returns:
        ctm_total_area: float
            Total catchment area
        non_glaciation_area:
            Total non-glaciated area within catchment.
        glaciation_area (optional):
            Total glaciated area within catchment, if any.

    """
    ctm_gdf: GeoDataFrame = gpd.read_file(catchment_filepath)
    ctm_gdf.set_crs(epsg="2056")
    glacier_gdf: GeoDataFrame = gpd.read_file(glacier_shape_path)
    glacier_gdf.set_crs(epsg='2056')
    glaciated_area_gdf, non_glaciated_area_gdf = glacier_extent(ctm_gdf=ctm_gdf, glacier_gdf=glacier_gdf)

    ctm_gdf_area = ctm_gdf.to_crs(epsg='32632')
    non_glaciated_area_gdf_area = non_glaciated_area_gdf.to_crs(epsg='32632')
    glaciated_area_gdf_area = glaciated_area_gdf.to_crs(epsg='32632')
    ctm_total_area: float = ctm_gdf_area["geometry"].area.sum() / 10 ** 6
    non_glaciation_area = non_glaciated_area_gdf_area["geometry"].area.sum() / 10 ** 6

    if glaciated_area_gdf_area.empty:
        return ctm_total_area, non_glaciation_area
    else:
        glaciation_area = glaciated_area_gdf_area["geometry"].area.sum() / 10 ** 6
        return ctm_total_area, non_glaciation_area, glaciation_area


# def dfdf(catchment_filepath: Path, glacier_shape_path: Path, dem_filepaths: list[Path],
#          dem_out_path: Path, ctm_ch_id: str):
#     ctm_gdf: GeoDataFrame = gpd.read_file(catchment_filepath)
#     ctm_gdf.set_crs(epsg="2056")
#     glacier_gdf: GeoDataFrame = gpd.read_file(glacier_shape_path)
#     glacier_gdf.set_crs(epsg='2056')
#     glaciated_area_gdf, non_glaciated_area_gdf = glacier_extent(ctm_gdf=ctm_gdf, glacier_gdf=glacier_gdf)
#     non_glaciated_centroid = weighted_centroid(non_glaciated_area_gdf)
#     non_glaciation_area = non_glaciated_area_gdf["geometry"].area.sum()
#     ctm_area = ctm_gdf["geometry"].area.sum()
#
#     non_gla_height = hru_height(hru_area_gdf=non_glaciated_area_gdf.to_crs(epsg='2056'), dem_filepaths=dem_filepaths,
#                                 dem_out_dir=dem_out_path, ctm_ch_id=ctm_ch_id)
#
#     if glaciated_area_gdf.empty:
#         return glaciation_ratio, gla_height, non_glaciation_ratio, non_gla_height, Point(0, 0), non_glaciated_centroid
#     else:
#
#         glaciacted_centroid = weighted_centroid(glaciated_area_gdf)
#         glaciation_area = glaciated_area_gdf["geometry"].area.sum()
#         glaciation_ratio = (glaciation_area * 1000000) / (ctm_area * 1000000)
#         non_glaciation_ratio: float = (non_glaciation_area * 1000000) / (ctm_area * 1000000)
#
#         gla_height = hru_height(hru_area_gdf=glaciated_area_gdf.to_crs(epsg='2056'), dem_filepaths=dem_filepaths,
#                                 dem_out_dir=dem_out_path, ctm_ch_id=ctm_ch_id, glacier=True)
#         return glaciation_ratio, gla_height, non_glaciation_ratio, non_gla_height, glaciacted_centroid, non_glaciated_centroid


def weighted_centroid_vector(feature_gdf: GeoDataFrame) -> tuple[float, float]:
    """Calculates weighted centroid of a GeoDataFrame feature

    Args:
        feature_gdf: GeoDataFrame
            The feature whose centroid to calculate.

    Returns:
        A tuple with the X and Y coordinates of the centroid.

    """
    feature_gdf.set_crs(epsg='2056')
    feature_gdf.to_crs(epsg='32632', inplace=True)
    area = np.array(feature_gdf.geometry.area)
    centroids_x = np.array(feature_gdf.geometry.centroid.x)
    weighted_centroid_x = (np.sum(centroids_x * area) / np.sum(area))
    centroids_y = np.array(feature_gdf.geometry.centroid.y)
    weighted_centroid_y = (np.sum(centroids_y * area) / np.sum(area))
    weighted_centroid_x, weighted_centroid_y = crs_old_to_new(weighted_centroid_x, weighted_centroid_y,
                                                              32632, 4326)
    weighted_centroid_point = Point(weighted_centroid_x, weighted_centroid_y)
    return weighted_centroid_y, weighted_centroid_x


def weighted_centroid_raster_points(raster_points_gdf: GeoDataFrame) -> tuple:
    """Calculates weighted centroid for raster

    Args:
        raster_points_gdf: GeoDataFrame
            Raster points

    Returns:

    """
    raster_points_gdf_new = raster_points_gdf.to_crs(epsg='32632')
    raster_points_gdf_diss = raster_points_gdf_new.dissolve()
    centroid = raster_points_gdf_diss.centroid[0]
    return crs_old_to_new(lat_old=centroid.y, lon_old=centroid.x, epsg_old=32632, epsg_new=4326)


# def hru_height(hru_area_gdf: GeoDataFrame, dem_filepaths: list[Path], dem_out_dir: Path, ctm_ch_id, glacier=False):
#     dem_im = clip_dem_file(dem_filepaths)
#
#     try:
#
#         dem_mean = dem_clipped.mean(dim=["x", "y"], skipna=True).to_numpy()
#         mean_height = dem_mean.tolist()
#         return mean_height
#     except ValueError:
#         logger.exception("There has been an error clipping the DEM")
#         pass


# def clip_dem_file(hru_area_gdf, dem_filepaths, glacier, dem_out_dir, ctm_ch_id):
#     import rasterio as rio
#     from rasterio import merge
#     import rioxarray as rxr
#     from rasterio.io import MemoryFile
#
#     memfile = MemoryFile()
#     file_handler = [rio.open(row) for row in dem_filepaths]
#     rio.merge.merge(datasets=file_handler, dtype='float32',
#                     dst_path=memfile.name)
#     dataset = memfile.open(driver='GTiff')
#     dem_im = rxr.open_rasterio(dataset.name, masked=True).squeeze()
#
#     dem_clipped = dem_im.rio.clip(hru_area_gdf.geometry.apply(mapping))
#     if glacier:
#         dem_out_filepath = Path(dem_out_dir, "dem_clipped", f"dem_{ctm_ch_id}_glacier.tif")
#     else:
#         dem_out_filepath = Path(dem_out_dir, "dem_clipped", f"dem_{ctm_ch_id}.tif")
#     dem_clipped.rio.to_raster(dem_out_filepath)


def dem_mean(filepath: Path):
    """Calculates mean value of a DEM tif file

    Args:
        filepath: Path
            Path to DEM tif file

    Returns:
        mean:
            Mean value of DEM
    """

    dem_im = rxr.open_rasterio(filepath, masked=True).squeeze()
    dem_mean_raster = dem_im.mean(dim=["x", "y"], skipna=True).to_numpy()
    mean = dem_mean_raster.tolist()
    return mean


def dem_mean_rasterio(dem: Union[Dataset, xr.DataArray, list[Dataset]]):
    """Calculates mean value of a DEM Dataset or DataArray

    Args:
        dem: Union[xr.Dataset, xr.DataArray, list[Dataset]]
            DEM from xarray.

    Returns:
        mean:
            Mean value of DEM
    """

    dem_mean_raster = dem.mean(dim=["x", "y"], skipna=True).to_numpy()
    mean = dem_mean_raster.tolist()
    return mean


def elevation_bands(ctm_ch_id: str, data_dir: Path, catchment_filepath: Path, save_to_tif: bool = False,
                    glacier_shape_path=None, basic_grid=None, has_glacier: bool = False):
    """Calculates elevation bands of a catchment and returns the grid file with the elevation bands overlayed.

    Args:
        ctm_ch_id: str
            Catchment id.
        data_dir: Path
            Data directory.
        catchment_filepath: Path
            Catchment extent shape file path.
        save_to_tif: bool
            Set True if elevation bands should also be saved to tif files.
        glacier_shape_path (optional):
            Path to glacier extent shape file. Leave empty, if no glacier in catchment.
        basic_grid: GeoDataFrame
            The grid GeoDataFrame from the netCDF file. Elevation bands will be written to it.
        has_glacier: bool
            Set True if catchment has glacier.

    Returns:
        grid_weights_dgdf: dgpd.GeoDataFrame
            Grid weights Dask GeoDataFrame with elevation bands.

    """
    rio_non_glacier = rxr.open_rasterio(Path(data_dir, "DEM", "out", f"dem_{ctm_ch_id}.tif"), masked=True,
                                        chunks=True).squeeze()
    if has_glacier:
        grid_weights_dgdf = extract_elevation_band_from_rio_dem(dem=rio_non_glacier,
                                                                ctm_ch_id=ctm_ch_id,
                                                                data_dir=data_dir,
                                                                catchment_filepath=catchment_filepath,
                                                                save_to_tif=save_to_tif,
                                                                glacier_shape_path=glacier_shape_path,
                                                                basic_grid=basic_grid,
                                                                has_glacier=has_glacier)
    else:
        grid_weights_dgdf = extract_elevation_band_from_rio_dem(dem=rio_non_glacier,
                                                                ctm_ch_id=ctm_ch_id,
                                                                data_dir=data_dir,
                                                                catchment_filepath=catchment_filepath,
                                                                save_to_tif=save_to_tif,
                                                                basic_grid=basic_grid,
                                                                has_glacier=has_glacier)
    # with open(Path(data_dir, f"DEM/hbv/non_glacier/area_ratios_{ctm_ch_id}.txt"), 'w') as f:
    #     for key, value in rio_elevation_band_dict.items():
    #         f.write(f"{key}:{value}\n")
    return grid_weights_dgdf


def raster_to_polygon(raster_file_path: Path):
    """

    Args:
        raster_file_path:

    Returns:

    """
    df = load_dem_tif_to_dataframe(raster_file_path)
    df.dropna(subset=["alti"], inplace=True)
    df['points'] = gpd.points_from_xy(df.x, df.y)
    # df['poly'] = Polygon(((df['x']-5, df['y']-5), (df['x']-5, df['y']+5), (df['x']+5, df['y']+5), (df['x']+5, df['y']-5), (df['x']-5, df['y']-5)))
    x_pt_list = df['x']
    y_pt_list = df['y']
    df.reset_index(inplace=True)
    for i in range(len(df)):
        x = df.loc[i, 'x']
        y = df.loc[i, 'y']
        df.loc[i, 'poly'] = Polygon(((x - 5, y - 5), (x - 5, y + 5), (x + 5, y + 5), (x + 5, y - 5), (x - 5, y - 5)))
    df_mean = pd.DataFrame(columns=['alti'])
    df_mean['alti'] = df.mean().loc['alti']
    boundary = gpd.GeoSeries(unary_union(df['poly']))
    gl = gpd.GeoDataFrame.from_file(
        Path("/home/sirian/Applications/Hydrology/RAVEN/data/glaciers", "SGI_2016_glaciers.shp"))
    ct = gpd.GeoDataFrame.from_file(
        Path("/home/sirian/Applications/Hydrology/RAVEN/data/Catchment/reproject_2056/CH-0105.shp"))
    geodf = gpd.GeoDataFrame(geometry=boundary, crs='epsg:2056')
    geodf.set_crs('epsg:2056')
    gl.set_crs('epsg:2056')
    ct.set_crs('epsg:2056')
    geodf.to_crs('epsg:32632', inplace=True)
    gl.to_crs('epsg:32632', inplace=True)
    ct.to_crs('epsg:32632', inplace=True)
    u = geodf.overlay(gl, how='difference')
    p = u.overlay(ct, how='intersection')
    p['area'] = p.area
    p.to_crs('epsg:2056', inplace=True)
    p['centroid_x'] = p.centroid.x
    p['centroid_y'] = p.centroid.y
    p['lower_bound'] = get_lower_band_limit_from_filepath(str(raster_file_path))
    p['alti'] = df.mean().loc['alti']
    geodf['alti'] = df.mean().loc['alti']
    geodf['centroid'] = boundary.centroid
    return p


def rio_elevation_band_dict_to_txt(rio_elevation_band_dict, ctm_id, df_non_glacier, df_glacier):
    """Writes elevation band data from rasterio to text file.

    Args:
        rio_elevation_band_dict: dict
            Contains elevation bands as rasterio rasters.
        ctm_id: str
            Catchment id.
        df_non_glacier: GeoDataFrame
            Non-glaciated extent.
        df_glacier: GeoDataFrame
            Glaciated extent.
    """
    ratio_dict = {}
    it = rio_elevation_band_dict.items()
    for key, elevation_band in it:
        df_band = load_rio_dataset_to_dataframe(elevation_band)
        ratio_dict[key] = ((df_non_glacier.count().y - df_band[df_band['alti'].isna()].count().y) / \
                           ((df_non_glacier.count().y - df_non_glacier[df_non_glacier['alti'].isna()].count().y) + \
                            (df_glacier.count().y - df_glacier[df_glacier['alti'].isna()].count().y)))
    with open(f"/media/mainman/Work/RAVEN/data/DEM/hbv/bands_ratio_{ctm_id}.txt", mode="w") as f:
        f.write(json.dumps(ratio_dict))


def round_up(x: float) -> int:
    """Rounds up number.

    Args:
        x: float
            Number to round up.

    Returns:
        rounded_up: int
            Rounded up number.
    """
    rounded_up = int(math.ceil(x / 100.0)) * 100
    return rounded_up


def round_down(x: float) -> int:
    """Rounds down number.

    Args:
        x: float
            Number to round down.

    Returns:
        rounded_down: int
            Rounded down number.
    """

    rounded_down = int(math.floor(x / 100.0)) * 100
    return rounded_down


def load_dem_tif_to_dataframe(dem_tif_filepath: Path) -> pd.DataFrame:
    """Creates DataFrame from raster tif file.

    Args:
        dem_tif_filepath: Path
            Path to tif file to convert into dataframe.

    Returns:
        dem_df: pd.DataFrame
            DataFrame with DEM data.
    """

    rio_dem = rxr.open_rasterio(dem_tif_filepath, masked=True).squeeze()
    rio_dem.name = "alti"
    dem_df = rio_dem.to_dataframe().reset_index()
    return dem_df


def load_dem_tif_to_dask_dataframe(dem_tif_filepath: Path) -> dd.DataFrame:
    """Creates Dask DataFrame from raster tif file.

    Args:
        dem_tif_filepath: Path
            Path to tif file to convert into Dask DataFrame.

    Returns:
        dem_ddf: dd.DataFrame
            Dask DataFrame with DEM data.
    """

    rio_dem = rxr.open_rasterio(dem_tif_filepath, masked=True, chunks=True).squeeze()
    rio_dem.name = "alti"
    dem_ddf = rio_dem.to_dask_dataframe()
    return dem_ddf


def load_rio_dataset_list_to_dataframe(rio_dataset_list) -> list[pd.DataFrame]:
    """Creates list containing DataFrames, each one of which converted from a rasterio DataSet.

    Args:
        rio_dataset_list: list
            List with rasterio DataSets to convert to DataFrames.
    Returns:
        df_list: list[pd.DataFrame]
            List with DataFrames
    """

    df_list: list[pd.DataFrame] = []
    for idx, i in enumerate(rio_dataset_list):
        i.name = "alti"
        df_list.append(i.to_dataframe().reset_index())
    return df_list


def load_rio_dataset_to_dataframe(rio_dataset) -> pd.DataFrame:
    """Convert a single rasterio DataSet to a DataFrame

    Args:
        rio_dataset:
            Rasterio DataSet
    Returns:
        df: pd.DataFrame
            Converted DataFrame
    """

    rio_dataset.name = "alti"
    df: pd.DataFrame = rio_dataset.to_dataframe().reset_index()
    return df


def load_rio_dataset_to_dask_dataframe(rio_dataset: xr.DataArray) -> dd.DataFrame:
    """Convert a single xarray DataArray to a Dask DataFrame

    Args:
        rio_dataset: xr.DataArray
            DataArray to convert.
    Returns:
        ddf: dd.DataFrame
            Converted Dask DataFrame
    """
    rio_dataset.name = "alti"
    ddf: dd.DataFrame = rio_dataset.to_dask_dataframe()
    return ddf


def extract_elevation_band_from_rio_dem(dem: xr.DataArray, ctm_ch_id: str, data_dir: Path, catchment_filepath: Path,
                                        save_to_tif: bool = False, glacier_shape_path=None,
                                        basic_grid=None, has_glacier: bool = False) -> \
        tuple[dict[Any, Any], dict, Any]:
    """Extracts elevation bands from DataArray and writes them to a grid file.

    Args:
        dem: xr.DataArray
            DEM elevation data.
        ctm_ch_id: str
            Catchment id.
        data_dir: Path
            Data directory
        catchment_filepath: Path
            Path to catchment extent shape file
        save_to_tif: bool
            Set True if elevation bands should be saved to tif files.
        glacier_shape_path: Path
            Path to glacier extent shape file
        basic_grid: gpd.GeoDataFrame
            Grid to which elevation band data is added.
        has_glacier: bool
            Set True if catchment has glaciers

    Returns:
        new_grid_dgdf: dgpd.GeoDataFrame
            Grid Dask GeoDataFrame to which elevation band data has been added.
    """

    catchment_extent = gpd.read_file(catchment_filepath)
    dem_clipped_to_ctm = dem.rio.clip(catchment_extent.geometry, crs="epsg:2056")
    lower = round_down(float(dem_clipped_to_ctm.min().compute()))
    upper = round_up(float(dem_clipped_to_ctm.max().compute()))
    elevation_band_limit_list = list(range(lower, upper, 100))
    basic_grid_dgdf = dgpd.from_geopandas(basic_grid, npartitions=2)
    if has_glacier:
        gla_extent, non_gla_extent = hru_extent_from_shp(ctm_shp=catchment_filepath,
                                                         hru_shp=glacier_shape_path)
        ctm_without_glacier_extent = create_dask_overlay(underlay=catchment_extent, overlay=gla_extent,
                                                         keep_geom_type=False, overlay_type='difference')
        dem_clipped_to_ctm_without_glacier = dem.rio.clip(ctm_without_glacier_extent.geometry, crs="epsg:2056")
        dem_df = dem_clipped_to_ctm_without_glacier.to_dataframe(name='alti').reset_index()
        dem_gdf = gpd.GeoDataFrame(dem_df, geometry=gpd.points_from_xy(dem_df.x, dem_df.y, crs=2056), crs=2056)
        dem_dgdf = dgpd.from_geopandas(dem_gdf, npartitions=8)
        new_grid_dgdf = basic_grid_dgdf
        new_grid_sjoin = new_grid_dgdf.sjoin(df=dem_dgdf, how='inner', predicate='intersects')
        for limit in elevation_band_limit_list:
            new_grid_dgdf = new_grid_dgdf.set_index('cell_id').join(
                other=new_grid_sjoin.drop(
                    columns=['polygons', 'index_right', 'y', 'x', 'band', 'spatial_ref']).rename(
                    columns={'alti': str(limit)}).where(
                    (new_grid_sjoin.alti >= limit) & (new_grid_sjoin.alti < limit + 100)).dropna().groupby(
                    'cell_id').count())
        gla_extent = dgpd.from_geopandas(gla_extent.explode().reset_index(), npartitions=1)
        g = dem.rio.clip(gla_extent.explode().geometry, crs="epsg:2056")
        g_df = g.to_dataframe(name='gla').reset_index()
        g_gdf = gpd.GeoDataFrame(g_df, geometry=gpd.points_from_xy(g_df.x, g_df.y, crs=2056), crs=2056)
        g_dgdf = dgpd.from_geopandas(g_gdf, npartitions=8)
        gla_grid_dgdf = basic_grid_dgdf.sjoin(df=g_dgdf, how='inner', predicate='intersects')
        new_grid_dgdf = new_grid_dgdf.join(
            other=gla_grid_dgdf.drop(columns=['polygons', 'index_right', 'y', 'x', 'band', 'spatial_ref']).rename(
                columns={'alti': 'gla'}).dropna().groupby('cell_id').count())
        # elevation_band_ratio_dict = {limit: (band.sum().values / dem_clipped_to_ctm.sum().values) for
        #                              (limit, band) in
        #                              elevation_band_dict.items()}
        if save_to_tif:
            for limit in elevation_band_limit_list:
                dem_clipped_to_ctm_without_glacier.where((dem_clipped_to_ctm_without_glacier.values >= limit) & (
                        dem_clipped_to_ctm_without_glacier.values < limit)).rio.to_raster(Path(data_dir,
                                                                                               f"DEM/hbv/non_glacier/dem_{ctm_ch_id}_{limit}_{limit + 99}.tif"))
    else:
        dem_df = dem_clipped_to_ctm.to_dataframe(name='alti').reset_index()
        dem_gdf = gpd.GeoDataFrame(dem_df, geometry=gpd.points_from_xy(dem_df.x, dem_df.y, crs=2056), crs=2056)
        dem_dgdf = dgpd.from_geopandas(dem_gdf, npartitions=8)
        new_grid_dgdf = basic_grid_dgdf
        new_grid_sjoin = new_grid_dgdf.sjoin(df=dem_dgdf, how='inner', predicate='intersects')
        for limit in elevation_band_limit_list:
            new_grid_dgdf = new_grid_dgdf.set_index('cell_id').join(
                other=new_grid_sjoin.drop(
                    columns=['polygons', 'index_right', 'y', 'x', 'band', 'spatial_ref']).rename(
                    columns={'alti': str(limit)}).where(
                    (new_grid_sjoin.alti >= limit) & (new_grid_sjoin.alti < limit + 100)).dropna().groupby(
                    'cell_id').count())
        if save_to_tif:
            for limit in elevation_band_limit_list:
                dem_clipped_to_ctm.where((dem_clipped_to_ctm.values >= limit) & (
                        dem_clipped_to_ctm.values < limit)).rio.to_raster(Path(data_dir,
                                                                               f"DEM/hbv/non_glacier/dem_{ctm_ch_id}_{limit}_{limit + 99}.tif"))
    return new_grid_dgdf


def crs_old_to_new(lat_old, lon_old, epsg_old: int, epsg_new: int):
    """Converts lat/lon coordinates to a new EPSG

    Args:
        lat_old:
        lon_old:
        epsg_old: int
        epsg_new: int

    Returns:
        lat_new:
        lon_new:

    """
    from pyproj import Transformer
    transformer = Transformer.from_crs(f"EPSG:{epsg_old}", f"EPSG:{epsg_new}")
    lat_new, lon_new = transformer.transform(yy=lat_old, xx=lon_old)
    return lat_new, lon_new


def create_elevation_band_tif_list(base_path: str, type: str) -> list:
    """Creates list of all elevation band limits from tif files

    Args:
        base_path: Path
        type: str
            Tif file type, e.g. 'aspect', 'slopes', etc.

    Returns:
        band_lower_list: list
            Lower limits of elevations bands.
    """

    asp_list = glob.glob(base_path + f'/{type}/dem_{ctm}_*_{type}.tif')
    band_lower_list = []
    for a in asp_list:
        rs = re.search(r'(?:_)(\d{4}?)(?:_)', a).group()
        band_lower_list.append(re.search(r'(\d{4})', re.search(r'(?:_)(\d{4}?)(?:_)', a).group()).group())
    band_lower_list.sort()
    return band_lower_list


def create_elevation_band_tif_list_hbv(data_dir: Path, tif_type: str, ctm_ch_id: str):
    """Creates list of all elevation band limits from tif files

    Args:
        data_dir: Path
            Data directory
        tif_type: str
            Tif file type, e.g. 'aspect', 'slopes', etc.
        ctm_ch_id: str
            Catchment id.

    Returns:
        band_lower_list: list
            Lower limits of elevation bands.
    """
    asp_list = glob.glob(str(data_dir) + f'/DEM/hbv/{tif_type}/dem_{ctm_ch_id}_*.tif')
    band_lower_list = []
    for a in asp_list:
        # rs = re.search(r'(?:_)(\d+?)(?:_)', str(a)).group()
        band_lower_list.append(int(re.search(r'(\d+)', (re.search(r'(?:_)(\d+)(?:_)', a).group())).group()))
    band_lower_list.sort()
    return band_lower_list


def get_lower_band_limit_from_filepath(filepath):
    """Returns lower elevation band limit from a file path.

    Args:
        filepath: Path
            Path to elevation band tif file.

    Returns:
        lower_band_limit: str
            Lower elevation band limit.
    """
    lower_band_limit = re.search(r'(\d{4})', re.search(r'(?:_)(\d{4}?)(?:_)', filepath).group()).group()
    return lower_band_limit


def add_elevation_to_netcdf(base_path: Path, ctm_ch_id: str, forcing_type: str):
    """Add elevation variable and data to netCDF file.

    Args:
        base_path: Path
            Directory where RAVEN folder resides.
        ctm_ch_id: str
            Catchment id.
        forcing_type: str
            Forcing type.
    """
    netcdf_filepath = Path(base_path,
                           f"RAVEN/data/MeteoSwiss_gridded_products/{forcing_type}/out/{forcing_type}_198101010000_202012310000_{ctm_ch_id}_clipped.nc")
    infile = Dataset(filename=netcdf_filepath, mode="r+", format="NETCDF4")
    rio_dem = rxr.open_rasterio(
        Path(base_path, f"RAVEN/data/DEM/out/dem_gridcells_{ctm_ch_id}.tif"),
        masked=True).squeeze()
    rio_pix = rxr.open_rasterio("/media/mainman/Work/RAVEN/data/DEM/shape/ASCII_GRID_1part/dhm25_grid_raster_2056.tif",
                                masked=True).squeeze()
    try:
        ele = infile.createVariable('elevation', np.float32, fill_value=-999.99, dimensions=('N', 'E'))
        ele.units = 'meters'
        ele.long_name = 'elevation above sea level'
        ele.setncatts({'grid_name': "ch01r.swiss.lv95",
                       'version': "v1.2",
                       'prod_date': "2022-10-01",
                       'coordinates': "lon lat",
                       'grid_mapping': u"swiss_lv95_coordinates"})
    except RuntimeError:
        pass
    np_dem = rio_dem.to_numpy()
    np_pix = rio_pix.to_numpy()
    np_dem[22][69] = np_pix
    infile['elevation'][:] = np.flipud(np_dem)
    # Close the dataset and write the file
    infile.close()


def dhm_transform():
    """Reprojects and slices dem.

    """
    dem_old = rxr.open_rasterio("/media/mainman/Work/RAVEN/data/DEM/shape/ASCII_GRID_1part/dhm25_grid_raster.asc",
                                masked=True).squeeze()
    dem_new = dem_old.rio.reproject(dst_crs='EPSG:2056', dst_resolution=(1000.0, 1000.0))
    geometries = [
        {
            'type': 'Polygon',
            'coordinates': [[

                [2740000, 1144000],
                [2740000, 1145000],
                [2741000, 1145000],
                [2741000, 1144000],
                [2740000, 1144000]
            ]]
        }
    ]
    slice_dict = {
        "x": slice("2740000", "2741000"),
        "y": slice("1144000", "1145000")
    }
    dem_slice = dem_new.rio.slice_xy(minx=2739000, miny=1143000, maxx=2742000, maxy=1146000)
    dem_slice.rio.to_raster("/media/mainman/Work/RAVEN/data/DEM/shape/ASCII_GRID_1part/slice_out.tif")


def elevation_band_slope_aspect(data_dir: Path, ctm_ch_id: str, band_list, hru_id, hru_info_df, grid: gpd.GeoDataFrame,
                                grid_weights_df):
    """Calculates slope and aspect for each elevation band and, if not already existent, creates corresponding tif files.

    Args:
        data_dir: Path
            Data directory.
        ctm_ch_id: str
            Catchment extent
        band_list: list
            Elevation band names.
        hru_id: list
            HRU ids.
        hru_info_df: pd.DataFrame
            DataFrame with info on all hrus
        grid: gpd.GeoDataFrame
            Grid from netCDF file
        grid_weights_df: pd.DataFrame

    Returns:
        slope_aspect_info_df: pd.DataFrame
            Slope and aspect for each elevation band.
        area_total: float
            Total catchment area.
        info: pd.DataFrame
    """
    if grid.empty:
        grid = gpd.read_file(Path(data_dir,
                                  f"MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/grid_weights_{ctm_ch_id}.shp"))
    slope_aspect_info_df = pd.DataFrame(columns=['cell_id'])
    for (bd, hru) in zip(band_list, hru_id):
        if bd != "gla":
            dem = rxr.open_rasterio(filename=Path(data_dir,
                                                  f"DEM/hbv/non_glacier/dem_{ctm_ch_id}_{bd}_{str(int(bd) + 99)}.tif"),
                                    masked=True, chunks=True).squeeze()
            dem_df = dem.rename('alti').to_dataframe().reset_index().dropna()
            dem_gdf = gpd.GeoDataFrame(dem_df, crs='epsg:2056',
                                       geometry=gpd.points_from_xy(dem_df.reset_index().x, dem_df.reset_index().y))
            ov = create_overlay(grd=grid, ctm_gdf=dem_gdf)
            cou = ov.value_counts('cell_id')
            # cou = pd.DataFrame()
            cou = cou.to_frame().reset_index()
            dem_aspect: rxr.raster_array
            filename_aspects = Path(data_dir, f"DEM/hbv/aspects/dem_{ctm_ch_id}_{bd}_{str(int(bd) + 99)}_aspects.tif")
            filename_slope = Path(data_dir, f"DEM/hbv/slopes/dem_{ctm_ch_id}_{bd}_{str(int(bd) + 99)}_slopes.tif")
            try:
                dem_aspect = rxr.open_rasterio(filename=filename_aspects, masked=True, chunks=True).squeeze()
                pass
            except rasterio.errors.RasterioIOError:
                filename_dem = Path(data_dir, f"DEM/hbv/non_glacier/dem_{ctm_ch_id}_{bd}_{str(int(bd) + 99)}.tif")
                rcode = subprocess.run(args=['raven_tools/gdal_slope_aspect.sh', str(filename_dem)])
                try:
                    dem_aspect = rxr.open_rasterio(
                        filename=filename_aspects, masked=True, chunks=True).squeeze()
                except:
                    break
            try:
                dem_slope = rxr.open_rasterio(filename=filename_slope, masked=True, chunks=True).squeeze()
            except:
                filename_dem = Path(data_dir, f"DEM/hbv/non_glacier/dem_{ctm_ch_id}_{bd}_{str(int(bd) + 99)}.tif")
                rcode = subprocess.call(['raven_tools/gdal_slope_aspect.sh', str(filename_dem)])
                try:
                    dem_slope = rxr.open_rasterio(
                        filename=filename_aspects, masked=True, chunks=True).squeeze()
                except:
                    break
            cou['Ctm'], cou['Band'], cou['hru_id'] = ctm_ch_id, bd, hru
            cou['aspect_mean'], cou['slope_mean'], cou['alti_mean'] = (dem_mean_rasterio(dem_aspect),
                                                                       dem_mean_rasterio(dem_slope),
                                                                       dem_mean_rasterio(dem))
            cou['centroid_lat'], cou['centroid_lon'] = weighted_centroid_raster_points(raster_points_gdf=dem_gdf)
            slope_aspect_info_df = pd.concat([slope_aspect_info_df, cou], ignore_index=True)
            # res = ov.value_counts('cell_id')

    slope_aspect_info_df.rename(columns={'count': 'area'}, inplace=True)
    slope_aspect_info_df.area = slope_aspect_info_df.area.astype(int)
    slope_aspect_info_df['hru_id'] = slope_aspect_info_df['hru_id'].astype(int)
    area_non_gla = float(hru_info_df[hru_info_df['Ctm'] == ctm_ch_id]['NonGlaArea'].iloc[0])
    area_gla = float(hru_info_df[hru_info_df['Ctm'] == ctm_ch_id]['GlaArea'].iloc[0])
    if math.isnan(area_gla):
        area_gla = 0
    area_total = area_non_gla + area_gla
    # slope_aspect_info_df['ratio'] = slope_aspect_info_df.area / slope_aspect_info_df.area.sum()
    # slope_aspect_info_df['grid_weight'] = slope_aspect_info_df['ratio'] * (area_non_gla / area_total)
    # tes = slope_aspect_info_df.value_counts('cell_id')
    slope_aspect_info_df.sort_values(['hru_id', 'cell_id'])
    info = pd.DataFrame()
    for id in hru_id[:-1]:
        info = pd.concat(
            [info,
             slope_aspect_info_df[slope_aspect_info_df['hru_id'] == id].iloc[0].to_frame().transpose().set_index(
                 'hru_id')])

    info.to_csv(Path(data_dir, "Catchment", f"slope_aspect_{ctm_ch_id}.txt"))
    return slope_aspect_info_df, area_total, info


def elevation_band_ratio(data_dir: Path, ctm_ch_id: str, band_list, hru_id, hru_info_df, grid: gpd.GeoDataFrame):
    """Returns ratio for each elevation band of the total catchment area.

    Args:
        data_dir: Path
            Data directory
        ctm_ch_id: str
            Catchment id
        band_list: list
            Elevation band limits
        hru_id: list
            HRU ids
        hru_info_df: pd.DataFrame
            Info on each HRU
        grid: gpd.GeoDataFrame
            Grid from netCDF file
    """
    if grid.empty:
        grid = gpd.read_file(Path(data_dir,
                                  f"MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/grid_weights_{ctm_ch_id}.shp"))
    for (bd, hru) in zip(band_list, hru_id):
        pass
    elevation_band_ratio_df = pd.DataFrame(columns=['cell_id'])
    elevation_band_ratio_df.rename(columns={'count': 'area'}, inplace=True)
    elevation_band_ratio_df.area = elevation_band_ratio_df.area.astype(int)
    elevation_band_ratio_df['hru_id'] = elevation_band_ratio_df['hru_id'].astype(int)
    area_non_gla = float(hru_info_df[hru_info_df['Ctm'] == ctm_ch_id]['NonGlaArea'].iloc[0])
    area_gla = float(hru_info_df[hru_info_df['Ctm'] == ctm_ch_id]['GlaArea'].iloc[0])
    if math.isnan(area_gla):
        area_gla = 0
    area_total = area_non_gla + area_gla
    elevation_band_ratio_df['ratio'] = elevation_band_ratio_df.area / elevation_band_ratio_df.area.sum()
    elevation_band_ratio_df['grid_weight'] = elevation_band_ratio_df['ratio'] * (area_non_gla / area_total)


def write_hru_info(data_dir, catchment_ch_id, area_total, res):
    """Writes HRU info for each elevation band into a file to be used in HBV .rvh file.

    Args:
        data_dir: Path
            Data directory
        catchment_ch_id: str
            Catchment id
        area_total: float
            Total area
        res: pd.DataFrame
            HRU info for each elevation band.
    """
    with open(Path(data_dir, f"Catchment/HBV/hrus_{catchment_ch_id}.txt"), "w") as f:
        for hruid in res['hru_id'].unique().tolist():
            hru_ratio = res.grid_weight.loc[res.hru_id == hruid].sum()
            hru_area = hru_ratio * area_total
            hru_alti = res.alti_mean.loc[res.hru_id == int(hruid)].unique()[0]
            hru_lat = res.centroid_lat.loc[res.hru_id == int(hruid)].unique()[0]
            hru_lon = res.centroid_lon.loc[res.hru_id == int(hruid)].unique()[0]
            hru_asp = res.aspect_mean.loc[res.hru_id == int(hruid)].unique()[0]
            hru_slo = res.slope_mean.loc[res.hru_id == int(hruid)].unique()[0]
            f.write(
                f"            {hruid}, {hru_area}, {hru_alti},{hru_lat}, {hru_lon}, 1, LU_ALL, VEG_ALL, DEFAULT_P, [NONE], [NONE], {hru_slo}, {hru_asp}\n")


def write_hbv_weights_to_file(grid_weights_file_path, grid_weights_dgdf):
    """Writes grid weights file for HBV with elevation bands.

    Args:
        grid_weights_file_path: Path
            Output file path.
        grid_weights_dgdf: dgpd.GeoDataFrame
            Dask GeoDataFrame with grid weight data.
    """
    hru_id_list = range(1, grid_weights_dgdf.shape[1])
    col_dict = {hru_id: col_name for hru_id, col_name in zip(hru_id_list, (col_name for col_name in
                                                                           grid_weights_dgdf.columns.to_list() if
                                                                           col_name != "polygons"))}
    # gd = dd.from_pandas(pd.DataFrame(), npartitions=8)
    for col_name in col_dict.values():
        grid_weights_dgdf[col_name] = grid_weights_dgdf[col_name] / grid_weights_dgdf[col_name].sum()
    g_list = [grid_weights_dgdf.loc[:, col_name].dropna().rename(hru_id).to_frame() for (hru_id, col_name) in
              col_dict.items()]
    # g_list = [grid_weights_dgdf.loc[:, col_name].rename(hru_id).to_frame() for (col_name, hru_id) in
    #           zip(col_dict, hru_id_list)]
    weight_dd = dd.concat(g_list, keys=col_dict.keys()).compute()
    # for col_num in hru_id_list:
    #     # print(f"{col_num}".join((weight_dd.loc[:, col_num].dropna().to_string(header=False).split(sep="\n"))))
    #     print(
    #         *[f"{col_num}   {it}\n" for it in
    #           weight_dd.loc[:, col_num].dropna().to_string(header=False).split(sep="\n")])
    # g = [weight_dd.loc[:, col_num].dropna().to_string(header=False).split(sep="\n") for col_num in col_dict]
    # j = ["1   ".join(line) for line in g[0]]
    with (open(grid_weights_file_path, "w") as ff):
        ff.write(':GridWeights                     \n')
        ff.write('   #                                \n')
        ff.write('   # [# HRUs]                       \n')
        ff.write(
            f'   :NumberHRUs       {len(col_dict)}            \n')  # Currently for GR4J, there's 1 HRU
        ff.write(f"   :NumberGridCells  {len(grid_weights_dgdf.index.unique())}            \n")
        ff.write('   #                                \n')
        ff.write('   # [HRU ID] [Cell #] [w_kl]       \n')
        ff.writelines(
            (f"{hru_id}   {it}\n" for (hru_id) in col_dict.keys()
             for it in weight_dd.loc[:, hru_id].dropna().to_string(header=False).split(sep="\n")))
        ff.write(':EndGridWeights \n')

    # g = [weight_dd.loc[:, col_num].dropna().to_string(header=False).split(sep='\n') for col_num, hru_id in
    #      col_dict.items()]
    # len(g[0])
    # with (open(grid_weights_file_path, "w") as f):
    #     f.writelines(
    #         (f"{hru_id}   {}\n" for hru_id,  in col_dict.items()))


if __name__ == '__main__':
    base_path = Path("/media/mainman/Work/")
    weights = pd.read_csv(
        "/media/mainman/Work/RAVEN/models/CH-0053/HBV/model/data_obs/RhiresD_v2.0_swiss.lv95/out/grid_weights_CH-0053.txt",
        sep='\s+', skiprows=7, skipfooter=1, names=['hru_id', 'cell_id', 'weight'])
    for i in weights['hru_id'].unique():
        print(i)
        print(weights[weights['hru_id'] == i].sum())
    # for forcing_type in raven_tools.config.variables.forcings_dirs:
    #     for ctm in raven_tools.config.variables.catchments.keys():
    #         add_elevation_to_netcdf(base_path=base_path, ctm_ch_id=ctm)

    # for forcing_type in raven_tools.config.variables.forcings_dirs:
    #     add_elevation_to_netcdf(base_path=base_path, ctm_ch_id="CH-0053", forcing_type=forcing_type)

    # # # catchments_by_id = [key for key in raven_tools.config.variables.catchments]
    # # # r = raster_to_polygon(Path("/home/sirian/Applications/Hydrology/RAVEN/data/DEM/hbv/aspects/dem_CH-0105_1400_1499_aspects.tif"))
    # #
    # catchments_by_id = ["CH-0053"]
    # base_path_prefix = "/media/mainman/Work/"
    # base_path = base_path_prefix + "RAVEN/data/DEM/hbv"
    # base_path_path = Path(base_path)
    # for ctm in catchments_by_id:
    #     # area_ratios, ratio_dict = elevation_bands(
    #     #     filepath_non_glacier=Path(base_path_path.parent, f"dem_clipped/non_glacier/dem_{ctm}.tif"),
    #     #     filepath_glacier=Path(base_path_path.parent,
    #     #                           f"dem_clipped/glacier/dem_{ctm}_glacier.tif"),
    #     #     ctm_id=ctm, save_to_tif=True, base_path_prefix=base_path_prefix)
    #     # # area_ratios = {}
    #     # with open(f"{base_path_prefix}RAVEN/data/DEM/hbv/non_glacier/area_ratios_{ctm}.txt", "r") as f:
    #     #     for line in f:
    #     #         s = line.strip().split(":")
    #     #         area_ratios[s[0]] = float(s[1])
    #     # elevation_band_areas = area_from_ratio_dem_props(area_ratios=area_ratios, ctm_id=ctm,
    #     #                                                  base_path_prefix=base_path_prefix)
    #     # dict_to_txt(dict=elevation_band_areas, ctm=ctm, base_path_prefix=base_path_prefix)
    #     res = pd.DataFrame(columns=['cell_id'])
    #     band_list = create_elevation_band_tif_list_hbv(base_path, 'non_glacier')
    #     hru_id = list(range(2, len(band_list) + 2))
    #     for (bd, hru) in zip(band_list, hru_id):
    #         grid = gpd.read_file(
    #             f"{base_path_prefix}RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/grid_weights_{ctm}.shp")
    #         dem = rxr.open_rasterio(
    #             f"{base_path_prefix}RAVEN/data/DEM/hbv/non_glacier/dem_{ctm}_{bd}_{str(int(bd) + 99)}.tif",
    #             masked=True).squeeze()
    #         dem_df = dem.to_dataframe(name='alti').dropna()
    #         dem_gdf = gpd.GeoDataFrame(dem_df, crs='epsg:2056',
    #                                    geometry=gpd.points_from_xy(dem_df.reset_index().x, dem_df.reset_index().y))
    #         dem_gdf.reset_index(inplace=True)
    #         centroid = dem_gdf.dissolve().centroid
    #         ov = create_overlay(grd=grid, ctm_gdf=dem_gdf)
    #         cou = ov.value_counts('cell_id')
    #         cou = cou.to_frame().reset_index()
    #         dem_aspect = rxr.open_rasterio(
    #             filename=f"{base_path_prefix}RAVEN/data/DEM/hbv/aspects/dem_{ctm}_{bd}_{str(int(bd) + 99)}_aspects.tif",
    #             masked=True).squeeze()
    #         dem_slope = rxr.open_rasterio(
    #             filename=f"{base_path_prefix}RAVEN/data/DEM/hbv/slopes/dem_{ctm}_{bd}_{str(int(bd) + 99)}_slopes.tif",
    #             masked=True).squeeze()
    #         cou['Ctm'] = ctm
    #         cou['Band'] = bd
    #         cou['hru_id'] = hru
    #         cou['aspect_mean'] = dem_mean_rasterio(dem_aspect)
    #         cou['slope_mean'] = dem_mean_rasterio(dem_slope)
    #         cou['alti_mean'] = dem_mean_rasterio(dem)
    #         cou['centroid_lat'], cou['centroid_lon'] = weighted_centroid_raster_points(raster_points_gdf=dem_gdf)
    #         lat, lon = crs_old_to_new(lat_old=dem_gdf.to_crs(epsg=32632).y.mean(),
    #                                   lon_old=dem_gdf.to_crs(epsg=32632).y.mean(), epsg_old=32632, epsg_new=4326)
    #         res = pd.concat([res, cou], ignore_index=True)
    #
    #         # res = ov.value_counts('cell_id')
    #     res.rename(columns={0: 'area'}, inplace=True)
    #     res.area = res.area.astype(int)
    #     res['hru_id'] = res['hru_id'].astype(int)
    #     hru_info = pd.read_csv(f"{base_path_prefix}RAVEN/data/Catchment/hru_info.csv", sep=",")
    #     area_non_gla = float(hru_info[hru_info['Ctm'] == ctm]['NonGlaArea'])
    #     area_gla = float(hru_info[hru_info['Ctm'] == ctm]['GlaArea'])
    #     area_total = area_non_gla + area_gla
    #     res['ratio'] = res.area / res.area.sum()
    #     res['grid_weight'] = res['ratio'] * (area_non_gla / area_total)
    #     tes = res.value_counts('cell_id')
    #     res.sort_values(['hru_id', 'cell_id'])
    #
    #     totar: float = 0
    #     with open(f"/tmp/hrus_{ctm}.txt", "w") as f:
    #         for hruid in res['hru_id'].unique().tolist():
    #             hru_ratio = res.grid_weight.loc[res.hru_id == hruid].sum()
    #             hru_area = hru_ratio * area_total
    #             hru_alti = res.alti_mean.loc[res.hru_id == int(hruid)].unique()[0]
    #             hru_lat = res.centroid_lat.loc[res.hru_id == int(hruid)].unique()[0]
    #             hru_lon = res.centroid_lon.loc[res.hru_id == int(hruid)].unique()[0]
    #             hru_asp = res.aspect_mean.loc[res.hru_id == int(hruid)].unique()[0]
    #             hru_slo = res.slope_mean.loc[res.hru_id == int(hruid)].unique()[0]
    #             f.write(
    #                 f"            {hruid}, {hru_area}, {hru_alti},{hru_lat}, {hru_lon}, 1, LU_ALL, VEG_ALL, DEFAULT_P, [NONE], [NONE], {hru_slo}, {hru_asp}\n")
    #
    #     weight_sum = 0
    #     with open(f"/tmp/weights_{ctm}.txt", "w") as f:
    #         for hruid in res['hru_id'].unique():
    #             for cellnum in res[res['hru_id'] == hruid]['cell_id'].sort_values():
    #                 weight = \
    #                     res.grid_weight.loc[(res.hru_id == hruid) & (res.cell_id.astype(int) == int(cellnum))].values[
    #                         0] / \
    #                     (res.grid_weight.loc[(res.hru_id == hruid)].sum())
    #                 weight_sum = weight_sum + weight
    #                 print(weight_sum)
    #                 f.write(
    #                     f"   {hruid}   {cellnum}   {weight}\n")

    # res = pd.concat([res, pd.DataFrame(ov)], ignore_index=True)
    # type = "aspects"
    # band_lower_list = create_elevation_band_tif_list(base_path=base_path, type=type)
    #
    # mean_asp_dict = {}
    # mean_slo_dict = {}
    # mean_alti_dict = {}
    # r = []
    # for b in band_lower_list:
    #     r.append(raster_to_polygon(Path(base_path, type, f"dem_{ctm}_{int(b)}_{int(b) + 99}_{type}.tif")))
    # r_concat = gpd.GeoDataFrame(pd.concat(r))
    # r_concat.drop(columns=["id", "ezgnr", "ar_1903", "watr_nm", "place", "geometry"]).set_index(
    #     "lower_bound").to_csv(Path(base_path, type, f"dem_{ctm}_{type}.txt"))
    # with open(Path(base_path, type, f"dem_{ctm}_{type}.txt"), "w") as f:
    #     f.writelines([f"{str(i['area'][0])}\n" for i in r])
    # mean_slo = dem_mean(Path(base_path, f"/slopes/dem_{ctm}_1400_1499_slopes.tif"))
    # mean_hei = dem_mean(Path(base_path, f"/non_glacier/dem_{ctm}_1400_1499.tif"))


def calc_hru_weight(basic_grid, hru_extent_gdf, hru_short_name: str = "hru_short_name") -> gpd.GeoDataFrame:
    hru_res_union = create_overlay(grd=basic_grid, ctm_gdf=hru_extent_gdf)
    hru_rel_area = calc_relative_area(hru_res_union, hru_short_name=hru_short_name)
    new_grid: gpd.GeoDataFrame = basic_grid.set_index('cell_id').join(other=hru_rel_area.set_index('cell_id'),
                                                                      rsuffix=f'_{hru_short_name}')
    return new_grid


def calc_hru_weight_to_grid(basic_grid, hru):
    pass
