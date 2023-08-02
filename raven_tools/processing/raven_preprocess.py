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
from typing import Union

import geopandas as gpd
import netCDF4
import numpy as np
import pandas as pd
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


def create_bbox_geometry(extent_shape_path: Path):
    """Create a bounding box Polygon for an input shape file.

    Args:
        extent_shape_path : Path
            Path to the shape file for which to create a bounding box.

    Returns:
        bbox_poly : Polygon
            A shapely polygon that defines the bounding box.
        ext_gdf : GeoDataFrame
            The GeoDataFrame with the data from the extent shape file.

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
        create_bbox_shp: bool
            True if create_bbox() should write to a shape file.
        extent_shape_file_path : Path
            Path to the shape file for which to create a bounding box (in EPSG2056).
        bb_file_path : Path
            Full path to the bounding box shape file to be written (in EPSG2056).

    Returns:
        bbox_gdf : GeoDataFrame
            Bounding box as a GeoDataFrame (in EPSG2056).
        ext_gdf : GeoDataFrame
            The GeoDataFrame with the data from the extent shape file (in EPSG2056)
        bb_file_path: str
            (optional) Path to the bounding box shape file (in EPSG2056)
    """

    # Create the bounding box Polygon
    bbox_poly, ext_gdf = create_bbox_geometry(extent_shape_file_path)
    # Create the bounding shape in a GeoDataFrame
    bbox_gdf: GeoDataFrame = gpd.GeoDataFrame(pd.DataFrame(['p1'], columns=['geom']),
                                              crs={'init': 'epsg:2056'},
                                              geometry=[bbox_poly])
    if create_bbox_shp:
        # This writes the bounding box as a shape file
        bbox_gdf.to_file(str(bb_file_path))
        return bbox_gdf, ext_gdf, bb_file_path
    else:
        return bbox_gdf, ext_gdf


def netcdf_to_dataset(netcdf_file_path: Path) -> xr.Dataset:
    """Reads a netCDF file into an xarray dataset

    Args:
        netcdf_file_path : Path
            Path to the netCDF file to clip. It has to be in CH1903+/LV95 (EPSG=2056).

    Returns:
        xds : Dataset
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
        catchment:
        xds_to_write : Dataset
            xarray Dataset to write.
        netcdf_file_path : Path
            netCDF file path to write to.

    """
    try:
        os.mkdir(f"{netcdf_file_path.parent.parent}/out/{catchment}")
    except FileExistsError:
        pass
    # Write the clipped netCDF file
    xds_to_write.to_netcdf(
        f"{netcdf_file_path.parent.parent}/out/{catchment}/{netcdf_file_path.stem}_{catchment}_clipped{netcdf_file_path.suffix}",
        mode="w")


def netcdf_clipper(netcdf_file_path: Path, extent_file_path: Path):
    """

    Args:
        netcdf_file_path:
        extent_file_path:
    """
    logger.debug("Trying to call extent.sh...")
    rcode = subprocess.call(['raven_tools/extent.sh', str(extent_file_path), str(netcdf_file_path)])
    logger.debug(f"extent.sh executed with return code: {rcode}")


def netcdf_pet_hamon(netcdf_file_path: Path, name_pattern: dict[str, str]):
    """

    Args:
        name_pattern : dict[str, str]
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


def pet_temp_monthly_ave(pet_filepath: Path, temp_filepath: Path):
    """Calculates monthly averages for temperature and PET from CSV values

    Args:
        pet_filepath: Path
        temp_filepath: Path

    Returns:
        pet_monthly:
        temp_monthly:
    """
    import pandas as pd
    from pathlib import Path
    pet_filepath = Path(pet_filepath)
    pet_monthly_from_order = pd.read_csv(pet_filepath, sep=";")

    pet_monthly_from_order['time'] = pd.to_datetime(pet_monthly_from_order['time'], format='%Y%m')
    pet_monthly_from_order['month'] = pet_monthly_from_order['time'].dt.month
    pet_monthly = pet_monthly_from_order.groupby(pet_monthly_from_order.time.dt.month)['ets150m0'].mean()
    print(str(pet_monthly_from_order['ets150m0'].mean()))
    print(str(pet_monthly.values.mean()))

    temp_filepath = Path(temp_filepath)
    temp_hourly_from_order = pd.read_csv(temp_filepath, sep=";")
    temp_hourly_from_order['time'] = pd.to_datetime(temp_hourly_from_order['time'], format='%Y%m%d%H')
    temp_monthly = temp_hourly_from_order.groupby(temp_hourly_from_order.time.dt.month)['tre200h0'].mean()
    print(temp_monthly.values)
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


def nc_merge(start_year: int, end_year: int, forcing_dir: Path, forcing_prefix: str, catchment: str = "Broye",
             legacy_mode: bool = False):
    """

    Args:
        legacy_mode: bool
            True if you want to use legacy mode
        catchment: str
        forcing_prefix: str
        start_year : int
            Start year
        end_year : int
            End Year
        forcing_dir : Path
            Root directory where forcing files are located

    """
    if legacy_mode:
        logger.debug("Trying to call nc_combine.sh...")
        rcode = subprocess.call(
            ['raven_tools/nc_combine.sh', str(start_year), str(end_year), str(forcing_dir), catchment])
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
        start_year :
            Used to slice data so as to save computation time
        forcing_name : str
            Used to select column in netCDF data
        out_path : Path
            If export_shp=True, write shape file to this Path
        netcdf_filepath : Path
            Path to netCDF file whose grid should be used (in EPSG2056)
        bounding_box_filename : Path
            Path to bounding box shape file (in EPSG2056)
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
    grid["polygons"] = polygons
    grid["cell_id"] = cell_id

    if export_shp:
        # Export the grid to a shape file
        grid.to_file(str(str(out_path) + ".shp"))
    logger.debug(f"grid weights shp file written: {out_path}")
    return grid


def create_overlay(grd: GeoDataFrame, ctm_gdf: GeoDataFrame):
    """Overlays a GeoDataFrame over another to create overlay Polygons

    Overlays two GeoDataFrame over each other and returns two new GeoDataFrames, one for the mode 'intersection',
    the second for the mode 'difference'

    Args:
        ctm_gdf: GeoDataFrame
            Catchment GDF
        grd : GeoDataFrame
            Grid as given by the netCDF file in EPSG2056

    Returns:
        res_u : GeoDataFrame
            Grid cells within the catchment area in EPSG2056

    """
    ctm_gdf.set_crs(epsg="2056")
    grd.set_crs(epsg='2056')
    ctm_gdf.to_crs(32632)
    grd.to_crs(32632)
    res_u: GeoDataFrame = ctm_gdf.overlay(grd, how='intersection')
    # res_u.set_index("cell_id", inplace=True)
    res_u.to_crs(2056)
    return res_u


def calc_relative_area(gdf: GeoDataFrame, hru_short_name: str) -> GeoDataFrame:
    """Calculates the relative area of each polygon in a GeoDataFrame.

    Calculates the relative area of each polygon in a GeoDataFrame with EPSG=2056, writes it into a new column and
    returns the GeoDataFrame with EPSG=2056.

    Args:
        hru_short_name: str
        gdf : GeoDataFrame
            GeoDataFrame in EPSG=2056

    Returns:
        gdf : GeoDataFrame
            GeoDataFrame with relative areas of each polygon

    """
    area_column_name: str = f"{hru_short_name}_area_rel"

    # Set the CRS to WGS84 to preserve areas
    gdf = gdf.to_crs(32632)

    # For each intersected feature, compute the area and write into a new field
    gdf[area_column_name] = (gdf["geometry"].area * 1000000) / np.sum(gdf.geometry.area * 1000000)
    # Convert the area to float
    area_sum: float = float(gdf["area"].sum())
    # for index, row in gdf.iterrows():
    #     # To minimize numerical errors, multiply and then divide the result by 1000000
    #     gdf.at[index, area_column_name] = (gdf.loc[index]["area"] / area_sum)
    #     gdf.at[index, "index"] = int(index)
    # Re-project back into the LV95 projection and export to a shape file
    gdf: GeoDataFrame = gdf.to_crs(2056)
    return gdf


def write_weights_to_file(grd: dict[GeoDataFrame], grid_dir_path: Path, glacier: bool = False):
    """Write grid weights to Raven compatible file

    Args:
        glacier: bool
            Set True if glacier
        grd : GeoDataFrame
            Grid derived from the netCDF file
        grid_dir_path : Path
            Path to the grid weights text file

    """
    # Write to GridWeights.txt
    # The following section has been adapted from Juliane Mai's derive_grid_weights.py script. It writes to a file
    # that is compatible with RAVEN
    grid_filepath = Path(str(grid_dir_path) + ".txt")
    write_grid_data_to_file(grd=grd, grid_filepath=grid_filepath, glacier=glacier)


def create_grid_data_list(grd: GeoDataFrame, glacier: bool = False):
    """Loops over each grid cell and extracts the grid weights.

    Args:
        glacier : bool
            True if glacier
        grd : GeoDataFrame
            Grid as derived from the netCDF file

    Returns:
        data_to_write : list[list[float]]
            List with the relative areas/grid weights of each cell.

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


def write_grid_data_to_file(grd, grid_filepath: Path):
    """Writes grid weight data to rvt file.

    Args:
        grd:
        grid_filepath: Path
    """
    with open(grid_filepath, 'w') as ff:
        ff.write(':GridWeights                     \n')
        ff.write('   #                                \n')
        ff.write('   # [# HRUs]                       \n')
        total_cells = len(grd['1'])
        if not grd['2'].empty:
            num_hrus = 2
        else:
            num_hrus = 1
        ff.write(f'   :NumberHRUs       {num_hrus}            \n')  # Currently for GR4J, there's 1 HRU
        ff.write(f'   :NumberGridCells  {total_cells}            \n')
        ff.write('   #                                \n')
        ff.write('   # [HRU ID] [Cell #] [w_kl]       \n')
        averages = {}
        df = grd['1'].reset_index()
        for cell_id in df['cell_id'].unique():
            tempdf = df[df['cell_id'] == cell_id]
            average = tempdf['non_gla_area_rel'].mean()
            averages[cell_id] = [average]
            weights = gpd.GeoDataFrame.from_dict(averages, orient='index', columns=['area'])
            weights.fillna(0, inplace=True)
        for index, weight in zip(weights.index, weights['area']):
            ff.write(f"   1   {index}   {weight}\n")
        if not grd['2'].empty:
            del weights
            averages = {}
            df = grd['1'].reset_index()
            for cell_id in df['cell_id'].unique():
                tempdf = df[df['cell_id'] == cell_id]
                average = tempdf['non_gla_area_rel'].mean()
                averages[cell_id] = [average]
                weights = gpd.GeoDataFrame.from_dict(averages, orient='index', columns=['area'])
                weights.fillna(0, inplace=True)
            for index, weight in zip(weights.index, weights['area']):
                ff.write(f"   2   {index}   {weight}\n")
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
        start_date: str
        end_date: str

    Returns:
        subset_dataframe: pd.DataFrame

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
        start_time: str
        df: pd.DataFrame
        out_path:

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


def glacier_extent(ctm_gdf: GeoDataFrame, glacier_gdf: GeoDataFrame):
    """Creates a GDF with the glaciated area in a catchment from GeoDataFrames.

    Args:
        ctm_gdf: GeoDataFrame
        glacier_gdf: GeoDataFrame

    Returns:
        ctm_glaciation: GeoDataFrame
        ctm_non_glaciation: GeoDataFrame

    """
    ctm_glaciation: GeoDataFrame = ctm_gdf.overlay(glacier_gdf, how='intersection')
    ctm_non_glaciation: GeoDataFrame = ctm_gdf.overlay(glacier_gdf, how='difference')
    # ctm_glaciation.set_crs(epsg='2056')
    # ctm_non_glaciation.set_crs(epsg='2056')
    # ctm_glaciation.set_index("cell_id")
    # ctm_glaciation.set_crs(epsg='2056')
    return ctm_glaciation, ctm_non_glaciation


def glacier_extent_from_shp(ctm_shp: Path, glacier_shp: Path):
    """Creates a GDF with the glaciated area in a catchment from shape files.

    Args:
        ctm_shp: Path
        glacier_shp: Path

    Returns:
        ctm_glaciation: GeoDataFrame
        ctm_non_glaciation: GeoDataFrame

    """
    ctm_gdf: GeoDataFrame = gpd.read_file(ctm_shp)
    glacier_gdf: GeoDataFrame = gpd.read_file(glacier_shp)
    ctm_glaciation: GeoDataFrame = ctm_gdf.overlay(glacier_gdf, how='intersection')
    ctm_non_glaciation: GeoDataFrame = ctm_gdf.overlay(glacier_gdf, how='difference')

    return ctm_glaciation, ctm_non_glaciation


def area_from_ratio_dem_props(area_ratios: dict, ctm_id: str):
    props = pd.read_csv("/media/mainman/Work/RAVEN/data/Catchment/hru_info.csv")
    non_gla_area = float(props[props['Ctm'] == ctm_id]['NonGlaArea'])
    gla_area = float(props[props['Ctm'] == ctm_id]['NonGlaArea'])
    total_area = non_gla_area + gla_area
    for lower_bound, ratio in area_ratios.items():
        area_ratios[lower_bound] = ratio * non_gla_area
    return area_ratios


def dict_to_txt(dict, ctm):
    with open(f"/media/mainman/Work/RAVEN/data/DEM/hbv/non_glacier/elevation_band_areas_{ctm}.txt", 'w') as f:
        for key, value in dict.items():
            f.write(f"{key}:{value}\n")


def area_ratio(catchment_filepath: Path, glacier_shape_path: Path):
    """Calculates the relative glaciated area of a catchment

    Args:
        catchment_filepath: Path
        glacier_shape_path: Path

    Returns:
        ctm_total_area: float
        non_glaciation_area:
        glaciation_area (optional):

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


def weighted_centroid_vector(feature_gdf: GeoDataFrame) -> tuple:
    """Calculates weighted centroid of a GDF feature

    Args:
        feature_gdf: The feature whose centroid to calculate.

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
    """Calculates mean height of a DEM

    Args:
        filepath: Path
            Path to DEM

    Returns:
        mean:

    """
    import rioxarray as rxr

    dem_im = rxr.open_rasterio(filepath, masked=True).squeeze()
    dem_mean_raster = dem_im.mean(dim=["x", "y"], skipna=True).to_numpy()
    mean = dem_mean_raster.tolist()
    return mean


def dem_mean_rasterio(dem: Union[Dataset, xr.DataArray, list[Dataset]]):
    dem_mean_raster = dem.mean(dim=["x", "y"], skipna=True).to_numpy()
    mean = dem_mean_raster.tolist()
    return mean


def elevation_bands(filepath_non_glacier: Path, filepath_glacier: Path, ctm_id: str):
    df_non_glacier = load_dem_tif_to_dataframe(filepath_non_glacier)
    df_glacier = load_dem_tif_to_dataframe(filepath_glacier)
    rio_non_glacier = rxr.open_rasterio(filepath_non_glacier, masked=True).squeeze()
    lower = round_down(df_non_glacier.alti.min())
    upper = round_up(df_non_glacier.alti.max())
    number_points_total: int = int(df_non_glacier['alti'].count())

    rio_elevation_band_dict = extract_elevation_band_from_rio_dem(dem=rio_non_glacier, lower=lower, upper=upper,
                                                                  save_to_tif=False, ctm_id=ctm_id)

    number_points = {}
    for i in rio_elevation_band_dict.keys():
        number_points[i] = (rio_elevation_band_dict[i].count().data.min() * 1000) / (number_points_total * 1000)

    # df_elevation_band_list = load_rio_dataset_list_to_dataframe(rio_elevation_band_list)

    rio_elevation_band_dict_to_txt(rio_elevation_band_dict, ctm_id, df_non_glacier, df_glacier)

    dem_glacier = rxr.open_rasterio(filepath_glacier, masked=True).squeeze()

    dem_im = rxr.open_rasterio(filepath_non_glacier, masked=True).squeeze()
    dem_new = dem_im.where(dem_im.values > 2200).where(dem_im.values < 2300)
    # dem_new.rio.to_raster("/home/sirian/Downloads/out.tif")
    # dem_gdf = gpd.read_file("/home/sirian/Downloads/out.tif")
    dem_new.name = "data"
    dem_im.name = "data"
    df_im = dem_im.to_dataframe().reset_index()
    df_new = dem_new.to_dataframe().reset_index()
    geometry = gpd.points_from_xy(df_new.x, df_new.y)
    gdf = gpd.GeoDataFrame(df_new, crs=dem_new.rio.crs, geometry=geometry)
    df_nan = df_new[df_new['data'].isna()]
    with open(f"/media/mainman/Work/RAVEN/data/DEM/hbv/non_glacier/area_ratios_{ctm_id}.txt", 'w') as f:
        for key, value in number_points.items():
            f.write(f"{key}:{value}\n")

    return number_points


def raster_to_polygon(raster_file_path: Path):
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
    ratio_dict = {}
    it = rio_elevation_band_dict.items()
    for key, elevation_band in it:
        df_band = load_rio_dataset_to_dataframe(elevation_band)
        ratio_dict[key] = ((df_non_glacier.count().y - df_band[df_band['alti'].isna()].count().y) / \
                           ((df_non_glacier.count().y - df_non_glacier[df_non_glacier['alti'].isna()].count().y) + \
                            (df_glacier.count().y - df_glacier[df_glacier['alti'].isna()].count().y)))
    with open(f"/media/mainman/Work/RAVEN/data/DEM/hbv/bands_ratio_{ctm_id}.txt", mode="w") as f:
        f.write(json.dumps(ratio_dict))


def round_up(x):
    return int(math.ceil(x / 100.0)) * 100


def round_down(x):
    return int(math.floor(x / 100.0)) * 100


def load_dem_tif_to_dataframe(dem_tif_filepath: Path):
    rio_dem = rxr.open_rasterio(dem_tif_filepath, masked=True).squeeze()
    rio_dem.name = "alti"
    return rio_dem.to_dataframe().reset_index()


def load_rio_dataset_list_to_dataframe(rio_dataset_list):
    df_list = []
    for idx, i in enumerate(rio_dataset_list):
        i.name = "alti"
        df_list.append(i.to_dataframe().reset_index())
    return df_list


def load_rio_dataset_to_dataframe(rio_dataset):
    rio_dataset.name = "alti"
    return rio_dataset.to_dataframe().reset_index()


def extract_elevation_band_from_rio_dem(dem: xr.DataArray, lower: int, upper: int, save_to_tif: False, ctm_id) -> dict:
    elevation_band_dict = {}
    for r in range(lower, upper, 100):
        elevation_band = dem.where(dem.values >= r).where(dem.values < r + 100)
        elevation_band_dict[r] = elevation_band
        if save_to_tif:
            elevation_band.rio.to_raster(f"/media/mainman/Work/RAVEN/data/DEM/hbv/dem_{ctm_id}_{r}_{r + 99}.tif")
    return elevation_band_dict


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


def create_elevation_band_tif_list(base_path: str, type: str):
    asp_list = glob.glob(base_path + f'/{type}/dem_{ctm}_*_{type}.tif')
    band_lower_list = []
    for a in asp_list:
        rs = re.search(r'(?:_)(\d{4}?)(?:_)', a).group()
        band_lower_list.append(re.search(r'(\d{4})', re.search(r'(?:_)(\d{4}?)(?:_)', a).group()).group())
    band_lower_list.sort()
    return band_lower_list


def create_elevation_band_tif_list_hbv(base_path: str, type: str):
    asp_list = glob.glob(base_path + f'/{type}/dem_{ctm}_*.tif')
    band_lower_list = []
    for a in asp_list:
        rs = re.search(r'(?:_)(\d{4}?)(?:_)', a).group()
        band_lower_list.append(re.search(r'(\d{4})', re.search(r'(?:_)(\d{4}?)(?:_)', a).group()).group())
    band_lower_list.sort()
    return band_lower_list


def get_lower_band_limit_from_filepath(filepath):
    lower_band_limit = re.search(r'(\d{4})', re.search(r'(?:_)(\d{4}?)(?:_)', filepath).group()).group()
    return lower_band_limit


if __name__ == '__main__':
    # catchments_by_id = [key for key in raven_tools.config.variables.catchments]
    # r = raster_to_polygon(Path("/home/sirian/Applications/Hydrology/RAVEN/data/DEM/hbv/aspects/dem_CH-0105_1400_1499_aspects.tif"))
    catchments_by_id = ["CH-0105"]
    base_path = "/media/mainman/Work/RAVEN/data/DEM/hbv"
    base_path_path = Path(base_path)
    for ctm in catchments_by_id:
        # area_ratios = elevation_bands(
        #     filepath_non_glacier=Path(base_path_path.parent, f"dem_clipped/non_glacier/dem_{ctm}.tif"),
        #     filepath_glacier=Path(base_path_path.parent,
        #                           f"dem_clipped/glacier/dem_{ctm}_glacier.tif"),
        #     ctm_id=ctm)
        area_ratios = {}
        with open(f"/media/mainman/Work/RAVEN/data/DEM/hbv/non_glacier/area_ratios_{ctm}.txt", "r") as f:
            for line in f:
                s = line.strip().split(":")
                area_ratios[s[0]] = float(s[1])
        elevation_band_areas = area_from_ratio_dem_props(area_ratios=area_ratios, ctm_id=ctm)
        dict_to_txt(dict=elevation_band_areas, ctm=ctm)
        res = pd.DataFrame(columns=['cell_id'])
        band_list = create_elevation_band_tif_list_hbv(base_path, 'non_glacier')
        hru_id = list(range(2, len(band_list) + 1))
        for (bd, hru) in zip(band_list, hru_id):
            grid = gpd.read_file(
                f"/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/grid_weights_{ctm}.shp")
            dem = rxr.open_rasterio(
                f"/media/mainman/Work/RAVEN/data/DEM/hbv/non_glacier/dem_{ctm}_{bd}_{str(int(bd) + 99)}.tif",
                masked=True).squeeze()
            dem_df = dem.to_dataframe(name='alti').dropna()
            dem_gdf = gpd.GeoDataFrame(dem_df, crs='epsg:2056',
                                       geometry=gpd.points_from_xy(dem_df.reset_index().x, dem_df.reset_index().y))
            dem_gdf.reset_index(inplace=True)
            centroid = dem_gdf.dissolve().centroid
            ov = create_overlay(grd=grid, ctm_gdf=dem_gdf)
            cou = ov.value_counts('cell_id')
            cou = cou.to_frame().reset_index()
            dem_aspect = rxr.open_rasterio(
                filename=f"/media/mainman/Work/RAVEN/data/DEM/hbv/aspects/dem_{ctm}_{bd}_{str(int(bd) + 99)}_aspects.tif",
                masked=True).squeeze()
            dem_slope = rxr.open_rasterio(
                filename=f"/media/mainman/Work/RAVEN/data/DEM/hbv/slopes/dem_{ctm}_{bd}_{str(int(bd) + 99)}_slopes.tif",
                masked=True).squeeze()
            cou['Ctm'] = ctm
            cou['Band'] = bd
            cou['hru_id'] = hru
            cou['aspect_mean'] = dem_mean_rasterio(dem_aspect)
            cou['slope_mean'] = dem_mean_rasterio(dem_slope)
            cou['alti_mean'] = dem_mean_rasterio(dem)
            cou['centroid_lat'], cou['centroid_lon'] = weighted_centroid_raster_points(raster_points_gdf=dem_gdf)
            lat, lon = crs_old_to_new(lat_old=dem_gdf.to_crs(epsg=32632).y.mean(),
                                      lon_old=dem_gdf.to_crs(epsg=32632).y.mean(), epsg_old=32632, epsg_new=4326)
            res = pd.concat([res, cou], ignore_index=True)

            # res = ov.value_counts('cell_id')
        res.rename(columns={0: 'area'}, inplace=True)
        res.area = res.area.astype(int)
        res['hru_id'] = res['hru_id'].astype(int)
        hru_info = pd.read_csv("/media/mainman/Work/RAVEN/data/Catchment/hru_info.csv", sep=",")
        area_non_gla = float(hru_info[hru_info['Ctm'] == ctm]['NonGlaArea'])
        area_gla = float(hru_info[hru_info['Ctm'] == ctm]['GlaArea'])
        area_total = area_non_gla + area_gla
        res['ratio'] = res.area / res.area.sum()
        res['grid_weight'] = res['ratio'] * (area_non_gla / area_total)
        tes = res.value_counts('cell_id')
        res.sort_values(['hru_id', 'cell_id'])

        totar: float = 0
        with open("/tmp/hrus.txt", "w") as f:
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

        with open("/tmp/weights.txt", "w") as f:
            for hruid in res['hru_id'].unique():
                for cellnum in res[res['hru_id'] == hruid]['cell_id'].sort_values():
                    weight = \
                        res.grid_weight.loc[(res.hru_id == hruid) & (res.cell_id.astype(int) == int(cellnum))].values[
                            0] / \
                        (res.grid_weight.loc[(res.hru_id == hruid)].sum())
                    f.write(
                        f"   {hruid}   {cellnum}   {weight}\n")
        res = dem_gdf.overlay(grid, how='intersection')

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
