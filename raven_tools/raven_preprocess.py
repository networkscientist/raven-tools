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
import shutil
import subprocess
from pathlib import Path

import geopandas as gpd
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from geopandas import GeoDataFrame
from netCDF4 import Dataset
from numpy import exp
from shapely.geometry import Polygon, mapping


# from xarray import Dataset


def create_bbox_geometry(extent_shape_path: Path):
    """Create a bounding box Polygon for an input shape file.

    Parameters
    ----------
    extent_shape_path : Path
        Path to the shape file for which to create a bounding box.

    Returns
    -------
    bbox_poly : Polygon
        A shapely polygon that defines the bounding box.
    ext_gdf : GeoDataFrame
        The GeoDataFrame with the data from the extent shape file.

    """

    ext_gdf = gpd.read_file(extent_shape_path)
    se = ext_gdf.geometry.total_bounds  # Get bounds and store in array
    bbox_poly = Polygon([[se[0], se[1]], [se[2], se[1]], [se[2], se[3]], [se[0], se[3]]])
    return bbox_poly, ext_gdf  # Return the Polygon and the GeoDataFrame


def create_bounding_shape(extent_shape_file_path: Path, bb_file_path: Path):
    """Create a bounding box shape file.

    Parameters
    ----------
    extent_shape_file_path : Path
        Path to the shape file for which to create a bounding box.
    bb_file_path : Path
        Full path to the bounding box shape file to be written.
    Returns
    -------
    bbox_gdf : GeoDataFrame
        Bounding box as a GeoDataFrame.
    ext_gdf : GeoDataFrame
        The GeoDataFrame with the data from the extent shape file
    """

    # Create the bounding box Polygon
    bbox_poly, ext_gdf = create_bbox_geometry(extent_shape_file_path)
    # Create the bounding shape in a GeoDataFrame
    bbox_gdf: GeoDataFrame = gpd.GeoDataFrame(pd.DataFrame(['p1'], columns=['geom']),
                                              crs={'init': 'epsg:21781'},
                                              geometry=[bbox_poly])
    # This writes the bounding box as a shape file
    bbox_gdf.to_file(str(bb_file_path))
    return bbox_gdf, ext_gdf


def netcdf_to_dataset(netcdf_file_path: Path):
    """Reads a netCDF file into an xarray dataset

    :param Path netcdf_file_path: Path to the netCDF file to clip
    :return xds: The netCDF data as an xarray dataset
    :rtype xds: xr.Dataset

    """

    # Read in the netCDF file into an xarray Dataset
    xds: xr.Dataset = xr.open_dataset(netcdf_file_path)
    # Define the dimensions
    xds.rio.set_spatial_dims(x_dim="E", y_dim="N", inplace=True)
    # Set the CRS projection EPSG to 2056
    # TODO: Why set CRS to EPSG=2056?
    xds.rio.write_crs("EPSG:2056", inplace=True)
    # delete the attribute 'grid_mapping' to prevent an error
    vars_list = list(xds.data_vars)
    for var in vars_list:
        del xds[var].attrs['grid_mapping']
    return xds


def dataset_to_netcdf(xds_to_write: xr.Dataset, netcdf_file_path: Path):
    """ Writes xarray dataset to netCDF

    :param xr.Dataset xds_to_write xds_to_write: xarray Dataset to write.
    :param Path netcdf_file_path: netCDF file path to write to.

    """

    # Write the clipped netCDF file
    xds_to_write.to_netcdf(
        f"{netcdf_file_path.parent.parent}/out/{netcdf_file_path.stem}_clipped{netcdf_file_path.suffix}", "w")


def netcdf_clipper(netcdf_file_path: Path, bbox_file_path: Path, ext_gdf: GeoDataFrame):
    """Clips a netCDF file according to a bounding box.

    For one netCDF file in a directory, clips it according to a bounding box shape file.

    :param GeoDataFrame ext_gdf: The GeoDataFrame with the data from the extent shape file
    :param Path netcdf_file_path: Path to the netCDF file to clip
    :param Path bbox_file_path: Path to the shape file of the bounding box to be used to clip.
    :return xds_clipped: The clipped Dataset
    :rtype xds_clipped: xr.Dataset

    """

    xds: xr.Dataset = netcdf_to_dataset(netcdf_file_path)
    bbox_gdf = gpd.read_file(bbox_file_path)
    # Clip the xarray Dataset according to the bounding box GeoDataFrame
    xds_clipped = xds.rio.clip(bbox_gdf.geometry.apply(mapping), ext_gdf.crs)
    # Saves the clipped file as shape file
    dataset_to_netcdf(xds_clipped, netcdf_file_path)
    return xds_clipped


def netcdf_clipper_multi(netcdf_dir_path: Path, bbox_file_path: Path, bbox_gdf: GeoDataFrame):
    """ Clips multiple netCDF files in a directory

    :param bbox_file_path: Full Path to the bounding box shape file
    :param netcdf_dir_path: Path to directory with netCDF files to clip.
    :param bbox_gdf: Bounding box GeoDataFrame created with create_bounding_shape()

    """

    for f in glob.glob(f"{netcdf_dir_path}/original_files/*.nc"):
        netcdf_clipper(Path(f), bbox_file_path, bbox_gdf)


def netcdf_pet_hamon(netcdf_file_path: Path, name_pattern: dict[str, str]):
    """

    :param name_pattern:
    :param Path netcdf_file_path: Path to the netCDF file to calculate PET from
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


def nc_merge(start_year: int, end_year: int, forcing_dir: str):
    """

    :param start_year: Start year
    :param end_year: End Year
    :param forcing_dir: Root directory where forcing files are located
    :type end_year: int
    :type start_year: int
    :type forcing_dir: str
    """
    subprocess.call(['raven_tools/nc_combine.sh', str(start_year), str(end_year), forcing_dir])

if __name__ == '__main__':
    print()