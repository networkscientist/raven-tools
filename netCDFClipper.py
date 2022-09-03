# import glob
from pathlib import Path

# import rioxarray
import geopandas as gpd
import pandas as pd
import xarray
from geopandas import GeoDataFrame
from pandas import DataFrame
from shapely.geometry import mapping, Polygon
from xarray import Dataset

# Paths
home_path = Path.home()
netcdf_files_path = f"{home_path}/Applications/Hydrology/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95"  # the original netCDF files are stored here
hydro_shp_path = f"{home_path}/Applications/Hydrology/RAVEN/data/Hydro"  # Hydrology shape files are stored here
out_path = f"{netcdf_files_path}/out"
extent_shape_file_path = f"{hydro_shp_path}/River_network.shp"  # The shape file that defines the extent
cdf_files = Path(netcdf_files_path).glob('*.nc')  # generator to capture all the netCDF files in a folder


# def bbox(lat,lng, margin):
#     return Polygon([[lng-margin, lat-margin],[lng-margin, lat+margin],
#     [lng+margin,lat+margin],[lng+margin,lat-margin]])

def create_bbox_geometry(extent_shape_path: str):
    """Create a bounding box Polygon for an input shape file.

    :param str extent_shape_path: Path to the shape file for which to create a bounding box.
    :return shapely.geometry.polygon.Polygon: A shapely polygon that defines the bounding box
    :return GeoDataFrame geodf: The GeoDataFrame with the data from the extent shape file
    """
    geodf = gpd.read_file(extent_shape_path)
    se = geodf.geometry.total_bounds  # Get bounds and store in array
    return Polygon([[se[0], se[1]], [se[2], se[1]],
                    [se[2], se[3]], [se[0], se[3]]]), geodf


def create_bounding_shape(ext_shape_file_path: str):
    """Create a bounding box shape file.

    :param str ext_shape_file_path: Path to the shape file for which to create a bounding box.
    :return GeoDataFrame bounding_shape: Bounding box as a GeoDataFrame.
    :return GeoDataFrame geodf:
    """
    bounding_box, geodf = create_bbox_geometry(ext_shape_file_path)
    # Create the bounding shape in a GeoDataFrame
    bounding_shape: GeoDataFrame = gpd.GeoDataFrame(pd.DataFrame(['p1'], columns=['geom']),
                                                    crs={'init': 'epsg:21781'},
                                                    geometry=[bounding_box])
    bounding_shape.to_file(f"{out_path}/bbox.shp")  # This is the bounding box as a shape file
    return bounding_shape, geodf


def netcdf_clipper(cdf_path: str, bb_shape_path: str, cdf_filename: str):
    """Clips one or more netCDF files according to a bounding box.

    :param str cdf_path: Path to the netCDF file to clip
    :param str bb_shape_path: Path to the shape file of the bounding box to be used to clip.
    :param cdf_filename: If only one file should be clipped, enter the file name here
    """

    # TODO: Write the xarray statements into a function
    # Read in the netCDF file into an xarray Dataset
    xds: Dataset = xarray.open_dataset(
        f"{netcdf_files_path}/RhiresD_ch01h.swiss.lv95_196101010000_196112310000.nc", )

    # xds = xds[['precipitationCal', 'precipitationCal_cnt']].transpose('time', 'lat', 'lon')

    xds.rio.set_spatial_dims(x_dim="E", y_dim="N", inplace=True)  # Define the dimensions
    xds.rio.write_crs("EPSG:2056", inplace=True)
    # delete the attribute 'grid_mapping' to prevent an error
    vars_list = list(xds.data_vars)
    for var in vars_list:
        del xds[var].attrs['grid_mapping']

    clipped = xds.rio.clip(bounding_shape.geometry.apply(mapping),
                           geodf.crs)  # Clip the xarray Dataset according to the bounding box GeoDataFrame
    clipped.to_netcdf(f"{out_path}/RhiresD_ch01h.swiss.lv95_196101010000_196112310000_clipped.nc",
                      "w")  # Write the clipped netCDF file


bounding_shape, geodf = create_bounding_shape(extent_shape_file_path)
for f in
netcdf_clipper(netcdf_files_path, f"{out_path}/bbox.shp")
