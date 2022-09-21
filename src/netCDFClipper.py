import glob
from pathlib import Path
import geopandas as gpd
import pandas as pd
import xarray
from geopandas import GeoDataFrame
from shapely.geometry import mapping, Polygon
from xarray import Dataset

# Paths
home_path = Path.home()
netcdf_files_path = Path(
    f"{home_path}/Applications/Hydrology/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95")  # the original netCDF files are stored here
hydro_shp_path = Path(f"{home_path}/Applications/Hydrology/RAVEN/data/Hydro")  # Hydrology shape files are stored here
out_path = Path(f"{netcdf_files_path}/out")  # The output files will be store here
extent_shape_file_path = Path(f"{hydro_shp_path}/River_network.shp")  # The shape file that defines the extent
cdf_files = Path(netcdf_files_path).glob('*.nc')  # generator to capture all the netCDF files in a folder
bounding_box_filename = Path("bbox.shp")


def create_bbox_geometry(extent_shape_path: Path):
    """Create a bounding box Polygon for an input shape file.

    :param str extent_shape_path: Path to the shape file for which to create a bounding box.
    :return shapely.geometry.polygon.Polygon: A shapely polygon that defines the bounding box
    :return GeoDataFrame gdf: The GeoDataFrame with the data from the extent shape file
    """
    gdf = gpd.read_file(extent_shape_path)
    se = gdf.geometry.total_bounds  # Get bounds and store in array
    return Polygon([[se[0], se[1]], [se[2], se[1]],
                    [se[2], se[3]], [se[0], se[3]]]), gdf  # Return the Polygon and the GeoDataFrame


def create_bounding_shape(ext_shape_file_path: Path):
    """Create a bounding box shape file.

    :param str ext_shape_file_path: Path to the shape file for which to create a bounding box.
    :return GeoDataFrame bounding_shape: Bounding box as a GeoDataFrame.
    :return GeoDataFrame geodf: The GeoDataFrame with the data from the extent shape file
    """
    # Create the bounding box Polygon
    bounding_box, gdf = create_bbox_geometry(ext_shape_file_path)
    # Create the bounding shape in a GeoDataFrame
    bb_shape: GeoDataFrame = gpd.GeoDataFrame(pd.DataFrame(['p1'], columns=['geom']),
                                              crs={'init': 'epsg:21781'},
                                              geometry=[bounding_box])
    bb_shape.to_file(f"{out_path}/{bounding_box_filename}")  # This is the bounding box as a shape file
    return bb_shape, gdf


def netcdf_clipper(cdf_path_in: Path, bb_shape_path: Path, gdf):
    """Clips a netCDF file according to a bounding box.

    For one netCDF file in a directory, clips it according to a bounding box .shp file.
    :param gdf: The GeoDataFrame with the data from the extent shape file
    :param str cdf_path_in: Path to the netCDF file to clip
    :param str bb_shape_path: Path to the shape file of the bounding box to be used to clip.
    """

    # Read in the netCDF file into an xarray Dataset
    xds: Dataset = xarray.open_dataset(
        cdf_path_in)

    # xds = xds[['precipitationCal', 'precipitationCal_cnt']].transpose('time', 'lat', 'lon')

    xds.rio.set_spatial_dims(x_dim="E", y_dim="N", inplace=True)  # Define the dimensions
    xds.rio.write_crs("EPSG:2056", inplace=True)
    # delete the attribute 'grid_mapping' to prevent an error
    vars_list = list(xds.data_vars)
    for var in vars_list:
        del xds[var].attrs['grid_mapping']

    clip_box = gpd.read_file(f"{out_path}/{bb_shape_path}")
    clipped = xds.rio.clip(clip_box.geometry.apply(mapping),
                           gdf.crs)  # Clip the xarray Dataset according to the bounding box GeoDataFrame
    clipped.to_netcdf(f"{cdf_path_in.parent.parent}/out/{cdf_path_in.stem}_clipped{cdf_path_in.suffix}",
                      "w")  # Write the clipped netCDF file


def netcdf_clipper_multi(ncdf_path, gdf):
    for f in glob.glob(f"{ncdf_path}/original_files/*.nc"):
        netcdf_clipper(Path(f), bounding_box_filename, gdf)


# Create the .shp file of the bounding box.
bounding_shape, geodf = create_bounding_shape(extent_shape_file_path)
# Clip all rainfall .nc files in the directory and store the clipped files in the output folder.
netcdf_clipper_multi(netcdf_files_path, geodf)

netcdf_files_path = Path(
    f"{home_path}/Applications/Hydrology/RAVEN/data/MeteoSwiss_gridded_products/SrelD_v2.0_swiss.lv95")  # the original netCDF files are stored here
netcdf_clipper_multi(netcdf_files_path, geodf)

netcdf_files_path = Path(
    f"{home_path}/Applications/Hydrology/RAVEN/data/MeteoSwiss_gridded_products/TabsD_v2.0_swiss.lv95")  # the original netCDF files are stored here
netcdf_clipper_multi(netcdf_files_path, geodf)

netcdf_files_path = Path(
    f"{home_path}/Applications/Hydrology/RAVEN/data/MeteoSwiss_gridded_products/TmaxD_v2.0_swiss.lv95")  # the original netCDF files are stored here
netcdf_clipper_multi(netcdf_files_path, geodf)

netcdf_files_path = Path(
    f"{home_path}/Applications/Hydrology/RAVEN/data/MeteoSwiss_gridded_products/TminD_v2.0_swiss.lv95")  # the original netCDF files are stored here
netcdf_clipper_multi(netcdf_files_path, geodf)