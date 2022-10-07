"""
This module contains tools to clip netCDF files according to the extent of a shape file.

Functions:
    create_bbox_geometry(Path)
    create_bounding_shape(Path)
    netcdf_to_dataset(Path)
    dataset_to_netcdf(Dataset)
    netcdf_clipper(Path,Path,GeoDataFrame)
    netcdf_clipper_multi(Path,GeoDataFrame)
"""

import glob
from pathlib import Path
import geopandas as gpd
import pandas as pd
import xarray
from geopandas import GeoDataFrame
from shapely.geometry import mapping, Polygon
from xarray import Dataset
from netCDF4 import Dataset


def create_bbox_geometry(extent_shape_path: Path):
    """Create a bounding box Polygon for an input shape file.

    :param str extent_shape_path: Path to the shape file for which to create a bounding box.
    :return pl: A shapely polygon that defines the bounding box.
    :return gdf: The GeoDataFrame with the data from the extent shape file.
    :rtype pl: shapely.geometry.polygon.Polygon
    :rtype gdf: GeoDataFrame

    """

    gdf = gpd.read_file(extent_shape_path)
    se = gdf.geometry.total_bounds  # Get bounds and store in array
    pl = Polygon([[se[0], se[1]], [se[2], se[1]],[se[2], se[3]], [se[0], se[3]]])
    return pl, gdf  # Return the Polygon and the GeoDataFrame


def create_bounding_shape(ext_shape_file_path: Path):
    """Create a bounding box shape file.

    :param str ext_shape_file_path: Path to the shape file for which to create a bounding box.
    :return bounding_shape: Bounding box as a GeoDataFrame.
    :return geodf: The GeoDataFrame with the data from the extent shape file
    :rtype bounding_shape: GeoDataFrame
    :rtype geodf: GeoDataFrame

    """

    # Create the bounding box Polygon
    bounding_box, gdf = create_bbox_geometry(ext_shape_file_path)
    # Create the bounding shape in a GeoDataFrame
    bb_shape: GeoDataFrame = gpd.GeoDataFrame(pd.DataFrame(['p1'], columns=['geom']),
                                              crs={'init': 'epsg:21781'},
                                              geometry=[bounding_box])
    bb_shape.to_file(f"{out_path}/{bounding_box_filename}")  # This is the bounding box as a shape file
    return bb_shape, gdf


def netcdf_to_dataset(cdf_path):
    """Reads a netCDF file into an xarray dataset

    :param Path cdf_path: Path to the netCDF file to clip
    :return xds: The netCDF data as an xarray dataset
    :rtype xds: Dataset

    """

    # Read in the netCDF file into an xarray Dataset
    xds: Dataset = xarray.open_dataset(
        cdf_path)
    xds.rio.set_spatial_dims(x_dim="E", y_dim="N", inplace=True)  # Define the dimensions
    xds.rio.write_crs("EPSG:2056", inplace=True)
    # delete the attribute 'grid_mapping' to prevent an error
    vars_list = list(xds.data_vars)
    for var in vars_list:
        del xds[var].attrs['grid_mapping']
    return xds


def calculate_hamon_pet(cdf_path):
    xds = netcdf_to_dataset(cdf_path)
    xdf = xds.to_dataframe()
    return xdf


def dataset_to_netcdf(dataset, cdf_path):
    """ Writes xarray dataset to netCDF

    :param dataset dataset: xarray Dataset to write.
    :param Path cdf_path: netCDF file path to write to.

    """

    dataset.to_netcdf(f"{cdf_path.parent.parent}/out/{cdf_path.stem}_clipped{cdf_path.suffix}",
                      "w")  # Write the clipped netCDF file


def netcdf_clipper(cdf_path_in: Path, bb_shape_path: Path, gdf):
    """Clips a netCDF file according to a bounding box.

    For one netCDF file in a directory, clips it according to a bounding box shape file.

    :param GeoDataFrame gdf: The GeoDataFrame with the data from the extent shape file
    :param str cdf_path_in: Path to the netCDF file to clip
    :param str bb_shape_path: Path to the shape file of the bounding box to be used to clip.
    :return clipped: Clipped DataSet??
    :rtype clipped: ??

    """

    xds = netcdf_to_dataset(cdf_path_in)
    clip_box = gpd.read_file(f"{out_path}/{bb_shape_path}")
    clipped = xds.rio.clip(clip_box.geometry.apply(mapping),
                           gdf.crs)  # Clip the xarray Dataset according to the bounding box GeoDataFrame
    # TODO: Check if it is a dataset?
    return clipped


def netcdf_clipper_multi(ncdf_path, gdf):
    """ Clips multiple netCDF files in a directory

    :param Path ncdf_path:
    :param GeoDataFrame gdf:

    """

    for f in glob.glob(f"{ncdf_path}/original_files/*.nc"):
        netcdf_clipper(Path(f), bounding_box_filename, gdf)


if __name__ == '__main__':
    home_path = "/media/mainman/Data/"
    netcdf_files_path = Path(
        f"{home_path}/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95")  # the original netCDF files are stored here
    hydro_shp_path = Path(f"{home_path}/RAVEN/data/Hydro")  # Hydrology shape files are stored here
    out_path = Path(f"{netcdf_files_path}/out")  # The output files will be store here
    extent_shape_file_path = Path(f"{hydro_shp_path}/River_network.shp")  # The shape file that defines the extent
    cdf_files = Path(netcdf_files_path).glob('*.nc')  # generator to capture all the netCDF files in a folder
    bounding_box_filename = Path("bbox.shp")

    # with Dataset(cdf_in_path,format="NETCDF4") as cdf_dataset_in, Dataset(cdf_out_path, "w",format="NETCDF4", clobber=True) as cdf_dataset_out:
#     # copy global attributes all at once via dictionary
#     cdf_dataset_out.setncatts(cdf_dataset_in.__dict__)
#     # copy dimensions
#     for name, dimension in cdf_dataset_in.dimensions.items():
#         cdf_dataset_out.createDimension(
#             name, (len(dimension) if not dimension.isunlimited() else None))
#     # copy all file data except for the excluded
#     for name, variable in cdf_dataset_in.variables.items():
#         x = cdf_dataset_out.createVariable(name, variable.datatype, variable.dimensions)
#         cdf_dataset_out[name][:] = cdf_dataset_in[name][:]
#         # copy variable attributes all at once via dictionary
#         cdf_dataset_out[name].setncatts(cdf_dataset_in[name].__dict__)

    cdf_in_path = "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/TabsD_v2.0_swiss.lv95/out/TabsD_ch01r.swiss.lv95_200001010000_200012310000_clipped.nc"
    cdf_out_path = "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/TabsD_v2.0_swiss.lv95/out/TabsD_ch01r.swiss.lv95_200001010000_200012310000_pet.nc"
    # shutil.copyfile(cdf_in_path, cdf_out_path)
    # xds = netcdf_to_dataset("/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/TabsD_v2.0_swiss.lv95/out/TabsD_ch01r.swiss.lv95_200001010000_200012310000_clipped.nc")
    # cdf_dataset_in = Dataset("/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/TabsD_v2.0_swiss.lv95/out/TabsD_ch01r.swiss.lv95_200001010000_200012310000_clipped.nc",format="NETCDF4")
    # cdf_dataset_out = Dataset("/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/TabsD_v2.0_swiss.lv95/out/TabsD_ch01r.swiss.lv95_200001010000_200012310000_pet.nc","r+",format="NETCDF4")

    # print(f"Dimensions: {cdf_dataset_out.dimensions}")
    # print(f"Variables: {cdf_dataset_out.variables}")
    # print(f"Groups: {cdf_dataset_out.groups}")
    # pet = cdf_dataset_out.createVariable('pet', np.float32, ('time', 'N', 'E'), fill_value=-999.99)
    # pet.units = 'mm/d'
    # pet.long_name = 'daily potential evapotranspiration (Hamon)'
    # cdf_dataset_out.close()


    # xdf = xds.to_dataframe()
    # xdf = calculate_hamon_pet(
    #     "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/TabsD_v2.0_swiss.lv95/out/TabsD_ch01r.swiss.lv95_200001010000_200012310000_clipped.nc")
    # lon, lat,spref, coord, tmean = [xdf[col] for col in xdf.columns]
    # df_spref = spref.reset_index()
    # df_spref.set_index(['E','N','time'], inplace=True)
    # df_tmean = tmean.reset_index(level=[0, 1])
    #
    # df_coord = coord.reset_index()
    # df_coord.set_index(['E','N','time'], inplace=True)
    # df_lon = lon.reset_index()
    # df_lon.set_index(['E','N','time'], inplace = True)
    # df_lat = lat.reset_index(level=[0, 1])
    #
    #
    # df_lat["lat"] = df_lat["lat"] * np.pi / 180

    # dl = meteo_utils.daylight_hours(df_tmean["TabsD"].index, df_lat["lat"])
    # test = (dl / 12) ** 2 * exp(df_tmean["TabsD"] / 16)
    # test = test.to_frame()
    # test.reset_index(inplace=True)
    # df_tmean = df_tmean.reset_index()
    # xdf.reset_index(inplace=True)
    # result = pd.concat([xdf,test], ignore_index=True, axis=1)
    # result.drop([5,6,8], axis=1, inplace=True)
    # columns = ['E','N','time','lon','lat','TabsD','PET']
    # result.columns= columns
    # result['E'].attrs = {'units': 'metre'}
    # result.set_index(['E','N','time'],inplace=True)
    # xarray.Dataset(result.to_xarray()).to_netcdf("/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/TabsD_v2.0_swiss.lv95/out/TabsD_ch01r.swiss.lv95_200001010000_200012310000_pet.nc")
    # # df_lat = df_lat.reset_index()
    # df_tmean.set_index(['E','N','time'], inplace=True)
    # # df_lat.set_index(['E','N','time'], inplace=True)
    # test.rename({0:"PET"},axis=1)

    #
    # test.set_index(['E','N','time'], inplace=True)
    # df_join = df_tmean.join([test])
    # result = pd.concat([df_lon, df_lat, df_spref, df_coord, df_tmean, test])

    # result = result = pd.concat([result, df_coord], axis=1)
    # result = test.rename({0:'PET'}, axis=1)
    # result2.reset_index(inplace=True)
    # result2 = result2.set_index(['E','N','time'])
    # tmean =


    # Create the .shp file of the bounding box.
    bounding_shape, geodf = create_bounding_shape(extent_shape_file_path)
    # Clip all rainfall .nc files in the directory and store the clipped files in the output folder.
    netcdf_clipper_multi(netcdf_files_path, geodf)

    netcdf_files_path = Path(
        f"{home_path}/RAVEN/data/MeteoSwiss_gridded_products/SrelD_v2.0_swiss.lv95")  # the original netCDF files are stored here
    netcdf_clipper_multi(netcdf_files_path, geodf)

    netcdf_files_path = Path(
        f"{home_path}/RAVEN/data/MeteoSwiss_gridded_products/TabsD_v2.0_swiss.lv95")  # the original netCDF files are stored here
    netcdf_clipper_multi(netcdf_files_path, geodf)

    netcdf_files_path = Path(
        f"{home_path}/RAVEN/data/MeteoSwiss_gridded_products/TmaxD_v2.0_swiss.lv95")  # the original netCDF files are stored here
    netcdf_clipper_multi(netcdf_files_path, geodf)

    netcdf_files_path = Path(
        f"{home_path}/RAVEN/data/MeteoSwiss_gridded_products/TminD_v2.0_swiss.lv95")  # the original netCDF files are stored here
    netcdf_clipper_multi(netcdf_files_path, geodf)
