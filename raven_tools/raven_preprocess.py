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

    This function takes a shape file, creates a bounding box around the features and returns this box as a
    GeoDataFrame. Additionally, returns the original extent as a GeoDataFrame.

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

    Parameters
    ----------
    netcdf_file_path : Path
        Path to the netCDF file to clip

    Returns
    -------
    xds : Dataset
        The netCDF data as a netCDF4 dataset

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

    Parameters
    ----------
    xds_to_write : Dataset
        xarray Dataset to write.
    netcdf_file_path : Path
        netCDF file path to write to.

    """

    # Write the clipped netCDF file
    xds_to_write.to_netcdf(
        f"{netcdf_file_path.parent.parent}/out/{netcdf_file_path.stem}_clipped{netcdf_file_path.suffix}", "w")


def netcdf_clipper(netcdf_file_path: Path, bbox_file_path: Path, ext_gdf: GeoDataFrame):
    """Clips a netCDF file according to a bounding box.

    For one netCDF file in a directory, clips it according to a bounding box shape file.

    Parameters
    ----------
    ext_gdf : GeoDataFrame
        The GeoDataFrame with the data from the extent shape file
    netcdf_file_path : Path
        Path to the netCDF file to clip
    bbox_file_path : Path
        Path to the shape file of the bounding box to be used to clip.

    Returns
    -------
    xds_clipped : Dataset
        The clipped Dataset

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

    Parameters
    ----------
    bbox_file_path : Path
        Full Path to the bounding box shape file
    netcdf_dir_path : Path
        Path to directory with netCDF files to clip.
    bbox_gdf : GeoDataFrame
        Bounding box GeoDataFrame created with create_bounding_shape()

    """

    for f in glob.glob(f"{netcdf_dir_path}/original_files/*.nc"):
        netcdf_clipper(Path(f), bbox_file_path, bbox_gdf)


def netcdf_pet_hamon(netcdf_file_path: Path, name_pattern: dict[str, str]):
    """

    Parameters
    ----------
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


def nc_merge(start_year: int, end_year: int, forcing_dir: str):
    """

    Parameters
    ----------
    start_year : int
        Start year
    end_year : int
        End Year
    forcing_dir : str
        Root directory where forcing files are located

    """
    subprocess.call(['raven_tools/nc_combine.sh', str(start_year), str(end_year), forcing_dir])


def create_grid(export_shp = True):
    # Since the netCDF has a cell size of 1000m, set this here
    length: int = 1000
    wide: int = 1000

    # Create the northing and easting coordinates of the bounding vertex for each cell
    cols: list = list(np.arange(xmin, xmax + wide, wide))
    rows: list = list(np.arange(ymin, ymax + length, length))

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

    if export_shp:
        # Export the grid to a shape file
        grid.to_file("/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/grid.shp")


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
    res_d: GeoDataFrame = grd.overlay(ctm, how='difference')
    res_d = res_d.rename_geometry('geometry')
    res_u.set_index("cell_id")
    res_d.set_index("cell_id")
    return res_u, res_d


def calc_relative_area(gdf: GeoDataFrame):
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
    gdf = gdf.to_crs(32632)
    # For each intersected feature, compute the area and write into a new field
    gdf["area"] = gdf["geometry"].area
    # Convert the area to float
    area_sum: float = float(gdf["area"].sum())
    for index, row in gdf.iterrows():
        # To minimize numerical errors, multiply and then divide the result by 1000000
        gdf.at[index, "area_rel"] = ((gdf.loc[index]["area"] * 1000000) / area_sum) / 1000000
        gdf.at[index, "index"] = int(index)
    # Re-project back into the LV95 projection and export to a shape file
    gdf: GeoDataFrame = gdf.to_crs(2056)
    return gdf


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
    data_to_write : list[list[float]]
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
        "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/GridWeights.txt")
    bounding_box_filename = Path(
        "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/bbox.shp")
    # The shape file that defines the extent of the catchment
    extent_shape_file_path = Path("/media/mainman/Data/RAVEN/data/Catchment/Broye_Payerne.shp")
    # The netCDF file that will be used to get the bounds/extent
    netCDF_file_path = Path(
        "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out"
        "/RhiresD_ch01h.swiss.lv95_196201010000_196212310000_clipped.nc")

    # Read in the bounding box shape file from the clipping folder as a GeoPandas DataFrame
    bbox: GeoDataFrame = gpd.read_file(bounding_box_filename)
    bbox.set_crs(21781, allow_override=True)
    bbox = bbox.to_crs("EPSG:2056")

    # Read in the clipped netCDF file into a xarray Dataset
    ds: Dataset = xr.open_dataset(netCDF_file_path, engine="netcdf4")

    # Since we're only interested in the grid (not the actual cell values), only use 1 day
    ds = ds.sel(time=slice('1962-01-01', '1962-01-02'))
    # Select the actual data value column
    xarr: xr.DataArray | xr.Dataset = ds['RhiresD']
    # Convert Dataset to DataFrame and reset the index
    df: pd.DataFrame = xarr.to_dataframe().reset_index()
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
    cols: list = list(np.arange(xmin, xmax + wide, wide))
    rows: list = list(np.arange(ymin, ymax + length, length))

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
    print()