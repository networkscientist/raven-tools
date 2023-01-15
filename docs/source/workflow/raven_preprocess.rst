==============
Pre-Processing
==============

Bounding Box
--------------

Creating a bounding box can be done in two ways:

    1. By creating a shapely Polygon:

        Use :func:`raven_preprocess.create_bbox_geometry`. This will not write the result to any file.

    2. By creating a shape file:

        Use :func:`raven_preprocess.create_bounding_shape`. This will write the result to a shape file in addition to returning it as a GeoDataFrame.


netCDF File Handling
--------------------

Merging netCDF files
^^^^^^^^^^^^^^^^^^^^

Merging netCDF files is best done with the Bash library *cdo*. Therefore, if you call :meth:`raven_preprocess.nc_merge` this will start an external Bash script *nc_combine.sh* to make use of *cdo*, which is way faster than a Python implementation. Simply put the desired files into a folder and pass its path and the time range to :meth:`raven_preprocess.nc_merge`. You will find the resulting file in the same folder.

Clipping netCDF files
^^^^^^^^^^^^^^^^^^^^^

To save disk space and computation time, you may want to clip your netCDF files' extents, typically with a bounding box of your catchment/research area. There are two corresponding functions to do so:

    1. To clip a single netCDF file, use :meth:`raven_preprocess.netcdf_clipper`. The original netCDF file will not be overwritten. Clipped files can be found in the *out/* folder. In addition, the function returns the Dataset.

    2. To clip multiple netCDF files in a directory, use :meth:`raven_preprocess.netcdf_clipper_multi`. This function uses globbing and processes any .nc file in the directory. The output files can be found in the *out/* folder.

PET calculations
^^^^^^^^^^^^^^^^

Although Raven provides built-in routines for PET estimation, you may want to create netCDF files on your own. One possible approach is implemented in :meth:`raven_preprocess.netcdf_pet_hamon`: It will append PET estimated according to Hamon to a netCDF file containing forcing data. This function creates a new file and does not overwrite the original one.

netCDF conversions
^^^^^^^^^^^^^^^^^^

To work with netCDF files in Python, there is the great `netCDF4 <https://pypi.org/project/netCDF4/>`_ module. Make use of it! To help you, there are two functions:

    1. :meth:`raven_preprocess.netcdf_to_dataset` to read in .nc files into a Dataset,

    2. :meth:`raven_preprocess.dataset_to_netcdf` to write a Dataset to a .nc file in the *out/* folder.

Grid Weights
------------

Discharge
---------
