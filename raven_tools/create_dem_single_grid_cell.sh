#!/bin/bash
gdalwarp -s_srs EPSG:21781 -t_srs EPSG:2056 -te 2740000 1144000 2741000 1145000 -te_srs EPSG:2056 -tr 1000 1000 -tap -srcnodata -9999 -r 'average' -overwrite "/media/mainman/Work/RAVEN/data/DEM/shape/ASCII_GRID_1part/dhm25_grid_raster.asc" "/media/mainman/Work/RAVEN/data/DEM/shape/ASCII_GRID_1part/dhm25_grid_raster_2056.tif"
