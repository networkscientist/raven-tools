#!/usr/bin/zsh
gdalwarp -srcnodata -9999.00 -r 'average' -tap -tr 1000 1000 dem_CH-0105.tif dem_CH-0105_resampled.tif
