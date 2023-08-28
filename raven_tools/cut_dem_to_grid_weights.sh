#!/bin/bash
ctm_list=("CH-0053" "CH-0057" "CH-0058" "CH-0083" "CH-0105" "CH-0118" "CH-0139" "CH-0140" "CH-0159" "CH-0161" "CH-0198")
for ctm in "${ctm_list[@]}"; do
  gdalwarp -of GTiff -cutline "/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/grid_weights_${ctm}.shp" -crop_to_cutline '/media/mainman/Work/RAVEN/data/DEM/out/dem.vrt' "/media/mainman/Work/RAVEN/data/DEM/out/dem_${ctm}.tif"
done
