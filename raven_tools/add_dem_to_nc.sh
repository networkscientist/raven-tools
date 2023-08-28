#!/bin/zsh
#cellid=$1
cd '/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out' || exit 91
#ogrinfo -al grid_weights_CH-0053.shp | grep cell_id | awk '{print $4}' >/media/mainman/Work/RAVEN/data/DEM/out/grid_cells_CH-0053/cell_list.txt
while read cellid; do
  {
    gdalwarp -of VRT -cutline grid_weights_CH-0053.shp -crop_to_cutline -cwhere "cell_id = '${cellid}'" /media/mainman/Work/RAVEN/data/DEM/out/dem.vrt /media/mainman/Work/RAVEN/data/DEM/out/grid_cells_CH-0053/dem_clip_CH-0053_"${cellid}".vrt
  } || {
    echo "${cellid}" >>/media/mainman/Work/RAVEN/data/DEM/out/grid_cells_CH-0053/fails.txt
  }
done </media/mainman/Work/RAVEN/data/DEM/out/grid_cells_CH-0053/cell_list.txt
