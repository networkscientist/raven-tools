#!/usr/bin/zsh
GRIDFILE_PREFIX="/media/mainman/Work/RAVEN/data/MeteoSwiss_gridded_products/RhiresD_v2.0_swiss.lv95/out/grid_weights_CH-"
GRIDFILE_SUFFIX=".shp"
NCFILE="/media/mainman/Work/RAVEN/data/DEM/out/dem.vrt"
CTMARRAY=("0053" "0057" "0058" "0083" "0105" "0118" "0139" "0140" "0159" "0161" "0198")

for ctm in "${CTMARRAY[@]}"; do
  X=($(ogrinfo -so -al "${GRIDFILE_PREFIX}${ctm}${GRIDFILE_SUFFIX}" | grep Extent: | grep -Po "\d+\.\d+"))
  gdalwarp -srcnodata -9999.00 -r 'average' -tap -tr 1000 1000 -te "${X[1]}" "${X[2]}" "${X[3]}" "${X[4]}" "${NCFILE}" "${NCFILE%.*}_gridcells_${ctm}.vrt"
  gdal_translate -of GTiff "${NCFILE%.*}_gridcells_${ctm}.vrt" "${NCFILE%.*}_gridcells_CH-${ctm}.tif"
done
# shellcheck disable=SC2086
