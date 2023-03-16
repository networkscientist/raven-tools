#!/bin/bash

X=( $(ogrinfo -so -al $1 | grep Extent | grep -Po '\d+\.\d+') )
echo "$1"
pwd
echo $X
echo "${X[1]}"
#The increment makes sure the grid covers the shape file extent completely
increment="0.01"
LonMin=$(echo "${X[0]} - ${increment}" | bc -l)
LonMax=$(echo "${X[2]} + ${increment}" | bc -l)
LatMin=$(echo "${X[1]} - ${increment}" | bc -l)
LatMax=$(echo "${X[3]} + ${increment}" | bc -l)
#LonMin="${X[0]}"
#LonMax="${X[2]}"
#LatMin="${X[1]}"
#LatMax="${X[3]}"
echo "LonMin = $LonMin"
echo "LonMax = $LonMax"
echo "LatMin = $LatMin"
echo "LatMax = $LatMax"
echo "\$1=$1"
echo "\$2=$2"
out_dir="$(dirname "${2}")"
out_file_nc_prefix="$(basename "${2}")"
out_file_nc_prefix="${out_file_nc_prefix%.*}"
out_file_shp_prefix="$(basename "${1}")"
out_file_shp_prefix="${out_file_shp_prefix%.*}"
echo "OutDir=$out_dir"
echo "OutFileNetCDFPrefix=$out_file_nc_prefix"
echo "OutFileSHPPrefix=$out_file_shp_prefix"

cdo sellonlatbox,"$LonMin","$LonMax","$LatMin","$LatMax" "$2" "${out_dir}/${out_file_nc_prefix}_${out_file_shp_prefix}_clipped.nc"
