#!/bin/bash

# The first approach uses linklist downloaded from swisstopo
dem_dir="/media/mainman/Work/RAVEN/data/DEM/"
file="linklist_CH-0058.csv"

#function convert_linklist () {
#  while read line
#  do
#    base="$(basename "$line")"
#    gdalwarp -tr 10 -10 "${dem_dir}${base}" "${dem_dir}out/${base}"
#  done < "$dem_dir$file"
#}

# The second approach converts all the dem files in the DEM folder
#counter=0
function convert_all_dem_files () {
  for f in "${dem_dir}"swissalti3d_*.tif
  do
    base="$(basename "$f")"
    gdalwarp -tr 10 -10 "${f}" "${dem_dir}out/${base}"
#    echo "${f}"
#    let "counter+=1"
  done
#  echo $counter
}

convert_all_dem_files
#