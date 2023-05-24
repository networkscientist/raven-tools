#!/bin/bash
#### Script that merges netCDF files of different years into 1 file
START_DATE="$1"
END_DATE="$2"
FORCING_NAME="$3"
FORCING_DIR="$4"

function merge_netcdf(){
  start_date=$1
  end_date=$2
  forcing_name=$3
  forcing_dir=$4
  out_dir="$(dirname "$forcing_dir")/out/"
  declare -a IN_PATHS
  for year in $(seq "$start_date" "$end_date")
  do
	  IN_PATHS+=(${forcing_dir}/*_${year}01010000_${year}12310000.nc)
  done

  cdo mergetime "${IN_PATHS[@]}" "${out_dir}${forcing_name}_${start_date}01010000_${end_date}12310000.nc"
  unset IN_FILES
}

merge_netcdf "$START_DATE" "$END_DATE" "$FORCING_NAME" "$FORCING_DIR"