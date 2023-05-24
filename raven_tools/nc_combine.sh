#!/bin/bash
#### Script that merges netCDF files of different years into 1 file.
FORCING_DIR=$3
CATCHMENT=$4
RAIN_DIR="${FORCING_DIR}/RhiresD_v2.0_swiss.lv95/"
RAIN_FILE="RhiresD_ch01h.swiss.lv95_"
T_MEAN_FILE="TabsD_ch01r.swiss.lv95_"
T_MAX_FILE="TmaxD_ch01r.swiss.lv95_"
T_MIN_FILE="TminD_ch01r.swiss.lv95_"
S_REL_FILE="SrelD_ch01r.swiss.lv95_"
T_MEAN_DIR="${FORCING_DIR}/TabsD_v2.0_swiss.lv95/"
T_MAX_DIR="${FORCING_DIR}/TmaxD_v2.0_swiss.lv95/"
T_MIN_DIR="${FORCING_DIR}/TminD_v2.0_swiss.lv95/"
S_REL_DIR="${FORCING_DIR}/SrelD_v2.0_swiss.lv95/"

DIRS=($RAIN_DIR $T_MEAN_DIR $T_MAX_DIR $T_MIN_DIR $S_REL_DIR)
FILES=($RAIN_FILE $T_MEAN_FILE $T_MAX_FILE $T_MIN_FILE $S_REL_FILE)

#### merge_netcdf BEGIN
# Merges netCDF files of different years into one file.
# GLOBALS:
#        FORCING_DIR
#        RAIN_DIR
#        RAIN_FILE
#        T_MEAN_FILE
#        T_MAX_FILE
#        T_MIN_FILE
#        S_REL_FILE
#        T_MEAN_DIR
#        T_MAX_DIR
#        T_MIN_DIR
#        S_REL_DIR
# ARGUMENTS:
#        Start year as an integer
#        End year as an integer
# OUTPUTS:
#        Merged netCDF files
# Return:
#        0 if success, non-zero otherwise.
### FUNCTION END
function merge_netcdf(){
  # Take in the start and end year as arguments
  start_date=$1
  end_date=$2


  # Create a counter to use in the next for loop
  COUNTER=0
  # Iterates over each directory in the DIRS array
  for i in "${DIRS[@]}"
  do
    # Create empty arrays
    declare -a IN_PATHS
    declare -a DATE_RANGE
    # For loop that iterates over the years
    for d in $(seq $start_date $end_date)
    do
      # Appends the year to the DATE_RANGE array
      DATE_RANGE+=("$d")
      # Appends created filename to the input file paths array IN_PATHS
      IN_PATHS+=("${i}out/${CATCHMENT}/${FILES[$COUNTER]}${d}01010000_${d}12310000_${CATCHMENT}_clipped.nc")
    done
    echo ${FILES[$COUNTER]}
    # uses the cdo command and expands the IN_PATHS to use every value as arguments
    cdo -mergetime "${IN_PATHS[@]}" "${i}merged/${FILES[$COUNTER]}${DATE_RANGE[0]}01010000_${DATE_RANGE[-1]}12310000_${CATCHMENT}_clipped.nc"
    ((COUNTER=COUNTER+1))
    unset DATE_RANGE
    unset IN_PATHS
  done
}

merge_netcdf $1 $2

