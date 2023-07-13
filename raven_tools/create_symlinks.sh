#!/bin/bash
#### Script that creates symlinks to data dirs 

FORCING_DIR_SRC="/storage/homefs/pz09y074/raven_master_files/RAVEN/data/MeteoSwiss_gridded_products"
MODEL_DIR="/home/sirian/Applications/Hydrology/RAVEN/models"
RAIN_FILE="RhiresD_ch01h.swiss.lv95"
T_MEAN_FILE="TabsD_ch01r.swiss.lv95"
T_MAX_FILE="TmaxD_ch01r.swiss.lv95"
T_MIN_FILE="TminD_ch01r.swiss.lv95"
S_REL_FILE="SrelD_ch01r.swiss.lv95"

FILES=($RAIN_FILE $T_MEAN_FILE $T_MAX_FILE $T_MIN_FILE $S_REL_FILE)
FILES_LEN=$((${#FILES[@]}-1))
  # Take in the start and end year as arguments
create_symlinks () {
  # Iterates over each directory in the DIRS array
  for ccid in $MODEL_DIR/*
  do
    for mod in $ccid/*
    do
      for forc in "${FILES[@]}"
      do
#        echo -e "SRC: $FORCING_DIR_SRC/$forc"
#        echo -e "DST: $mod/model/data_obs/$forc"
        ln -sfn "$FORCING_DIR_SRC/$forc" "$mod/model/data_obs/$forc"
      done
    done
  done
}

create_symlinks

