#!/bin/bash

FILENAME=$1
PARENT_FOLDER="$(dirname "${FILENAME}")"
PARENT_FOLDER="$(dirname "${PARENT_FOLDER}")"
base="$(basename "${FILENAME}")"
{
  gdaldem slope "${FILENAME}" "${PARENT_FOLDER}/slopes/${base%.*}_slopes.tif" &&
    gdaldem aspect "${FILENAME}" "${PARENT_FOLDER}/aspects/${base%.*}_aspects.tif" &&
    exit 0
} || {
  exit 91
}
