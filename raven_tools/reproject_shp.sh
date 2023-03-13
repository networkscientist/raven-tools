#!/bin/bash
ctm_folder="/media/mainman/Work/RAVEN/data/Catchment"
cd "$ctm_folder" || exit
for f in CH*.shp
do
  echo "$f"
  echo reproject_2056/"$f"
  ogr2ogr -t_srs EPSG:2056 reproject_2056/"$f" "$f"
done