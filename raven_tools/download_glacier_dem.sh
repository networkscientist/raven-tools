#!/bin/bash

dem_dir="/media/mainman/Work/RAVEN/data/DEM/"
file="linklist_CH-0058.csv"

while read line
do
  wget -nc  $line
done < "$dem_dir$file"