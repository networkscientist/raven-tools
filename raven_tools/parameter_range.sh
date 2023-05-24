#!/bin/bash
MODEL="MOHYSE"
for f in {01..10}
do
  echo -e "f\"{params['${MODEL}']['names']['${MODEL}_Param_${f}']}		random		{params['${MODEL}']['lower']['${MODEL}_Param_${f}']}		{params['${MODEL}']['upper']['${MODEL}_Param_${f}']}	none   none 	none\","
done