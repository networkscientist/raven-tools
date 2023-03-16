#!/bin/bash
X=( $(echo "$string" | grep -Po '?<=^| |\()\d+\.\d+(?=$| |,|\)') )
