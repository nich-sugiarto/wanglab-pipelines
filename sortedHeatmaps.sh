#!/bin/bash

mkdir -p PBS
mkdir -p log

# Check to makes sure that the appropriate arguments were passed in
if [ -z "$2" ]; then 
  echo ERROR: META FILE WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in two positional argument:
  echo 	\$1 - The mode that you wanted to run it with
  echo Choices: default	diffBind	proxDistal
  echo If you don't know what those mean, then you're not ready to use this pipeline
  echo \$2 - The meta text you want to pass in
  exit 1
fi

if [ "$1" == default ]; then
    sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/backend/sortedHeatmaps/sortDefault.sh $2
  elif [ "$1" == diffBind ]; then
    sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/backend/sortedHeatmaps/sortDiffBind.sh $2
  elif [ "$1" == proxDistal ]; then
    sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/backend/sortedHeatmaps/sortProximalDistal.sh $2
else
  echo ERROR: PROPER NORMALIZATION WAS NOT SPECIFIED
  echo You have three options: 
  echo default	diffBind	proxDistal
fi
