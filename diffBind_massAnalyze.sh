#!/bin/bash

# Calls most recent version of diffBind_analyze.sh, but also takes in a meta file to mass process the samples

# This pipeline requires one positional argument:
# 	$1 - The meta file that has all the comparisons that need to be made
# The file should be tab-delimited, and contain the following fields:
# categoryForComparison sampleGroup controlGroup  rObjToRead

# Requires the ATACQC pipeline (will rename in the future)

# TODO: Merge with diffBind_analyze.sh such that if only one field is provided, assume that it is a meta file
# But, if three categories are provided, use those three categories to run as it is right now. 
# Break on two

#  All the parameters are not provided
if [ -z "$1" ]; then 
  echo ERROR: NOT ALL PARAMETERS HAVE BEEN SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo $1 - The meta file that has all the comparisons that need to be made
  echo The file should be tab-delimited, and contain the following fields:
  echo categoryForComparison sampleGroup controlGroup  rObjToRead
fi

libraryText=$1  # Passed in file name

cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/diffBind_analyze.sh ./

OLDIFS=$IFS
while IFS=$'\t' read -r -a varArray; do
	category=${varArray[0]}  # Category
	sample=${varArray[1]}  # The treatment group
	control=${varArray[2]}  # The control group
  objName=${varArray[3]}  # Name of object
	IFS=$OLDIFS

  sh diffBind_analyze.sh ${category} ${sample} ${control} ${objName}
done < ${libraryText}