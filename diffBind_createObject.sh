#!/bin/bash

# Pipeline that generates a DiffBind object, and saves it once DiffBind finishes the longest step of counting reads. 

# This pipeline requires two positional arguments:
# 		$1 - location of the meta text in question
# 		$2 - The name of the generated DiffBind object (default is dbObj)
# Requires the workflowr environment only for its installation of DiffBind

# Version Doc: https://docs.google.com/document/d/1CNl0aheyB66F-n6T-ILwqdSCW15p8n4eO7_81I-a-A0/edit

folder=$(cd "$(dirname "$0")";pwd)  # Stores current folder path as directory

# If a name is not provided
if [ -z "$1" ]; then 
  echo ERROR: TARGET FOLDER WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in two positional arguments:
  echo 	\$1 - target folder
  echo 	"\$2 - The name of the generated DiffBind object (default is dbObj)"
  exit 1
fi

#  If a name is not provided
if [ -z "$2" ]; then 
	NAME="dbObj"  # Set to default name
else
	NAME="$2"  # Otherwise use specified name
fi

# Create R Script
cat >${folder}/PBS/diffBind_all'.R' <<EOF
library(DiffBind)

getwd()
samples <- read.csv('$1')
dbObj <- dba(sampleSheet=samples)

dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel=FALSE)  # Does not run in parallel due to file size

saveRDS(dbObj, file = "$NAME.rds")
EOF

# Create slurm file to run said R script
cat >${folder}/PBS/diffBind_all.sbatch <<EOF
#!/bin/bash -l
#SBATCH --job-name=diffBind_all
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${NAME}_%j.txt -e ${NAME}_%j.err.txt

#------- END OF HEADER -------#
source activate workflowr

cd ${folder}

Rscript ./PBS/diffBind_all.R
EOF
cd log
sbatch ../PBS/diffBind_all.sbatch