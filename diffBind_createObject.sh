#!/bin/bash

# Pipeline that generates a DiffBind object, and saves it once DiffBind finishes the longest step of counting reads. 

# This pipeline requires two positional arguments:
# 		$1 - location of the meta text in question
# 		$2 - The name of the generated DiffBind object (default is dbObj)
# Requires the ATACQC because apparently that's where DiffBind is
# TODO: create separate DiffBind environment?

mkdir -p diffBind
mkdir -p PBS
mkdir -p log

folder=$(pwd)  # Stores current folder path as directory

# If a name is not provided
if [ -z "$1" ]; then 
  echo ERROR: TARGET CSVs AND NORMALIZATION METHOD WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in two positional arguments:
  echo 	"\$1 - target csv"
  echo 	"\$2 - Normalization Strategy"
  echo 	"\$3 - (Optional) The name of the generated DiffBind object (default is dbObj)"
  exit 1
fi

#  If a name is not provided
if [ -z "$3" ]; then 
	NAME="dbObj"  # Set to default name
else
	NAME="$3"  # Otherwise use specified name
fi

norm=""
off=FALSE
if [ "$2" == Default ]; then
    norm=""
  elif [ "$2" == TMM ]; then
    norm=DBA_NORM_TMM
  elif [ "$2" == RLE ]; then
    norm=DBA_NORM_RLE
  elif [ "$2" == Loess ]; then
    norm=DBA_NORM_OFFSETS
    off=TRUE
else
  echo ERROR: PROPER NORMALIZATION WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in three positional arguments:
  echo 	"\$1 - target folder"
  echo 	"\$2 - Normalization Strategy (Default, TMM, RLE, Loess)"
  echo 	"\$3 - (Optional) The name of the generated DiffBind object (default is dbObj)"
  exit 1
fi

# Create R Script
cat >${folder}/PBS/diffBind_${NAME}_${2}'.R' <<EOF
library(DiffBind)

samples <- read.csv('$1')
dbObj <- dba(sampleSheet=samples)

dbObj <- dba.blacklist(dbObj, blacklist = DBA_BLACKLIST_GRCH38)
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel=FALSE)  # Does not run in parallel due to file size
dbObj <- dba.normalize(dbObj, method = DBA_ALL_METHODS, normalize=${norm}, offsets=${off})

dba.save(dbObj, file = paste0("${NAME}","_","${2}"), dir = "diffBind")

pdf("diffBind/multiCorr_${NAME}_${2}.pdf")
  plot(dbObj)
dev.off()

pdf("diffBind/pca_${NAME}_${2}.pdf")
  dba.plotPCA(dbObj,DBA_FACTOR,label=DBA_ID)
dev.off()
EOF

# Create slurm file to run said R script
cat >${folder}/PBS/diffBind_${NAME}_${2}.sbatch <<EOF
#!/bin/bash -l
#SBATCH --job-name=diffBind_all
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${NAME}_%j.txt -e ${folder}/log/${NAME}_%j.err.txt

#------- END OF HEADER -------#
source activate ATACQC

cd ${folder}

Rscript ./PBS/diffBind_${NAME}_${2}.R
EOF

sbatch ${folder}/PBS/diffBind_${NAME}_${2}.sbatch


