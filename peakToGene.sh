#!/bin/bash

# Takes in the annotated output of the annotatepeaks.pl on ATAC data as well as the DESEQ2 output of rna_twoTwo.sh or rna_twoThree.sh
# and outputs DESEQ regions that are found in both

# USAGE:
# This pipeline takes in one positional argument:
# 	$1 - target meta file containing sample information
# Where the meta file is tab-delimited, and formated as follows:
# homerAnnotateOutput	deseqOutput	nameOfFile

# FIXME: Create new, more detailed spreadsheet instead.

folder=$(cd "$(dirname "$0")";pwd)  # Save current folder as a variable

# Setup - generate required folders
mkdir -p PBS
mkdir -p log
mkdir -p linkedGenes

# Error message if the files aren't provideds
if [ -z "$1" ]; then 
  echo ERROR: META FILE WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target meta file containing sample information
  exit 1
fi

folder=$(cd "$(dirname "$0")";pwd)  # Save location of current folder
libraryText=$1  # Passed in file name

OLDIFS=$IFS
while IFS=$'\t' read -r -a varArray; do
	peakTable=${varArray[0]}  # The control group
	geneTable=${varArray[1]}  # The treatment group
	name=${varArray[2]}  # Name of resulting spreadsheet
	IFS=$OLDIFS

	echo ${name}

	cat >${folder}/PBS/${name}_linker'.R' <<EOF
library(tidyverse)
library(dplyr)
library(janitor)

getwd()

genes <- read.csv("${geneTable}", header = TRUE) 
glimpse(genes)

peaks <- as.list(read.table("${peakTable}", sep = "\t", header = TRUE, fill = TRUE) %>% clean_names())
glimpse(peaks)

peakNames <- unique(peaks\$"gene_name")

longNames <- peakNames

altNames <- unique(peaks\$"gene_alias")
for (gList in altNames) {
	longNames <- append(longNames, str_split(gList, pattern = "|"))
}

longNames <- unique(longNames)

print(peakNames[1:50])
df <- genes[ genes\$"gene" %in% peakNames, ]
df2 <- genes[ genes\$"gene" %in% longNames, ]

print(paste0("There are ", (length(longNames) - length(peakNames)), " more items after adding alternative gene names."))
print(paste0("There are ", (nrow(df2) - nrow(df)), " more rows after adding alternative gene names."))

write.csv(df, "./linkedGenes/${name}_linkedPeakGenes.csv")
write.csv(df2, "./linkedGenes/${name}_linkedPeakGenes_long.csv")
EOF

	cat >${folder}/PBS/${name}_linker.pbs <<EOF
#!/bin/bash
#SBATCH --job-name=${name}_linker
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${name}_linker_%j.txt -e ${folder}/log/${name}_linker_%j.err.txt
cd ${folder}
#------- End of header -------#
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate ChIPseeker

Rscript ${folder}/PBS/${name}_linker.R
EOF

	sbatch ${folder}/PBS/${name}_linker.pbs
done < ${libraryText}