#!/bin/bash

# Takes in the annotated output of the annotatepeaks.pl on ATAC data as well as the DESEQ2 output of rna_twoTwo.sh or rna_twoThree.sh
# and outputs DESEQ regions that are found in both

# USAGE:
# This pipeline takes in one positional argument:
# 	$1 - target meta file containing sample information
# Where the meta file is tab-delimited, and formated as follows:
# homerAnnotateOutput	deseqOutput	nameOfFile

# FIXME: Create new, more detailed spreadsheet instead.

# Error message if the files aren't provideds
if [ -z "$1" ]; then 
  echo ERROR: META FILE WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target meta file containing sample information
  exit 1
fi

# Setup - generate required folders
mkdir -p PBS
mkdir -p log
mkdir -p linkedGenes

folder=$(cd "$(dirname "$0")";pwd)  # Save location of current folder
libraryText=$1  # Passed in file name

OLDIFS=$IFS
while IFS=$'\t' read -r -a varArray; do
	linkedGeneList=${varArray[0]}  # The control group
	geneTable=${varArray[1]}  # The treatment group
	name=${varArray[2]}  # Name of resulting spreadsheet
	IFS=$OLDIFS

	echo ${name}

	cat >${folder}/PBS/${name}_linker'.R' <<EOF
library(tidyverse)
library(dplyr)
library(janitor)
library(biomaRt)

getwd()

genes <- read.csv("${geneTable}", header = TRUE) 
glimpse(genes)

# Remove ' as a quote because of 5' and 3' annotations
peaks <- read.table("${peakTable}", sep = "\t", header = TRUE, quote = "", fill = TRUE) %>% clean_names()

mRNA <- peaks[str_sub(peaks\$"nearest_promoter_id", 2, 2)  == "M", ]

matchedGenes <- genes[genes\$"gene" %in% mRNA\$"gene_name", ]
matchedPeaks <- mRNA[mRNA\$"gene_name" %in% genes\$"gene", ]

print(nrow(matchedGenes))
print(nrow(matchedPeaks))

geneNames <- genes\$gene

write.csv(merge(x = matchedPeaks, y = matchedGenes, by.x = "gene_name", by.y = "gene"), "./linkedGenes/${name}_linkedPeakGenes.csv")

print(paste0("There are ", nrow(peaks), " peaks. Of which, ", 
	nrow(mRNA), " are linked to " , length(unique(mRNA\$"gene_name")) ," protein-coding genes. ", 
	(nrow(peaks) - nrow(mRNA)), " had no gene symbol attached to their name and were therefore tossed out."))
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