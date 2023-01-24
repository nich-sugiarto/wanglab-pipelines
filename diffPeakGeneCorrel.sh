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
library(biomaRt)
library(ggplot2)
library(viridis)

getwd()

genes <- read.csv("${geneTable}", header = TRUE) 
glimpse(genes)

# Remove ' as a quote because of 5' and 3' annotations
peaks <- read.csv("${peakTable}", header = TRUE, row.names = 1) %>% clean_names()

matchedGenes <- genes[genes\$"gene" %in% peaks\$"gene_name", ]
matchedPeaks <- peaks[peaks\$"gene_name" %in% genes\$"gene", ]

print(nrow(matchedGenes))
print(nrow(matchedPeaks))

diffCorr <- merge(x = matchedPeaks, y = matchedGenes, by.x = "gene_name", by.y = "gene")
str(diffCorr)

write.csv(diffCorr, "corrPeakGenes/${name}.csv")

q1 <- diffCorr %>%
	filter(fold > 0 & log2FoldChange > 0)

write.table(unique(q1\$gene_name), "corrPeakGenes/quadrants/${name}_q1Genes.csv", sep = ",", col.names = FALSE, row.names = FALSE)

q2 <- diffCorr %>%
	filter(fold < 0 & log2FoldChange > 0)

write.table(unique(q2\$gene_name), "corrPeakGenes/quadrants/${name}_q2Genes.csv", sep = ",", col.names = FALSE, row.names = FALSE)

q3 <- diffCorr %>%
	filter(fold < 0 & log2FoldChange < 0)

write.table(unique(q3\$gene_name), "corrPeakGenes/quadrants/${name}_q3Genes.csv", sep = ",", col.names = FALSE, row.names = FALSE)

q4 <- diffCorr %>%
	filter(fold > 0 & log2FoldChange < 0)

write.table(unique(q4\$gene_name), "corrPeakGenes/quadrants/${name}_q4Genes.csv", sep = ",", col.names = FALSE, row.names = FALSE)

pdf("corrPeakGenes/${name}.pdf")
corrPlot <- ggplot(diffCorr, aes(x=fold, y=log2FoldChange)) + 
  geom_bin2d(bins = 120) + theme_bw() + 
  geom_smooth(method=lm) +  
  labs(y = "RNA log2fold change", x = "ATAC log2fold change") + 
  ggtitle(paste("${name} Pearson Correlation:", cor(diffCorr\$fold, diffCorr\$log2FoldChange))) + 
  geom_text(x = -3, y = 3, label = sum(diffCorr\$fold < 0 & diffCorr\$log2FoldChange > 0), colour = "red") + 
  geom_text(x = 3, y = -3, label = sum(diffCorr\$fold > 0 & diffCorr\$log2FoldChange < 0), colour = "red") + 
  geom_text(x = 3, y = 3, label = sum(diffCorr\$fold > 0 & diffCorr\$log2FoldChange > 0), colour = "red") + 
  geom_text(x = -3, y = -3, label = sum(diffCorr\$fold < 0 & diffCorr\$log2FoldChange < 0), colour = "red") + 
  scale_fill_continuous(type = "viridis")
print(corrPlot)
dev.off()
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