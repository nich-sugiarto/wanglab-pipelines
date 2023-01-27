#!/bin/bash/

# Creates bar charts from bulk RNA sequencing

# USAGE:
# This pipeline takes in one positional argument:
# 	$1 - target meta file containing all the csvs to pull from

# There is only one field in this meta sheet
# and that is the location of the csv

# Error message if the files aren't provideds
if [ -z "$1" ]; then 
  echo ERROR: META FILE WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target meta file containing all the csvs to pull from
  exit 1
fi

folder=$(cd "$(dirname "$0")";pwd)  # Save current folder as a variable


# Setup - generate required folders
mkdir -p PBS
mkdir -p log

libraryText=$1

while IFS=$'\t' read -r -a varArray; do
	csv=${varArray[0]}  # The control group
	IFS=$OLDIFS

	name=$(basename "${csv}" ".csv")
	mkdir -p barcharts/$name

	cat >${folder}/PBS/${name}_barplots'.R' <<EOF
library(tidyverse)
library(dplyr)
library(janitor)
library(biomaRt)

genes <- read.csv("${csv}", header = TRUE) 

stem <- c("LEFTY1", "LGR5", "RGMB", "ASCL2", "SMOC2", "OLFM4")
enterocyte <- c("FABP1", "IL32", "CEACAM5", "CEACAM7", "CEACAM1")
goblet <- c("MUC2", "SPINK4", "ZG16", "REG4", "FCGBP", "TFF1", "TFF3")
tuft <- c("LRMP", "SH2D6", "HCK", "PTGS1", "HPGDS", "AZGP1", "BMX")
enteroendocrine <- c("CHGA", "NEUROD1", "GCG", "TPH1", "PCSK1N", "SCGN")
best4 <- c("CA7", "BEST4", "SPIBM", "HES4", "LYPD8")
cTA <- c("TUBB", "STMN1", "CKS1B", "CENPM", "CENPMW", "HMGB2")
enterocytes_progen <- c("FABP1", "CA1", "SLC26A3", "CA2", "KRT20")
im_enterocytes1 <- c("GUCA2A", "AQP8", "CLCA4", "ANPEP", "PRAP1", "CEACAM7")
im_enterocytes2 <- c("FABP1", "SLC51B", "CEACAM1", "CEACAM7", "CA4", "SLC26A3")
im_goblet <- c("MUC2", "SPINK4", "RETNLB", "SERPINA1", "CLCA1", "KLK1")
m <- c("SOX8", "CD74", "MSLN", "LTB", "TNFAIP2", "SPIB")
sec_ta <- c("RETNLB", "SERPINA1", "CLCA1", "ITLN1", "SLC12A2", "KLK1")
ta1 <- c("MIR3648", "CDC20B", "MIR663A", "C9orf38", "REG1A", "REG1B")
ta2 <- c("LEFTY1", "GOLM1", "LCN2", "APRT", "MLEC", "C1QBP")

overallList <- list(stem, enterocyte, goblet, tuft, enteroendocrine, best4, cTA, enterocytes_progen, 
	im_enterocytes1, im_enterocytes2, im_goblet, m, sec_ta, ta1, ta2)

names(overallList) <- c("Stem cells", "Enterocytes", "Goblet cells", "Tuft cells", 
	"Entereoendocrine cells", "Best4+ Enterocytes", "Cycling TA cells", "Enterocyte progenitors", 
	"Immature Enterocytes 1", "Immature Enterocytes 2", "Immature goblet cells", "M cells", 
	"Secretory TA cells", "TA1 cells", "TA2 cells")

for (i in 1:length(overallList)) {
	print(overallList)
	marks <- overallList[[i]]
	name <- names(overallList)[i]

	fMarks <- genes[genes\$gene %in% marks, ] %>% arrange(desc(log2FoldChange))
	bLvls <- unique(as.list(fMarks\$gene))

	ggplot(data=fMarks, aes(x=factor(gene, level = bLvls), y = log2FoldChange, fill = -log10(padj))) + 
		geom_bar(stat = "identity") + theme_bw() + ylab(expression("log"[2]*"foldChange")) + xlab("Gene") + 
		labs(fill = expression("-log"[10]*"pAdj")) + theme_bw() + 
		ggtitle(paste("${name}", name, "genes"))

	uName <- sub(" ", "_", name)
	ggsave(paste0("barcharts/$name/${name}_", uName, ".pdf"))
}
EOF

	cat >${folder}/PBS/${name}_barplots.pbs <<EOF
#!/bin/bash
#SBATCH --job-name=${name}_barplots
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${name}_barplots_%j.txt -e ${folder}/log/${name}_barplots_%j.err.txt
cd ${folder}
#------- End of header -------#
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate ChIPseeker

Rscript ${folder}/PBS/${name}_barplots.R
EOF
	sbatch ${folder}/PBS/${name}_barplots.pbs

done < ${libraryText}