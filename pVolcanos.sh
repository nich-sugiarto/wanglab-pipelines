#!/bin/bash/

# Reads in the csvs from two factor deseq comparison, as well as a list of genes. 
# Creates two volcano plots in the nicer aesthetic: 
# One that has no genes labelled, and one that has the selected genes labelled

# USAGE:
# This pipeline takes in one positional argument:
# 	$1 - target meta file containing sample information
# Where the meta file is tab-delimited, and formated as follows:
# locationOfDifferentialGeneCSV	listOfGenesToMark

# Note that the list of genes should NOT have a header

# Error message if the files aren't provideds
if [ -z "$1" ]; then 
  echo ERROR: META FILE WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target meta file containing sample information
  echo Where the meta file is tab-delimited, and formated as follows:
  echo locationOfDifferentialGeneCSV	listOfGenesToMark
  exit 1
fi

libraryText=$1

while IFS=$'\t' read -r -a varArray; do
	geneTable=${varArray[0]}  # Location for csv of all genes
	markGenes=${varArray[1]}  # Location for csv for genes to mark
	IFS=$OLDIFS

	cGene=$(basename "${geneTable}" "_all.csv")
	cMark=$(basename "${markGenes}" ".csv")

	folder=$(cd "$(dirname "$0")";pwd)  # Save current folder as a variable
	cat >${folder}/PBS/${cGene}_${cMark}'_volcano.R' <<EOF
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(knitr)

gTable <- read.csv("${geneTable}", header = TRUE)

mList <- read.csv("${markGenes}", header = FALSE)
str(mList)

mList <- mList\$V1

gTable <- gTable %>% 
	mutate( threshold_KD = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "Upregulated",
	log2FoldChange <= -1 & padj <= 0.05 ~ "Downregulated",
	TRUE ~ "Unchanged"))

top <- 25
top_genes_p <- bind_rows(
  gTable %>% 
    filter(threshold_KD == 'Upregulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top),
  gTable %>% 
    filter(threshold_KD == 'Downregulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top)
)

p2 <- ggplot(gTable, aes(log2FoldChange, -log10(padj))) +
		geom_point(aes(color = threshold_KD), size = 2/5) +
		xlab(expression("log"[2]*"FC")) +  
		ylab(expression("-log"[10]*"pAdj")) +
		scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
		guides(colour = guide_legend(override.aes = list(size=1.5))) +
		theme_bw() + theme(legend.position = "none")

pdf("pVolcanos/${cGene}_volcano.pdf")
	print(p2)
dev.off()

list <- gTable[gTable\$gene %in% mList, ]
p3 <-  p2 +
	geom_label_repel(data = list,
	mapping = aes(log2FoldChange, -log10(padj), 
		label = gene), min.segment.length = 0.0000001, size = 2)

pdf("pVolcanos/${cGene}_${cMark}_volcano.pdf")
	print(p3)
dev.off()

p3 <-  p2 +
	geom_label_repel(data = top_genes_p,
	mapping = aes(log2FoldChange, -log10(padj), 
		label = gene), min.segment.length = 0.0000001, size = 2)

print(p3)
ggsave("pVolcanos/${cGene}_volcano_t50.pdf")

p2 <- ggplot(gTable, aes(log2(baseMean), log2FoldChange)) +
		geom_point(aes(color = threshold_KD), size = 2/5) +
		xlab(expression("log"[2]*"baseMean")) +  
		ylab(expression("log"[2]*"foldChange")) +
		scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
		theme_bw() + theme(legend.position = "none")
print(p2)

ggsave("pVolcanos/${cGene}_MA.pdf")

p3 <-  p2 +
	geom_label_repel(data = list,
	mapping = aes(log2(baseMean), log2FoldChange, 
	label = gene), min.segment.length = 0.0000001, size = 2)

print(p3)
ggsave("pVolcanos/${cGene}_${cMark}_MA.pdf")


p3 <-  p2 +
	geom_label_repel(data = top_genes_p,
	mapping = aes(log2(baseMean), log2FoldChange, 
	label = gene), min.segment.length = 0.0000001, size = 2)
print(p3)
ggsave("pVolcanos/${cGene}_MA_t50.pdf")
EOF

	cat >${folder}/PBS/${cGene}_${cMark}'_volcano.PBS' <<EOF
#!/bin/bash
#SBATCH --job-name=${cGene}_${cMark}_volcano
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${cGene}_${cMark}_volcano_%j.txt -e ${folder}/log/${cGene}_${cMark}_volcano_%j.err.txt
cd ${folder}
#------- End of header -------#
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate deseq

Rscript ${folder}/PBS/${cGene}_${cMark}_volcano.R
EOF

	sbatch ${folder}/PBS/${cGene}_${cMark}'_volcano.PBS'
done < ${libraryText}