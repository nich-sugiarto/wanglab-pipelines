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
library(pheatmap)
library(RColorBrewer)
library(wesanderson)

options(ggrepel.max.overlaps = Inf)

# Read in and prepare inputs 
gTable <- read.csv("${geneTable}", header = TRUE)  # Read in gene table

mList <- read.csv("${markGenes}", header = FALSE)  # Read in the genes that will be annotated in various visualizations

mList <- mList\$V1  # Turn it into a vector

# Create a column to denote which genes are upregulated, downregulated, or not significantly changed
gTable <- gTable %>% 
	mutate( threshold_KD = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "Upregulated",
	log2FoldChange <= -1 & padj <= 0.05 ~ "Downregulated",
	TRUE ~ "Unchanged"))

# Volcano Plots
## Grab the top 25 upregulated genes, and the top 25 downregulated genes by padj value
## FIXME: For some reason, it seems like more significant genes are being skipped... something to look into?

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

## Create the volcano plot w/o any labels
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

## Create the labels for the genes present in the list
list <- gTable[gTable\$gene %in% mList, ]
p3 <-  p2 +
	geom_label_repel(data = list,
	mapping = aes(log2FoldChange, -log10(padj), 
		label = gene), min.segment.length = 0.0000001, size = 2)

pdf("pVolcanos/${cGene}_${cMark}_volcano.pdf")
	print(p3)
dev.off()

## Create the labels for top 50 genes (25 up/down ea.)
p3 <-  p2 +
	geom_label_repel(data = top_genes_p,
	mapping = aes(log2FoldChange, -log10(padj), 
		label = gene), min.segment.length = 0.0000001, size = 2)

print(p3)
ggsave("pVolcanos/${cGene}_volcano_t50.pdf")

# MA Plots
## Create base graph w/o any labels
p2 <- ggplot(gTable, aes(log2(baseMean), log2FoldChange)) +
		geom_point(aes(color = threshold_KD), size = 2/5) +
		xlab(expression("log"[2]*"baseMean")) +  
		ylab(expression("log"[2]*"foldChange")) +
		scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
		theme_bw() + theme(legend.position = "none")
print(p2)

ggsave("pVolcanos/${cGene}_MA.pdf")

## Add labels for the genes present in the provided list
p3 <-  p2 +
	geom_label_repel(data = list,
	mapping = aes(log2(baseMean), log2FoldChange, 
	label = gene), min.segment.length = 0.0000001, size = 2)

print(p3)
ggsave("pVolcanos/${cGene}_${cMark}_MA.pdf")

## Add labels for the top 50 genes (25 up/down ea.)
p3 <-  p2 +
	geom_label_repel(data = top_genes_p,
	mapping = aes(log2(baseMean), log2FoldChange, 
	label = gene), min.segment.length = 0.0000001, size = 2)
print(p3)
ggsave("pVolcanos/${cGene}_MA_t50.pdf")

# Heatmaps
## FIXME: Currently makes a couple pretty big assumptions.
## The first is that it assumes that the normalized_counts file is named very similarly to the provided .csv
## Where it's the exact same, but _all.csv is replaced by the appropriate suffix

normCounts <- read.table("${geneTable/_all.csv/_normalized_counts.txt}", header = TRUE, row.names = 1)

## Get 200 most significant genes by p value
t200 <- bind_rows(
  gTable %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(200))

## Turn into dataframe of top 200 genes by p value
hmap <- normCounts %>% 
	dplyr::filter(row.names(normCounts) %in% t200\$gene) %>%
	data.frame()

annoLbls <- row.names(hmap)
annoLbls[!(annoLbls %in% mList)] <- ""

topLbls <- row.names(hmap)
topLbls[!(topLbls %in% top_genes_p\$gene)] <- ""

lfc <- filter(gTable, abs(log2FoldChange) >= 1.5, padj <= 0.05,na.rm = TRUE) %>% 
	data.frame()

normMap <- normCounts %>% 
	data.frame()
str(normMap)

normMap <- normMap[row.names(normMap) %in% lfc\$gene, ]
str(normMap)

pheatmap(normMap,
	cluster_rows = T,
	cluster_cols = T,
	color = colorRampPalette(rev(brewer.pal(n = 7, name =
  		"RdBu")))(100),
	show_rownames = F,
	border_color = "white",
	fontsize = 4,
	scale = "row",
	treeheight_row = 0,
	fontsize_row = 5,
	fontsize_col = 5,
	cellwidth = 10, 
	filename = "pVolcanos/${cGene}_heatmap.pdf")

pheatmap(hmap,
	cluster_rows = T,
	cluster_cols = T,
	color = colorRampPalette(rev(brewer.pal(n = 7, name =
  		"RdBu")))(100),
	show_rownames = T,
	border_color = "white",
	fontsize = 4,
	labels_row = annoLbls,
	scale = "row",
	fontsize_row = 5,
	fontsize_col = 5,
	cellwidth = 10, 
	filename = "pVolcanos/${cGene}_${cMark}_heatmap.pdf")

pheatmap(hmap,
	cluster_rows = T,
	cluster_cols = T,
	color = colorRampPalette(rev(brewer.pal(n = 7, name =
  		"RdBu")))(100),
	show_rownames = T,
	border_color = "white",
	fontsize = 4,
	labels_row = topLbls,
	scale = "row",
	fontsize_row = 5,
	fontsize_col = 5,
	cellwidth = 10, 
	filename = "pVolcanos/${cGene}_t50_heatmap.pdf")


pal <- wes_palette("Zissou1", 100, type = "continuous")

pheatmap(normMap,
	cluster_rows = T,
	cluster_cols = T,
	color = pal,
	show_rownames = F,
	border_color = "white",
	fontsize = 4,
	scale = "row",
	treeheight_row = 0,
	fontsize_row = 5,
	fontsize_col = 5,
	cellwidth = 10, 
	filename = "pVolcanos/${cGene}_heatmap_wesAnd.pdf")

pheatmap(hmap,
	cluster_rows = T,
	cluster_cols = T,
	color = pal,
	show_rownames = T,
	border_color = "white",
	fontsize = 4,
	labels_row = annoLbls,
	scale = "row",
	fontsize_row = 5,
	fontsize_col = 5,
	cellwidth = 10, 
	filename = "pVolcanos/${cGene}_${cMark}_heatmap_wesAnd.pdf")

pheatmap(hmap,
	cluster_rows = T,
	cluster_cols = T,
	color = pal,
	show_rownames = T,
	border_color = "white",
	fontsize = 4,
	labels_row = topLbls,
	scale = "row",
	fontsize_row = 5,
	fontsize_col = 5,
	cellwidth = 10, 
	filename = "pVolcanos/${cGene}_t50_heatmap_wesAnd.pdf")

EOF

	cat >${folder}/PBS/${cGene}_${cMark}'_volcano.PBS' <<EOF
#!/bin/bash
#SBATCH --job-name=${cGene}_${cMark}_volcano
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
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