#!/bin/bash/

# Takes in a csv of metascape output (formatted the same way as the one that Luke gave me)
# Makes a corresponding GO plot of the data, and for each entry description, plots the log fold change and 
# adjusted p value.

# USAGE:
# This pipeline takes in one positional argument:
# 	$1 - target meta file containing sample information
# Where the meta file is tab-delimited, and formated as follows:
# locationOfMetaCSV	DESeqResultTable

# Error message if the files aren't provideds
if [ -z "$1" ]; then 
  echo ERROR: META FILE WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target meta file containing sample information
  echo Where the meta file is tab-delimited, and formated as follows:
  echo locationOfMetaCSV	DESeqResultTable
  exit 1
fi

mkdir -p goPlots
mkdir -p PBS
mkdir -p log

folder=$(pwd)  # Save current folder as a variable

libraryText=$1

while IFS=$'\t' read -r -a varArray; do
	metascape=${varArray[0]}  # Location for metascape table
	de=${varArray[1]}

	mName=$(basename "$metascape" ".csv")

	cat >${folder}/PBS/${mName}'_GO.R' <<EOF
library(ggplot2)
library(dplyr)
library(tidyverse)
library(janitor)
library(viridis)

gTable <- read.csv("$de", header = TRUE)

go <- read.csv("$metascape", header = TRUE) %>%
	clean_names() %>%
	unique() %>%
	arrange(x_in_go)

goLvls <- as.list(go\$description)

# Modify upper bound of xlim as needed
goPlot <- ggplot(data = go , aes(x = x_in_go , y =factor(description, level = unique(goLvls)), 
	color = log_p, size = x_gene_in_go_and_hit_list)) + 
	geom_point() + theme_bw() + ylab("") + xlab("% in GO") + 
	labs(color = "pAdj", size = "# Genes in GO and Hit List") + xlim(0, 8) +
	scale_color_viridis(option="viridis", direction = 1) + 
	ggtitle("GO enrichment analysis")

print(goPlot)
ggsave("goPlots/${mName}_GO.pdf")
ggsave("goPlots/${mName}_BiggerGO.pdf", width = 12, height = 7)

for (row in 1:nrow(go)) {
	hits <- as.list(strsplit(go[row, 21], "\\\\|"))  # Ask Nick if this doesn't make sense to you
	
	trunc <- gTable[gTable\$gene %in% hits[[1]], ] %>% arrange(desc(log2FoldChange))

	bLvls <- unique(as.list(trunc\$gene))
	barPlot <- ggplot(data=trunc, aes(x=factor(gene, level = bLvls), y = log2FoldChange, fill = -log10(padj))) + 
		geom_bar(stat = "identity") + theme_bw() + ylab(expression("log"[2]*"foldChange")) + xlab("Gene") + 
		labs(fill = expression("-log"[10]*"(pAdj)")) + theme_bw() + 
		theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1)) + 
		ggtitle(paste(go[row, 9], "gene hits"))

	name <- sub(" ", "_", go[row, 9])
	pw <- nrow(trunc)*0.5
	ggsave(paste0("goPlots/", name, ".pdf"), height = 7, width = pw, limitsize = FALSE)
}
EOF

	cat >${folder}/PBS/${nName}'_GO.pbs' <<EOF
#!/bin/bash
#SBATCH --job-name=${mName}_GO
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${mName}_GO_%j.txt -e ${folder}/log/${mName}_GO_%j.err.txt
cd ${folder}
#------- End of header -------#
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate ChIPseeker

Rscript ${folder}/PBS/${mName}_GO.R
EOF

	sbatch PBS/${nName}'_GO.pbs'
done < ${libraryText}