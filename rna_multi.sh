#!/bin/bash

# Takes all the rnaseq data from the user's 'fastq' file and
# visualizes them all in relation to each other. Adapted from the previous iteration called
# rnaseq_multi.sh

# Requires the deseq environment

# Version Doc: https://docs.google.com/document/d/1lq1KnU1qZx2OIL3b9nn_3gTlF2xTfb9kaZwjQZ46PI4/edit#

mkdir -p multideseq
folder=$(cd "$(dirname "$0")";pwd)

cat >${folder}/multideseq/multideseq'.R' <<EOF
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(cowplot)
library(clusterProfiler)
library(DEGreport)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(tximport)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(ensembldb)
library(SummarizedExperiment)
EOF

cd counts

# Create the meta text

cat > ${folder}/multideseq/meta.txt<<EOF
    sampletype
EOF

i=0
for file in *_featurecounts_Count.txt
do
	echo "${file%_feature*}  ${file%_R*}" >> ${folder}/multideseq/meta.txt
	echo ${file} added
	cat >>${folder}/multideseq/multideseq'.R' <<EOF

data${i} <- read.table("../counts/${file}",header = TRUE,skip=1)
names(data${i})[names(data${i}) == "aligned.${file%_feature*}_sorted.bam"] <- "${file%_feature*}"
EOF
dcount+="data${i},"
((i++))
done

cat >>${folder}/multideseq/multideseq'.R' <<EOF
name1 <- read.table("../counts/${file%_featurecounts_Count.txt}_featurecounts_Name.txt",header = TRUE,skip=1)
data <- data.frame(${dcount}row.names = name1\$Geneid)
meta <- read.table("meta.txt", header=T, row.names=1)

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

rld <- rlog(dds, blind=TRUE)
expmatrix <- SummarizedExperiment::assay(rld)
plotPCA(rld, intgroup="sampletype")
ggsave("PCA.png",dpi=300)

rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

xx <- pheatmap(rld_cor)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x\$gtable)
  dev.off()
}

save_pheatmap_pdf(xx, "Hierarchical_Clustering.pdf")

# Commented out because this is from a time that I was young and stupid
# # Not sure why I had to recreate the dds here, but it only worked this way. Otherwise I got a weird error
# dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
dds <- dds[rowSums(DESeq2::counts(dds)) > 1, ]  # Removed empty rows to make it faster
data_DESeq <- DESeq2::DESeq(dds)

# Taken from https://shiring.github.io/rna-seq/deseq2/teaching/2016/09/29/DESeq2-course

expmatrix_DESeq <- DESeq2::rlog(data_DESeq, fitType="local")
expmatrix <- SummarizedExperiment::assay(expmatrix_DESeq)

select <- order(rowMeans(expmatrix),decreasing=TRUE)[1:50]

df <- data.frame(sampletype = SummarizedExperiment::colData(data_DESeq)[,c("sampletype")], row.names = rownames(SummarizedExperiment::colData(data_DESeq)))
xx <- pheatmap::pheatmap(expmatrix[select,], 
	cluster_rows=TRUE,
	show_rownames=TRUE,
	cluster_cols=FALSE,
	annotation_col=df,
	border_color = NA,
	fontsize = 5,
	scale = "row",
	fontsize_row = 5,
	fontsize_col = 5,
	cellwidth = 10)

save_pheatmap_pdf(xx, "50siggenes.pdf")

xx <- pheatmap::pheatmap(expmatrix, 
	cluster_rows=TRUE,
	show_rownames=FALSE,
	cluster_cols=FALSE,
	annotation_col=df,
	border_color = NA,
	fontsize = 5,
	treeheight_row = 0,
	treeheight_col = 0,
	scale = "row",
	fontsize_row = 5,
	fontsize_col = 5,
	cellwidth = 10)

save_pheatmap_pdf(xx, "heatmap.pdf")
write.table(expmatrix, file="heatmap.csv",row.names = TRUE,quote = FALSE,sep = ",")
EOF

cat >${folder}/PBS/'multideseq.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=multideseq # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=16GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=10:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
################################
# Enter your code to run below #
################################
source activate deseq

cd ${folder}/multideseq
Rscript multideseq.R
EOF
cd ${folder}/log
sbatch ${folder}/PBS/multideseq.pbs