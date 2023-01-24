#!/bin/bash

# Compares RNA-seq results to chosen control to identify differntially expressed genes.
# Compares TWO conditions, each with TWO technical replicates

# Requires the deseq environment

# The meta text should be formatted as follows:
# controlGroup	sampleGroup	repDelimiter
# Where the rep delimiter is how reps are distinguised (e.g. R1 and R2 vs Rep1 and Rep2)

# Setup - Generate required folders
mkdir -p PBS
mkdir -p log
mkdir -p twoFactorComparison

echo WARNING:
echo I modified some things recently, and I have no idea if they work.
echo "If you don't see MA plots or volcano plots, please let Nick know"

folder=$(cd "$(dirname "$0")";pwd)  # Save current folder location

if [ -z "$1" ]; then 
  echo ERROR: TARGET FOLDER WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target meta file containing sample information
  echo Where the file is tab delimited and contains the following fields:
  echo controlGroup	sampleGroup	repDelimiter
  exit 1
fi
libraryText=$1

OLDIFS=$IFS

while IFS=$'\t' read -r -a varArray; do
	control=${varArray[0]}  # The control group
	treatment=${varArray[1]}  # The treatment group
	delim=${varArray[2]}  # How reps are distinguised
	IFS=$OLDIFS

	mkdir -p ${folder}/twoFactorComparison/${control}_over_${treatment}_deseq
	cat > ${folder}/twoFactorComparison/${control}_over_${treatment}_deseq/meta.txt<<EOF
    sampletype
${control}_${delim}1  ${control}
${control}_${delim}2  ${control}
${treatment}_${delim}1  ${treatment}
${treatment}_${delim}2  ${treatment}
EOF

cat >${folder}/twoFactorComparison/${control}_over_${treatment}_deseq/${treatment}_over_${control}_deseq'.R' <<EOF

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
  
  data1 <- read.table("../../counts/${control}_${delim}1_featurecounts_Count.txt",header = TRUE,skip=1)
  names(data1)[names(data1) == "aligned.${control}_${delim}1_sorted.bam"] <- "${control}_${delim}1"
  data2 <- read.table("../../counts/${control}_${delim}2_featurecounts_Count.txt",header = TRUE,skip=1)
  names(data2)[names(data2) == "aligned.${control}_${delim}2_sorted.bam"] <- "${control}_${delim}2"

  data3 <- read.table("../../counts/${treatment}_${delim}1_featurecounts_Count.txt",header = TRUE,skip=1)
  names(data3)[names(data3) == "aligned.${treatment}_${delim}1_sorted.bam"] <- "${treatment}_${delim}1"
  data4 <- read.table("../../counts/${treatment}_${delim}2_featurecounts_Count.txt",header = TRUE,skip=1)
  names(data4)[names(data4) == "aligned.${treatment}_${delim}2_sorted.bam"] <- "${treatment}_${delim}2"
  
  name1 <- read.table("../../counts/${treatment}_${delim}1_featurecounts_Name.txt",header = TRUE,skip=1)

  data <- data.frame(data1,data2,data3,data4,row.names = name1\$Geneid)
  meta <- read.table("meta.txt", header=T, row.names=1)
  dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

  dds <- estimateSizeFactors(dds)
  normalized_counts <- counts(dds, normalized=TRUE)
  write.table(normalized_counts, file="${control}_${treatment}_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

  rld <- rlog(dds, blind=TRUE)
  plotPCA(rld, intgroup="sampletype") + ylim(-10, 10)
  ggsave("${control}_${treatment}_PCA.png",dpi=300)

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
  save_pheatmap_pdf(xx, "${control}_${treatment}_Hierarchical_Clustering.pdf")


  dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
  dds <- DESeq(dds)
  plotDispEsts(dds)
  dev.copy(png,'${control}_${treatment}_dispersion.png')
  dev.off()


  contrast_kd <-  c("sampletype", "${control}","${treatment}")
  res_tableKD_unshrunken <- results(dds, contrast=contrast_kd, alpha = 0.05)
  res_tableKD <- lfcShrink(dds, contrast=contrast_kd, res=res_tableKD_unshrunken,type = "ashr")
  plotMA(res_tableKD_unshrunken, ylim=c(-2,2))
  dev.copy(png,'${control}_${treatment}_unshrunken_dispersion.png')
  dev.off()
  plotMA(res_tableKD, ylim=c(-2,2))
  dev.copy(png,'${control}_${treatment}_shrunken_dispersion.png')
  dev.off()
  sink("${control}_${treatment}_summary.txt")
  summary(res_tableKD)
  sink()


  padj.cutoff <- 0.05
  lfc.cutoff <- 0.58


  res_tableKD_tb <- res_tableKD %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()

  sigKD <- res_tableKD_tb %>%
    dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
  write.table(sigKD,file="${control}_${treatment}_fc1.5.csv",row.names = FALSE,quote = FALSE,sep = ",")
  
  total <- res_tableKD_tb %>%
    dplyr::filter(padj < padj.cutoff)
  write.table(total,file="${control}_${treatment}.csv",row.names = FALSE,quote = FALSE,sep = ",")

  up <- res_tableKD_tb %>%
    dplyr::filter(padj < padj.cutoff & log2FoldChange > 0)
  write.table(up,file="${control}_${treatment}_up.csv",row.names = FALSE,quote = FALSE,sep = ",")

  down <- res_tableKD_tb %>%
    dplyr::filter(padj < padj.cutoff & log2FoldChange < 0)
  write.table(down,file="${control}_${treatment}_down.csv",row.names = FALSE,quote = FALSE,sep = ",")
  
  up15 <- res_tableKD_tb %>%
    dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff & log2FoldChange > 0)
  write.table(up15,file="${control}_${treatment}_fc1.5_up.csv",row.names = FALSE,quote = FALSE,sep = ",")
  
  down15 <- res_tableKD_tb %>%
    dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff  & log2FoldChange < 0)
  write.table(down15,file="${control}_${treatment}_fc1.5_down.csv",row.names = FALSE,quote = FALSE,sep = ",")

  lfc2 <- res_tableKD_tb %>%
    dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > 1)
  write.table(lfc2,file="${control}_${treatment}_fc2.csv",row.names = FALSE,quote = FALSE,sep = ",")
 
  up2 <- res_tableKD_tb %>%
    dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > 1 & log2FoldChange > 0)
  write.table(up2,file="${control}_${treatment}_fc2_up.csv",row.names = FALSE,quote = FALSE,sep = ",")
  
  down2 <- res_tableKD_tb %>%
    dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > 1  & log2FoldChange < 0)
  write.table(down2,file="${control}_${treatment}_fc2_down.csv",row.names = FALSE,quote = FALSE,sep = ",")

  write.table(res_tableKD,file="${control}_${treatment}_all.csv",row.names = FALSE,quote = FALSE,sep = ",")

meta <- meta %>%
rownames_to_column(var="samplename") %>%
as_tibble()

normalized_counts <- normalized_counts %>%
data.frame() %>%
rownames_to_column(var="gene") %>%
as_tibble()



norm_sig <- normalized_counts %>%
dplyr::filter(gene %in% sigKD\$gene) %>%
data.frame() %>%
column_to_rownames(var = "gene")

annotation <- meta %>%
  dplyr::select(samplename, sampletype) %>%
  data.frame(row.names = "samplename")

xx <- pheatmap(norm_sig,
cluster_rows = T,
show_rownames = F,
annotation = annotation,
border_color = NA,
fontsize = 5,
scale = "row",
fontsize_row = 5,
fontsize_col = 5,
cellwidth = 10)
save_pheatmap_pdf(xx, "${control}_${treatment}_fc1.5_Heatmap.pdf")
write.table(norm_sig,file="${control}_${treatment}_fc1.5_heatmap.csv",row.names = TRUE,quote = FALSE,sep = ",")
         

norm_fc2 <- normalized_counts %>%
  dplyr::filter(gene %in% lfc2\$gene) %>%
  data.frame() %>%
  column_to_rownames(var = "gene")

xx <- pheatmap(norm_fc2,
cluster_rows = T,
show_rownames = F,
annotation = annotation,
border_color = NA,
fontsize = 5,
scale = "row",
fontsize_row = 5,
fontsize_col = 5,
cellwidth = 10)
save_pheatmap_pdf(xx, "${control}_${treatment}_fc2_Heatmap.pdf")
write.table(norm_fc2,file="${control}_${treatment}_fc2_heatmap.csv",row.names = TRUE,quote = FALSE,sep = ",")


top50_sig_genes <- res_tableKD_tb %>%
  arrange(padj) %>%     #Arrange rows by padj values
  pull(gene) %>%         #Extract character vector of ordered genes
  head(n=50)         #Extract the first 20 genes

top50_sig_norm <- normalized_counts %>%
  dplyr::filter(gene %in% top50_sig_genes) %>%
  data.frame() %>%
  column_to_rownames(var = "gene")
  
  write.table(top50_sig_norm,file="${control}_${treatment}_top50_sig.csv",row.names = TRUE,quote = FALSE,sep = ",")
  
xx <- pheatmap(top50_sig_norm,
               cluster_rows = T,
               cluster_cols = T,
               show_rownames = T,
               annotation = annotation,
               border_color = NA,
               fontsize = 5,
               scale = "row",
               fontsize_row = 5,
               fontsize_col = 5,
               cellwidth = 10)

save_pheatmap_pdf(xx, "${control}_${treatment}_Top50siggene.pdf")


res_tableKD_tb <- res_tableKD_tb %>% 
	mutate( threshold_KD = case_when(log2FoldChange >= 1.5 & padj <= 0.05 ~ "Upregulated",
	log2FoldChange <= -1.5 & padj <= 0.05 ~ "Downregulated",
	TRUE ~ "Unchanged"))
  

  top <- 25
  top_genes_p <- bind_rows(
    res_tableKD_tb %>% 
      filter(threshold_KD == 'Upregulated') %>% 
      arrange(padj, desc(abs(log2FoldChange))) %>% 
      head(top),
    gTable %>% 
      filter(threshold_KD == 'Downregulated') %>% 
      arrange(padj, desc(abs(log2FoldChange))) %>% 
      head(top)
  )

  ggplot(res_tableKD_tb) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_KD)) +
    ggtitle("${treatment}_${control}") +
    xlab(expression("log"[2]*"FC")) +  
		ylab(expression("-log"[10]*"pAdj")) +
		scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
    theme_bw + 
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) +
  	geom_label_repel(data = top_genes_p,
  	mapping = aes(log2FoldChange, -log10(padj), 
  	label = gene), min.segment.length = 0.0000001, size = 2)


  ggsave("${treatment}_${control}_fc1.5_VolcanoPlot.png",dpi=300) 
  
res_tableKD_tb <- res_tableKD_tb %>% 
	mutate( threshold_KD = case_when(log2FoldChange >= 2 & padj <= 0.05 ~ "Upregulated",
	log2FoldChange <= -2 & padj <= 0.05 ~ "Downregulated",
	TRUE ~ "Unchanged"))

  top_genes_p <- bind_rows(
    res_tableKD_tb %>% 
      filter(threshold_KD == 'Upregulated') %>% 
      arrange(padj, desc(abs(log2FoldChange))) %>% 
      head(top),
    gTable %>% 
      filter(threshold_KD == 'Downregulated') %>% 
      arrange(padj, desc(abs(log2FoldChange))) %>% 
      head(top)
  )
  
  ggplot(res_tableKD_tb) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_KD)) +
    ggtitle("${treatment}_${control}") +
    xlab(expression("log"[2]*"FC")) +  
		ylab(expression("-log"[10]*"pAdj")) +
		scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
    theme_bw + 
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))
    geom_label_repel(data = top_genes_p,
  	mapping = aes(log2FoldChange, -log10(padj), 
  	label = gene), min.segment.length = 0.0000001, size = 2)

  ggsave("${treatment}_${control}_fc2_VolcanoPlot.png",dpi=300)

  p2 <- ggplot(gTable, aes(log2(baseMean), log2FoldChange)) +
		geom_point(aes(color = threshold_KD), size = 2/5) +
		xlab(expression("log"[2]*"baseMean")) +  
		ylab(expression("log"[2]*"foldChange")) +
		scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
		theme_bw() + theme(legend.position = "none")
    print(p2)

  ggsave("pVolcanos/${treatment}_${control}_MA.pdf")
  
  p3 <-  p2 +
  	geom_label_repel(data = top_genes_p,
  	mapping = aes(log2(baseMean), log2FoldChange, 
  	label = gene), min.segment.length = 0.0000001, size = 2)
  print(p3)
  ggsave("pVolcanos/${treatment}_${control}_MA_t50.pdf")
EOF

cat >${folder}/PBS/${treatment}_over_${control}'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${treatment}_over_${control} # Name of the job

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
#SBATCH --output=${folder}/log/%x.%j.out
#SBATCH --error=${folder}/log/%x.%j.err
################################
# Enter your code to run below #
################################
source activate deseq

cd ${folder}/twoFactorComparison/${control}_over_${treatment}_deseq
Rscript ${treatment}_over_${control}_deseq'.R'
EOF

sbatch PBS/${treatment}_over_${control}'.pbs'

done < ${libraryText}