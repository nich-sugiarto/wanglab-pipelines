suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(GetoptLong))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(NbClust))


args = commandArgs(trailingOnly=TRUE)
factor = args[1]
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
#Importing CSVs of Normalized Diffbind Scores
scores1=as_tibble(read.csv(list.files("diffBind/Combined/",pattern = "Scores",full.names = T), sep = "", header = T))

#Store Chromosomal positions
positions <- as_tibble(cbind(scores1[,c(1:3)]))
positions <- unite(positions[,c(1:3)], index, sep="-")
positions <- as.matrix(positions)

scores <- cbind(scores1[,c(4:11)])

#Average and Make new matrix for EdgeR and DESeq
UT <- as.matrix(select(scores,matches('UT')) %>% rowMeans)
H4 <- as.matrix(select(scores,matches('_4h')) %>% rowMeans)
H24 <- as.matrix(select(scores,matches('_24h')) %>% rowMeans)
H72 <- as.matrix(select(scores,matches('_72h')) %>% rowMeans)
edger <- cbind(UT,H4,H24,H72)
colnames(edger) <- c("UT","H4","H24","H72")

#Indexing rownames with index
row.names(edger) <- positions[,1]

#Generating Matrix containing L2FC data, 
#Note that there are actually quite a few "0" entries resulting in NA / Inf after L2FC, this simply removes the entire row... 
#Not sure if this is the best solution, but it's all Ive got at the moment.
edger_L2FC <- cbind(log2(edger[,c(2:4)]/edger[,1]))
edger_L2FC <- edger_L2FC[is.finite(rowSums(edger_L2FC[,c(1:3)])),]
edger_L2FC <- as.matrix(edger_L2FC)

#Generating Matrix Containg L10 Transformed Scores
edger <- cbind(log10(edger[,c(1:4)]))
edger <- edger[is.finite(rowSums(edger[,c(1:4)])),]
edger <- as.matrix(edger)

a <- fviz_nbclust(edger_L2FC, kmeans, method = "silhouette", k.max = 24) + theme_minimal() + ggtitle("Silhouette Plot")
b <- a$data
#Get num clusters for max Silhouette score
clusters <- as.numeric(rownames(b)[which(b[,2] == max(b[,2]), arr.ind = TRUE)])

#Setting up colors for L2FC heatmaps
col_fun = diverging_hcl(n = 20, h = c(255, 12), c = c(60, 100), l = c(15, 95), power = c(0.5, 1))
#Plotting L2FC Heatmaps
d <- dist(edger_L2FC[,c(1:3)], method = "maximum") # distance matrix
a <- fviz_nbclust(d, kmeans, method = "silhouette", k.max = 24) + theme_minimal() + ggtitle("Silhouette Plot")

fit <- hclust(d, method = "ward.D2")
res.nbclust <- NbClust(fit, distance = "euclidean",
                       min.nc = 2, max.nc = 9, 
                       method = "complete", index ="all")

saveFile = paste0("diffBind_Heatmaps/",factor,"_L2FC.pdf")
name = paste0(factor,"_L2FC")

clusters = 
heatmap=pheatmap(edger_L2FC[,c(1:3)],
    cluster_rows = fit,
    cluster_cols = F,
    cutree_rows = clusters,
    show_rownames=F,
    main = name,
    color = col_fun)
pdf(file=saveFile)
draw(heatmap)
dev.off()

#Writes Clusters to Bed
for (i in 1:clusters){
  cluster <- data.frame(row.names(edger_L2FC[unlist(row_order(heatmap)[i]),]))
  colnames(cluster) <- "index"
  cluster$chr <- unlist(strsplit(as.character(cluster$index),"-"))[c(T,F,F)]
  cluster$pos1 <- unlist(strsplit(as.character(cluster$index),"-"))[c(F,T,F)]
  cluster$pos2 <- unlist(strsplit(as.character(cluster$index),"-"))[c(F,F,T)]
  cluster_bed <- cluster[,c("chr","pos1","pos2")]
  #Write Cluster to Bed
  cluster_name <- paste0("diffBind_Heatmaps/Clusters/",factor,"_cluster",i,".bed")
  write.table(cluster_bed,cluster_name,quote=F,sep = "\t",row.names = F, col.names = F)
}



#Setting up colors for scores heatmaps
col_fun = sequential_hcl(n = 40, h = c(300, 75), c = c(40, NA, 95), l = c(10, 90), power = c(1.75, 0.75))

# #Generating Scores Heatmaps
# d <- dist(edger[,c(1:4)], method = "maximum") # distance matrix
# fit <- hclust(d, method = "ward.D2")
# saveFile = paste0("diffBind_Heatmaps/",factor,"_Scores.pdf")
# name = paste0(factor,"_Score")
# clusters = length(fit)-3
# heatmap=pheatmap(edger[,c(1:4)],
#                  cluster_rows = fit,
#                  cutree_rows = clusters,
#                  cluster_cols = F,
#                  show_rownames=F,
#                  main = name,
#                  color = col_fun)
# pdf(file=saveFile)
# draw(heatmap)
# dev.off()
