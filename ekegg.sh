if [ -z "$1" ]; then 
  echo ERROR: TARGET FOLDER WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target folder
  exit 1
fi

folder=$(cd "$(dirname "$0")";pwd)

for file in $1/*.bed; do
  base=$(basename "$file" ".bed")
  mkdir -p ${folder}/ChIPseeker/${base}
  cat > ${folder}/PBS/${base}_ChIPseeker.r <<EOF

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(annotables)
library(org.Hs.eg.db)
library(DOSE)
library(dplyr)
library(R.utils)

samplefiles <- list.files("$folder/$1/", pattern= "${base}.bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("${base}")

# Assign annotation db
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-2000, 2000), verbose=FALSE)

annot <- as.data.frame(peakAnnoList[["${base}"]]@anno)
# Get unique entrez gene Ids
entrezids <- unique(annot\$geneId)
# Get hg38 entrez to gene symbol mappings

entrez2gene <- grch38 %>% filter(entrez %in% entrezids) %>% dplyr::select(entrez, symbol)

# Match to each annotation dataframe
m <- match(annot\$geneId, entrez2gene\$entrez)
dim(annot)
annot <- cbind(annot[,1:ncol(annot)], geneSymbol=entrez2gene\$symbol[m])

entrezids <- annot\$geneId %>%
  as.character() %>%
  unique()

ego <- enrichGO(gene = entrezids,
                keyType = "ENTREZID",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

egoSim <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
cluster_summary <- data.frame(ego)
cluster_summary_sim <- data.frame(egoSim)

pdf(file = "$folder/ChIPseeker/${base}/${base}_GO.pdf")
dotplot(ego,showCategory = 20,font.size=6)
dev.off()

R.utils::setOption("clusterProfiler.download.method","wget")
ekegg <- enrichKEGG(gene = entrezids,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
pdf(file = "$folder/ChIPseeker/${base}/${base}_KEGG.pdf")
dotplot(ekegg,showCategory = 20,font.size=6)
dev.off()

do = enrichDO(entrezids)
pdf(file = "$folder/ChIPseeker/${base}/${base}_DO.pdf")
dotplot(do, showCategory=20,font.size=6)
dev.off()

R.utils::setOption("clusterProfiler.download.method","wget")
ekegg <- enrichKEGG(gene = entrezids,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
pdf(file = "$folder/ChIPseeker/${base}/${base}_KEGG.pdf")
dotplot(ekegg,showCategory = 20,font.size=6)
dev.off()

do = enrichDO(entrezids)
pdf(file = "$folder/ChIPseeker/${base}/${base}_DO.pdf")
dotplot(do, showCategory=20,font.size=6)
dev.off()

EOF

Rscript ${folder}/PBS/${base}_ChIPseeker.r

done