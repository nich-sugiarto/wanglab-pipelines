#!/bin/bash

folder=$(cd "$(dirname "$0")";pwd)

mkdir ChIPseeker

for file in epic2/*.bed; #########IgG???
do
base=$(basename "$file" ".bed")
if [[ $base != *_IgG* ]]; then
mkdir ${folder}/ChIPseeker/${base}
cat > ${folder}/PBS/${base}_ChIPseeker.r <<EOF
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(annotables)
library(org.Hs.eg.db)
library(DOSE)
library(dplyr)

samplefiles <- list.files("$folder/epic2/", pattern= "${base}.bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("${base}")

# Assign annotation db
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-2000, 2000), verbose=FALSE)

sink("$folder/ChIPseeker/${base}/${base}_annotationsummary.txt")
peakAnnoList
sink()

pdf(file = "$folder/ChIPseeker/${base}/${base}_AnnoPie.pdf")
plotAnnoPie(peakAnnoList[["${base}"]])
dev.off()

pdf(file = "$folder/ChIPseeker/${base}/${base}_AnnoBar.pdf")
plotAnnoBar(peakAnnoList[["${base}"]])
dev.off()

pdf(file = "$folder/ChIPseeker/${base}/${base}_relativetoTSS.pdf")
plotDistToTSS(peakAnnoList, title="Distribution of binding loci relative to TSS")
dev.off()

annot <- as.data.frame(peakAnnoList[["${base}"]]@anno)
# Get unique entrez gene Ids
entrezids <- unique(annot\$geneId)
# Get hg38 entrez to gene symbol mappings

entrez2gene <- grch38 %>% filter(entrez %in% entrezids) %>% dplyr::select(entrez, symbol)

# Match to each annotation dataframe
m <- match(annot\$geneId, entrez2gene\$entrez)
annot <- cbind(annot[,1:21], geneSymbol=entrez2gene\$symbol[m])
write.table(annot,file="$folder/ChIPseeker/${base}/${base}_annotation.txt", sep="\t", quote=F, row.names=F)

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
write.csv(cluster_summary, "$folder/ChIPseeker/${base}/${base}_clusterProfiler.csv")
write.csv(cluster_summary_sim, "$folder/ChIPseeker/${base}/${base}_clusterProfilerSim.csv")

pdf(file = "$folder/ChIPseeker/${base}/${base}_GO.pdf")
dotplot(ego,showCategory = 20,font.size=6)
dev.off()

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

cat >${folder}/PBS/${base}_ChIPseeker'.pbs' <<thisprogram
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${base}_ChIPSeeker # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=32GB

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
cd ${folder}
source activate choccy

Rscript ${folder}/PBS/${base}_ChIPseeker.r

thisprogram

cd ${folder}/log

sbatch ${folder}/PBS/${base}_ChIPseeker.pbs
fi
done

