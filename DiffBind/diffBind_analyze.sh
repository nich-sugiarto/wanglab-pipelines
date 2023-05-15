#!/bin/bash

# The log fold change will be calculated as log(samp/control), or log(samp) - log(control)

# Script that creates a diffBind comparison, and completes the standard analysis

# This pipeline requires four positional arguments:
# 	$1 - The DiffBind category on which to make the comparison
# 	$2 - First category for comparison (sample group)
# 	$3 - Second category for comparison (control group)
# 	$4 - (Optional) name of the object from diffBind_createObject.sh to read in (default is 'dbObj.rds')

# NOTE: This assumes that there are only two reps. If there are more than two reps, then you must manually change that parameter within the script

# Requires the ATACQC pipeline (will rename in the future)

compare=$1
samp=$2
control=$3
dbObj=$4

caps=$(echo "$compare" | tr '[:lower:]' '[:upper:]')

#  All the parameters are not provided
if [ -z "$3" ]; then 
  echo ERROR: NOT ALL PARAMETERS HAVE BEEN SPECIFIED
  echo USAGE:
  echo This pipeline takes in a minimum of three positional argument:
  echo '$1 - The DiffBind category on which to make the comparison'
  echo '$2 - First category for comparison (sample group)'
  echo '$3 - Second category for comparison (control group)'
  echo "\$4 - (Optional) name of the object from diffBind_createObject.sh to read in (default is 'dbObj.rds')"
  exit 1
fi

if [ -z "$4" ]; then 
  dbObj="dbObj"
fi

off=FALSE
normType=${dbObj#*_}
if [ "$normType" == Default ]; then
    norm=""
  elif [ "$normType" == TMM ]; then
    norm=DBA_NORM_TMM
  elif [ "$normType" == RLE ]; then
    norm=DBA_NORM_RLE
  elif [ "$normType" == Loess ]; then
    norm=DBA_NORM_OFFSETS
    off=TRUE
else
  echo ERROR: PROPER NORMALIZATION WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in three positional arguments:
  echo 	"\$1 - target folder"
  echo 	"\$2 - Normalization Strategy (Default, TMM, RLE, Loess)"
  echo 	"\$3 - (Optional) The name of the generated DiffBind object (default is dbObj)"
  exit 1
fi

# Set up necessary files
mkdir -p PBS
mkdir -p log
subfolder="diffBind/${dbObj%_*}/${samp}_over_${control}_${dbObj#*_}"
mkdir -p ${subfolder}

folder=$(pwd)

cat >${folder}/PBS/${samp}over${control}_$dbObj'.R' <<EOF
# Load libraries
library(DiffBind)
library(dplyr)
library(ggplot2)

dbObj <- dba.load(file = "${dbObj}", dir = "diffBind")

mask <- dba.mask(dbObj,DBA_$caps, c("$samp","$control"), combine = "or")

dbObj <- dba(dbObj, mask=mask)
dba.normalize(dbObj, method = DBA_ALL_METHODS, normalize=$norm, offsets=$off)

dbObj <- dba.contrast(dbObj, contrast=c("${compare}","${samp}","${control}"), design = "~${compare}", minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)

png(filename = "${subfolder}/DESeq2_volcano.png", height = 1080, width = 1080)
dba.plotVolcano(dbObj, method = DBA_DESEQ2)
dev.off()

png(filename = "${subfolder}/edgeR_volcano.png", height = 1080, width = 1080)
dba.plotVolcano(dbObj, method = DBA_EDGER)
dev.off()

res_all <- dba.report(dbObj, method=DBA_ALL_METHODS, contrast = 1, th=1)
res_all

out <- as.data.frame(res_all) 
write.table(out, file="${subfolder}/${samp}over${control}.txt", sep="\t", quote=F, row.names=F)

res_deseq <- as.data.frame(dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1))

# Create a column to denote which genes are upregulated, downregulated, or not significantly changed
res_deseq <- res_deseq %>% 
	mutate( threshold_KD = case_when(Fold > 0 & FDR <= 0.05 ~ "Upregulated",
	Fold < 0 & FDR <= 0.05 ~ "Downregulated",
	TRUE ~ "Unchanged"))

## Create the volcano plot w/o any labels
p2 <- ggplot(res_deseq, aes(Fold, -log10(FDR))) +
		geom_point(aes(color = threshold_KD), size = 2/5) +
		xlab(expression("log"[2]*"FC")) +  
		ylab(expression("-log"[10]*"FDR")) +
		scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
		guides(colour = guide_legend(override.aes = list(size=1.5))) +
		theme_bw() + theme(legend.position = "none")

pdf("${subfolder}/DESeq2_pVolcano.pdf")
	print(p2)
dev.off()

res_edger <- as.data.frame(dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1))

# Create a column to denote which genes are upregulated, downregulated, or not significantly changed
res_edger <- res_edger %>% 
	mutate( threshold_KD = case_when(Fold > 0 & FDR <= 0.05 ~ "Upregulated",
	Fold < 0 & FDR <= 0.05 ~ "Downregulated",
	TRUE ~ "Unchanged"))

## Create the volcano plot w/o any labels
p2 <- ggplot(res_edger, aes(Fold, -log10(FDR))) +
		geom_point(aes(color = threshold_KD), size = 2/5) +
		xlab(expression("log"[2]*"FC")) +  
		ylab(expression("-log"[10]*"FDR")) +
		scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
		guides(colour = guide_legend(override.aes = list(size=1.5))) +
		theme_bw() + theme(legend.position = "none")

pdf("${subfolder}/edgeR_pVolcano.pdf")
	print(p2)
dev.off()

unchanged <- as.data.frame(res_all) %>%
  filter(FDR > 0.05) %>% 
  select(seqnames, start, end)

write.table(unchanged, file="${subfolder}/unchanged_long.bed", sep="\t", quote=F, row.names=F)

allDB <- out %>% 
  filter((Fold < 0 & FDR < 0.05) | (Fold > 0 & FDR < 0.05)) %>%
  select(seqnames, start, end, Fold, p.value, FDR)

allDBbed <- allDB %>% select(seqnames, start, end)

write.table(allDBbed, file="${subfolder}/allDB.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(allDB, file="${subfolder}/allDB.txt", sep="\t", quote=F, row.names=F, col.names=F)

KO_enrich <- out %>% 
  filter(Fold < 0 & FDR < 0.05) %>% 
  select(seqnames, start, end)
  
# Write to bed file
write.table(KO_enrich, file="${subfolder}/downregulated.bed", sep="\t", quote=F, row.names=F, col.names=F)

WT_enrich <- out %>% 
  filter(Fold > 0 & FDR < 0.05) %>% 
  select(seqnames, start, end)
  
# Write to bed file
write.table(WT_enrich, file="${subfolder}/upregulated.bed", sep="\t", quote=F, row.names=F, col.names=F)

# Note: These have been moved to the bottom because if one group is zero, then the program automatically terminates
png(filename = "${subfolder}/venn_DBvsNonDB.png", height = 1080, width = 1080)
dba.plotVenn(dbObj, contrast=1, main="${samp} over ${control}", bDB=TRUE, bNotDB=TRUE)
dev.off()

png(filename = "${subfolder}/venn.png", height = 1080, width = 1080)
dba.plotVenn(dbObj, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)
dev.off()
EOF

cat >${folder}/PBS/${samp}over${control}_$dbObj.sbatch <<EOF
#!/bin/bash -l
#SBATCH --job-name=${samp}over${control}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${samp}over${control}_%j.txt -e ${folder}/log/${samp}over${control}_%j.err.txt

#------- END OF HEADER -------#
source activate ATACQC

cd ${folder}

Rscript ./PBS/${samp}over${control}_$dbObj'.R'

tail -n +2 "${subfolder}/unchanged_long.bed" > ${subfolder}/unchanged.bed
rm ${subfolder}/unchanged_long.bed

sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/backend/downstreamDiffBind/diffBind_HomerMotif.sh ${subfolder}
sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/backend/downstreamDiffBind/diffBind_ChIPseeker.sh ${subfolder}
EOF

sbatch ${folder}/PBS/${samp}over${control}_$dbObj.sbatch 