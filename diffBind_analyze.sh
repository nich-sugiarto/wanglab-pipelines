#!/bin/bash

# The log fold change will be calculated as log(samp/control), or log(samp) - log(control)

# Script that creates a diffBind comparison, and completes the standard analysis

# This pipeline requires four positional arguments:
# 	$1 - The DiffBind category on which to make the comparison
# 	$2 - First category for comparison
# 	$3 - Second category for comparison
# 	$4 - (Optional) name of the object from diffBind_createObject.sh to read in (default is 'dbObj.rds')

# NOTE: This assumes that there are only two reps. If there are more than two reps, then you must manually change that parameter within the script
# Version control doc avaible at https://docs.google.com/document/d/12GkUyHpHb-jCVwmvbipUhdPgJDdXWWneKmI7nTRCfYc/edit

# Requires the ATACQC pipeline (will rename in the future)

compare=$1
samp=$2
control=$3
dbObj=$4

#  All the parameters are not provided
if [ -z "$3" ]; then 
  echo ERROR: NOT ALL PARAMETERS HAVE BEEN SPECIFIED
  echo USAGE:
  echo This pipeline takes in a minimum of three positional argument:
  echo '\$1 - The DiffBind category on which to make the comparison'
  echo '\$2 - First category for comparison'
  echo '\$3 - Second category for comparison'
  echo "\$4 - (Optional) name of the object from diffBind_createObject.sh to read in (default is 'dbObj.rds')"
  exit 1
fi

if [ -z "$4" ]; then 
  dbObj="dbObj.rds"
fi

folder=$(cd "$(dirname "$0")";pwd)

mkdir -p diffBind
mkdir -p diffBind/${samp}_vs_${control}

subfolder="diffBind/${samp}_vs_${control}"

cat >${folder}/PBS/${samp}vs${control}'.R' <<EOF
# Load libraries
library(DiffBind)
library(dplyr)

dbObj <- readRDS("${dbObj}")

dbObj <- dba.contrast(dbObj, contrast=c("${compare}","${samp}","${control}"), design = "~${compare} + Replicate", minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)

png(filename = "${subfolder}/venn.png", height = 1080, width = 1080)
dba.plotVenn(dbObj, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)
dev.off()

png(filename = "${subfolder}/venn_nonDB.png", height = 1080, width = 1080)
dba.plotVenn(dbObj, contrast=1, main="${samp} vs ${control}", bDB=TRUE, bNotDB=TRUE)
dev.off()

png(filename = "${subfolder}/DESeq2_volcano.png", height = 1080, width = 1080)
dba.plotVolcano(dbObj, method = DBA_DESEQ2)
dev.off()

png(filename = "${subfolder}/edgeR_volcano.png", height = 1080, width = 1080)
dba.plotVolcano(dbObj, method = DBA_EDGER)
dev.off()

res_deseq <- dba.report(dbObj, method=DBA_ALL_METHODS, contrast = 1, th=1)
res_deseq
out <- as.data.frame(res_deseq)
write.table(out, file="results/${samp}vs${control}.txt", sep="\t", quote=F, row.names=F)

KO_enrich <- out %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  select(seqnames, start, end, Fold, p.value, FDR)
  
# Write to bed file
write.table(KO_enrich, file="${subfolder}/downregulated.bed", sep="\t", quote=F, row.names=F, col.names=T)

WT_enrich <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  select(seqnames, start, end, Fold, p.value, FDR)
  
# Write to bed file
write.table(WT_enrich, file="${subfolder}/upregulated.bed", sep="\t", quote=F, row.names=F, col.names=T)
EOF

cat >${folder}/PBS/${samp}vs${control}.sbatch <<EOF
#!/bin/bash -l
#SBATCH --job-name=${samp}vs${control}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${samp}vs${control}_%j.txt -e ${folder}/log/${samp}vs${control}_%j.err.txt

#------- END OF HEADER -------#
source activate ATACQC

cd ${folder}

Rscript ./PBS/${samp}vs${control}'.R'
EOF

sbatch ${folder}/PBS/${samp}vs${control}.sbatch