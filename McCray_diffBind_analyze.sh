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
  dbObj="dbObj.rds"
fi

# Set up necessary files
mkdir -p PBS
mkdir -p log
mkdir -p diffBind/${samp}_over_${control}

folder=$(cd "$(dirname "$0")";pwd)
subfolder="diffBind/${samp}_over_${control}"

cat >${folder}/PBS/${samp}over${control}'.R' <<EOF
# Load libraries
suppressMessages(library(DiffBind))
suppressMessages(library(dplyr))
suppressMessages(library(csaw))
suppressMessages(library(profileplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(GetoptLong))

dbObj <- dba.load(file = "${dbObj}", dir = "diffBind/")
scores <- dba.peakset(dbObj, bRetrieve=T, DataType=DBA_DATA_FRAME)
write.table(scores, file="${folder}/diffBind/peakScores.bed", sep="\t", quote=F, row.names=F)

png(filename = "${subfolder}/Correlations.png", height = 1080, width = 1080)
plot(dbObj)
dev.off()

png(filename = "${subfolder}/Heatmap.png", height = 1080, width = 1080)
dba.plotHeatmap(dbObj, correlations=FALSE, scale="row")
dev.off()

dbObj <- dba.contrast(dbObj, contrast=c("${compare}","${samp}","${control}"), design = "~${compare} + Replicate", minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_EDGER)

png(filename = "${subfolder}/DiffCorrelations.png", height = 1080, width = 1080)
plot(dbObj, contrast=1, method=DBA_EDGER)
dev.off()

png(filename = "${subfolder}/volcano.png", height = 1080, width = 1080)
dba.plotVolcano(dbObj, method=DBA_EDGER)
dev.off()

  res_deseq <- dba.report(dbObj, method=DBA_EDGER, bNormalized = TRUE, contrast = 1, th=1)
  res_deseq
  out <- as.data.frame(res_deseq) 
  write.table(out, file=paste0("${subfolder}","/${samp}over${control}.txt"), sep="\t", quote=F, row.names=F)

#   unchanged <- as.data.frame(res_deseq) %>%
#     filter(FDR > 0.05) %>% 
#     select(seqnames, start, end)
#   write.table(unchanged, file=paste0("${subfolder}","/unchanged_long.bed"), sep="\t", quote=F, row.names=F)

  allDB <- out %>% 
    filter((Fold < 0 & FDR < 0.05) | (Fold > 0 & FDR < 0.05)) %>%
    select(seqnames, start, end, Fold, p.value, FDR)
  write.table(allDB, file=paste0("${subfolder}","/allDB.bed"), sep="\t", quote=F, row.names=F, col.names=F)

  KO_enrich <- out %>% 
    filter(Fold < 0 & FDR < 0.05) %>% 
    select(seqnames, start, end, Fold, p.value, FDR)
  write.table(KO_enrich, file=paste0("${subfolder}","/downregulated.bed"), sep="\t", quote=F, row.names=F, col.names=F)

  WT_enrich <- out %>% 
    filter(Fold > 0 & FDR < 0.05) %>% 
    select(seqnames, start, end, Fold, p.value, FDR)
  write.table(WT_enrich, file=paste0("${subfolder}","/upregulated.bed"), sep="\t", quote=F, row.names=F, col.names=F)
  
  #These plots cancel the job if diff sites is zero so put down here
  if (nrow(KO_enrich) == 0 || nrow(WT_enrich) == 0){
  } else {
      png(filename = paste0("${subfolder}","/venn.png"), height = 1080, width = 1080)
      dba.plotVenn(dbObj, method=DBA_EDGER, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)
      dev.off()

      png(filename = paste0("${subfolder}","/venn_DBvsNonDB.png"), height = 1080, width = 1080)
      dba.plotVenn(dbObj, method=DBA_EDGER, contrast=1, main="${samp} over ${control}", bDB=TRUE, bNotDB=TRUE)
      dev.off()

      png(filename = "${subfolder}/DiffHeatmap.png", height = 1080, width = 1080)
      dba.plotHeatmap(dbObj, contrast=1, method=DBA_EDGER, correlations=FALSE, scale="row")
      dev.off()

      png(filename = "${subfolder}/profile.png", height = 1080, width = 1080)
      dba.plotProfile(dbObj)
      dev.off()
  }


EOF

cat >${folder}/PBS/${samp}over${control}.sbatch <<EOF
#!/bin/bash -l
#SBATCH --job-name=${samp}over${control}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${samp}over${control}_%j.txt -e ${folder}/log/${samp}over${control}_%j.err.txt

#------- END OF HEADER -------#
source activate ATACQC

cd ${folder}

Rscript ./PBS/${samp}over${control}'.R'

# tail -n +2 "${subfolder}/unchanged_long.bed" > ${subfolder}/unchanged.bed
# rm ${subfolder}/unchanged_long.bed

cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/downstreamDiffBind/diffBind_HomerMotif.sh ${folder}
sh diffBind_HomerMotif.sh ${subfolder}
cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/downstreamDiffBind/diffBind_ChIPseeker.sh ${folder}
sh diffBind_ChIPseeker.sh ${subfolder}
EOF

sbatch ${folder}/PBS/${samp}over${control}.sbatch