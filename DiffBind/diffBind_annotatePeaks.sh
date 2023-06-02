# Annotates only the differentially bound peaks identified via diffBind

# This pipeline requires one positional argument:
# 	$1 - The meta file that has all the comparisons that need to be made
# The file should be tab-delimited, and contain one field:
# nameOfFolderUnderDiffBind

# Requires the ATACQC pipeline (will rename in the future)

# TODO: Merge with diffBind_analyze.sh such that if only one field is provided, assume that it is a meta file
# But, if three categories are provided, use those three categories to run as it is right now. 
# Break on two

#  All the parameters are not provided
if [ -z "$1" ]; then 
  echo ERROR: NOT ALL PARAMETERS HAVE BEEN SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo \$1 - The meta file that has all the comparisons that need to be made
  echo The file should be tab-delimited, and contain the following fields:
  echo Name of DiffBind folder
fi

folder=$(pwd)

libraryText=$1  # Passed in file name

mkdir -p PBS
mkdir -p log
mkdir -p linkedGenes

OLDIFS=$IFS
while IFS=$'\t' read -r -a varArray; do
	location=${varArray[0]}  # Category
	IFS=$OLDIFS
  mkdir -p geneLists/diffBind/${location}

  cat >${folder}/PBS/${location}_linker'.R' <<EOF
library(tidyverse)
library(dplyr)
library(janitor)

getwd()

# Remove ' as a quote because of 5' and 3' annotations
peaks <- read.table("geneLists/diffBind/${location}/${location}_noBlanks.txt", sep = "\t", header = TRUE, quote = "", fill = TRUE) %>% clean_names()
mRNA <- peaks[str_sub(peaks\$"nearest_promoter_id", 2, 2)  == "M", ]
colnames(mRNA)[1] <- "peakID"

dbVal <- read.table("${folder}/diffBind/${location}/allDB.txt", sep = "\t", header = FALSE) %>% clean_names()
colnames(dbVal) <- c("seqnames", "start", "end", "Fold", "p.value", "FDR")

str(dbVal)
str(mRNA)

mRNA_fin <- mRNA[,-1]
rownames(mRNA_fin) <- mRNA[,1]

finDf <- merge(mRNA_fin, dbVal, by='row.names')
str(finDf)
write.csv(finDf, "./linkedGenes/${location}_linkedPeakGenes.csv")

diffList <- finDf %>%
  select(gene_name, Fold, p.value, FDR)

diffList <- diffList[!diffList\$"gene_name" == "",]
write.csv(diffList, "./linkedGenes/${location}_foldTables.csv")
EOF

	cat >${folder}/PBS/${location}_diffBindPileup'.pbs' <<EOF
#!/bin/bash
#SBATCH --job-name=${location}_annoPeaks
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH -o ${folder}/log/${location}_diffBindPileup_%j.txt -e ${folder}/log/${location}_diffBindPileup_%j.err.txt
cd ${folder}
################################
# Enter your code to run below #
################################
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate ChIPseeker

annotatePeaks.pl ${folder}/diffBind/${location}/allDB.bed hg38 > ${folder}/geneLists/diffBind/${location}/dbRegions_annotated.txt

cat ${folder}/geneLists/diffBind/${location}/dbRegions_annotated.txt | \
	awk -F '\t' -v OFS='\t' ' \$11 != "NA" \
	{ print \$1 "	" \$2 "	" \$3 "	" \$4 "	" \$5 "	" \$6 "	" \$7 "	" \$8 "	" \$9 "	" \$10 "	" \$11 "	" \$16}' \
	> ${folder}/geneLists/diffBind/${location}/${location}_noBlanks.txt

Rscript ${folder}/PBS/${location}_linker'.R'
EOF

  sbatch ${folder}/PBS/${location}_diffBindPileup'.pbs'
done < ${libraryText}