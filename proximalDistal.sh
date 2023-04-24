#!/bin/bash

# Find the proximal and distal regions of merged bed filess

# Sample body for the meta text
# WT_A1A	#1e99f2
# WT_K27Ac	#a73ea7
# WT_K4me3	#5252b0

# Setup - generate required folders
mkdir -p PBS
mkdir -p log
mkdir -p distProx
mkdir -p TSSheatmap

# Error message if meta file is not provided
if [ -z "$1" ]; then 
  echo ERROR: META FILE WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target meta file containing sample information
  echo The meta file should be tab delimited with two columns, and formatted as follows:
  echo mark	colorOfHeatmap
  exit 1
fi

folder=$(cd "$(dirname "$0")";pwd)
libraryText=$1

OLDIFS=$IFS
while IFS=$'\t' read -r -a varArray; do
	markLoc=${varArray[0]}  # The control group
	color=${varArray[1]}  # The treatment group
	IFS=$OLDIFS

mark=$(basename "$markLoc" ".bed")
echo ${mark}

# Generate distal regions by removing the proximal sites
bedtools intersect \
	-v -a ${markLoc} \
	-b /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/TSSsites_ARID1A/hg38_TSS_UCSC_paddedUnique.bed \
	> distProx/${mark}_TSSdistal.bed

# Generate distal regions by removing the proximal sites
bedtools intersect \
	-wa -a ${markLoc} \
	-b /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/TSSsites_ARID1A/hg38_TSS_UCSC_paddedUnique.bed \
	> distProx/${mark}_TSSproximal.bed

cat >${folder}/PBS/${mark}_proximalDistal.pbs <<EOF
#!/bin/bash
#SBATCH --job-name=${mark}_proximalDistal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${mark}_%j.txt -e ${folder}/log/${mark}_%j.err.txt
cd ${folder}
#------- End of header -------#
source activate alignment

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R distProx/${mark}_TSSproximal.bed distProx/${mark}_TSSdistal.bed \
	-S $folder/normalized_bw/${mark}_normalized.bw \
	-o $folder/TSSheatmap/${mark}_TSS.gz --missingDataAsZero -p max --smartLabels

plotHeatmap -m $folder/TSSheatmap/${mark}_TSS.gz \
	-out $folder/TSSheatmap/${mark}_TSS.pdf --colorList 'white,${color}'
EOF

sbatch ${folder}/PBS/${mark}_proximalDistal.pbs
done < ${libraryText}