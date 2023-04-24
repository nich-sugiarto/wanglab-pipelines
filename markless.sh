#!/bin/bash

# Creates a heatmap from regions that lacks co-binding with another mark

# Sample body for the meta text
# regionToExamine	regionsToRemove	leftHeatmapColumn	midHeatmapColumn	rightHeatmapColumn	color

mkdir -p PBS
mkdir -p log
mkdir -p markless
mkdir -p no_A1A_heatmaps

if [ -z "$1" ]; then 
	echo ERROR: META FILE WAS NOT SPECIFIED
	echo USAGE:
	echo This pipeline takes in one positional argument:
	echo 	\$1 - target meta file containing sample information
	echo The meta file should be tab delimited with two columns, and formatted as follows:
	echo regionToExamine	regionsToRemove	leftColumn	midColumn	rightColumn	color
	exit 1
fi

folder=$(pwd)
libraryText=$1

# Read in library text, and parse accordingly
OLDIFS=$IFS
while IFS=$'\t' read -r -a varArray; do
	markFile=${varArray[0]}  # Top row
	noMark=${varArray[1]}  # Bottom row
	lColumn=${varArray[2]}  # Left column
	mColumn=${varArray[3]}  # Middle column
	rColumn=${varArray[4]}  # Right column
	color=${varArray[5]}  # Right column
	IFS=$OLDIFS

mark=$(basename "$markFile" ".bed")

# Remove regions
bedtools intersect \
	-v -a ${markFile} \
	-b ${noMark} \
	> markless/${mark}_mark.bed

# Find the distal regions of that removed 
bedtools intersect \
	-v -a markless/${mark}_markless.bed \
	-b /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/TSSsites_ARID1A/hg38_TSS_UCSC_paddedUnique.bed \
	> markless/${mark}_TSSdistal.bed

bedtools intersect \
	-wa -a markless/${mark}_markless.bed \
	-b /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/TSSsites_ARID1A/hg38_TSS_UCSC_paddedUnique.bed \
	> markless/${mark}_TSSproximal.bed

lc=$(basename "$lColumn" "_merged.bw")
lc=$(basename "$lc" "_normalized.bw")

mc=$(basename "$mColumn" "_merged.bw")
mc=$(basename "$mc" "_normalized.bw")
echo ${lc} ${mc}

cat >${folder}/PBS/${lc}_${mc}_markless.pbs <<EOF
#!/bin/bash
#SBATCH --job-name=${lc}_${mc}_proximalDistal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${lc}_${mc}_%j.txt -e ${folder}/log/${lc}_${mc}_%j.err.txt
cd ${folder}
#------- End of header -------#
source activate alignment

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R markless/${mark}_TSSproximal.bed markless/${mark}_TSSdistal.bed \
	-S ${lColumn} ${mColumn} ${rColumn} \
	-o no_A1A_heatmaps/${lc}_${mc}_TSS.gz --missingDataAsZero -p max --smartLabels

plotHeatmap -m no_A1A_heatmaps/${lc}_${mc}_TSS.gz \
	-out no_A1A_heatmaps/${lc}_${mc}_${color}_TSS.pdf --colorList 'white,${color}' 'white,${color}' 'white, #036603'
EOF

sbatch ${folder}/PBS/${lc}_${rc}_markless.pbs
done < ${libraryText}