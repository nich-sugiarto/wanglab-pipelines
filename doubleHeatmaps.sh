#!/bin/bash

# Creates heatmaps for both distal and proximal regions (that are passed in)

# Sample body for the meta text
# prox	dist	left	right	color

mkdir -p PBS
mkdir -p log
mkdir -p A1A_heatmaps

# Check to makes sure that the appropriate arguments were passed in
if [ -z "$1" ]; then 
  echo ERROR: META FILE WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target meta file containing sample information
  exit 1
fi

folder=$(pwd)
libraryText=$1

# Parses through the provided file
OLDIFS=$IFS
while IFS=$'\t' read -r -a varArray; do
	proxRegion=${varArray[0]}  # Top row
	distRegion=${varArray[1]}  # Bottom row
	lColumn=${varArray[2]}  # Left column
	rColumn=${varArray[3]}  # Right column
	color=${varArray[4]}  # Right column
	IFS=$OLDIFS

lc=$(basename "$lColumn" "_merged.bw")
lc=$(basename "$lc" "_normalized.bw")

rc=$(basename "$rColumn" "_merged.bw")
rc=$(basename "$rc" "_normalized.bw")
echo ${lc} ${rc}

cat >${folder}/PBS/${lc}_${rc}_proximalDistal.pbs <<EOF
#!/bin/bash
#SBATCH --job-name=${lc}_${rc}_proximalDistal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${lc}_${rc}_%j.txt -e ${folder}/log/${lc}_${rc}_%j.err.txt
cd ${folder}
#------- End of header -------#
source activate alignment

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R ${proxRegion} ${distRegion} \
	-S ${lColumn} ${rColumn} \
	-o A1A_heatmaps/${lc}_${rc}_TSS.gz --missingDataAsZero -p max --smartLabels

plotHeatmap -m A1A_heatmaps/${lc}_${rc}_TSS.gz \
	-out A1A_heatmaps/${lc}_${rc}_${color}_TSS.pdf --colorList 'white,${color}'
EOF

sbatch ${folder}/PBS/${lc}_${rc}_proximalDistal.pbs
done < ${libraryText}