#!/bin/bash

# Creates a heatmap of the upregulated, unchanged, and downregulated regions an the appropriate signal
# Creates both single-region heatmaps, as well as a compiled heatmap with three rows

# Format for the meta text:
# diffBindSubfolder	left	right	color	subfolderLocation(optional)

# Setup - Create required folders
mkdir -p PBS
mkdir -p log
mkdir -p triple_heatmaps

if [ -z "$1" ]; then 
  echo ERROR: META FILE WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target meta file containing sample information
  echo Meta file shouled be tab-delimited, and formatted as follows:
  echo EXAMPLE:
  echo diffBindSubfolder	left	right	color	subfolderLocation \(optional\)
  exit 1
fi

folder=$(cd "$(dirname "$0")";pwd)
libraryText=$1

# Parse through libraryText, extracting the relevant information
OLDIFS=$IFS
while IFS=$'\t' read -r -a varArray; do
	dbFolder=${varArray[0]}  # Top row
	lColumn=${varArray[1]}  # Left column
	rColumn=${varArray[2]}  # Right column
	color=${varArray[3]}  # Right column
	subfolder=${varArray[4]}  # Location of subfolder to dump heatmaps into
	IFS=$OLDIFS

	mkdir -p triple_heatmaps/${subfolder}

	# Strip suffixes to make it nicer
	lc=$(basename "$lColumn" "_merged.bw")
	lc=$(basename "$lc" "_normalized.bw")

	rc=$(basename "$rColumn" "_merged.bw")
	rc=$(basename "$rc" "_normalized.bw")
	name=$(basename ${dbFolder%'/'} "")_tripleHeatmaps

cat >${folder}/PBS/${name}.sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${name}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${lc}_${rc}_%j.txt -e ${folder}/log/${lc}_${rc}_%j.err.txt
cd ${folder}
#------- End of header -------#
source activate alignment

# Single region; Upregulated
computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R ${dbFolder}upregulated.bed \
	-S ${lColumn} ${rColumn} \
	-o triple_heatmaps/${subfolder}/${lc}_${rc}_up.gz --missingDataAsZero -p max --smartLabels

plotHeatmap -m triple_heatmaps/${subfolder}/${lc}_${rc}_up.gz \
	-out triple_heatmaps/${subfolder}/${lc}_${rc}_${color}_up.pdf --colorList 'white,${color}'

# Single region; Unchanged
computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R ${dbFolder}unchanged.bed \
	-S ${lColumn} ${rColumn} \
	-o triple_heatmaps/${subfolder}/${lc}_${rc}_stays.gz --missingDataAsZero -p max --smartLabels

plotHeatmap -m triple_heatmaps/${subfolder}/${lc}_${rc}_stays.gz \
	-out triple_heatmaps/${subfolder}/${lc}_${rc}_${color}_stays.pdf --colorList 'white,${color}'

# Single region; Downregulated
computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R ${dbFolder}downregulated.bed \
	-S ${lColumn} ${rColumn} \
	-o triple_heatmaps/${subfolder}/${lc}_${rc}_down.gz --missingDataAsZero -p max --smartLabels

plotHeatmap -m triple_heatmaps/${subfolder}/${lc}_${rc}_down.gz \
	-out triple_heatmaps/${subfolder}/${lc}_${rc}_${color}_down.pdf --colorList 'white,${color}'

# All regions
computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R ${dbFolder}upregulated.bed ${dbFolder}unchanged.bed ${dbFolder}downregulated.bed \
	-S ${lColumn} ${rColumn} \
	-o triple_heatmaps/${subfolder}/${lc}_${rc}_triple.gz --missingDataAsZero -p max --smartLabels

plotHeatmap -m triple_heatmaps/${subfolder}/${lc}_${rc}_triple.gz \
	-out triple_heatmaps/${subfolder}/${lc}_${rc}_${color}_triple.pdf --colorList 'white,${color}'

rm -r ${subfolder}/${lc}_${rc}*.gz
EOF

	sbatch ${folder}/PBS/${lc}_${rc}_proximalDistal.pbs
done < ${libraryText}