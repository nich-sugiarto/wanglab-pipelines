#!/bin/bash

# The first line should be two tab-delimited fields long: The region file and the sorted 
# file whose sorting will determine the order for the rest of the heatmaps. This will only
# provide the .bed file. The second line and onwards should have the three following fields:
# bigWigFiles (multiple should be space delimited)	color	name

# Script automatically generated using the "slouch" command on Mon Apr 17 15:19:11 EDT 2023

folder=$(cd "$(dirname "$0")";pwd)  # Stores current folder as a variable

# Set up required files
mkdir -p PBS
mkdir -p log
mkdir -p distProx
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

libraryText=$1
if [ -z "$2" ]; then
	firstLine=$(head $1 -n 1)
	IFS=$'\t' read -ra sortFiles <<< "$firstLine"
	bigWig=${sortFiles[0]}
	bedFile=${sortFiles[1]}
	subfolder=${sortFiles[2]}
	cName=$(basename $bedFile ".bed")
	echo $cName

	mkdir -p distProx/$subfolder

	# Generate distal regions by removing the proximal sites
	bedtools intersect \
		-v -a ${bedFile} \
		-b /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/TSSsites_ARID1A/hg38_TSS_UCSC_paddedUnique.bed \
		> distProx/$subfolder/TSSdistal.bed

	# Generate distal regions by removing the proximal sites
	bedtools intersect \
		-wa -a ${bedFile} \
		-b /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/TSSsites_ARID1A/hg38_TSS_UCSC_paddedUnique.bed \
		> distProx/$subfolder/TSSproximal.bed

	cat >${folder}/PBS/${cName}_makeSortedMatrix.pbs <<EOF
#!/bin/bash
#SBATCH --job-name=${cName}_makeSortedMatrix
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${cName}_makeSortedMatrix_%j.txt -e ${folder}/log/${cName}_makeSortedMatrix_%j.err.txt
cd ${folder}
#------- End of header -------#
source activate alignment

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R distProx/$subfolder/TSSproximal.bed \
	--sortRegions descend \
	-S $bigWig \
	-o distProx/$subfolder/TEMP_proximal.gz --outFileSortedRegions distProx/$subfolder/TSSproximal_sorted.bed --missingDataAsZero -p max --smartLabels

plotHeatmap -m distProx/TEMP_proximal.gz \
	--sortRegions no \
	-out distProx/$subfolder/TEMP_proximal.pdf --colorList 'white,red'

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R distProx/$subfolder/TSSdistal.bed \
	--sortRegions descend \
	-S $bigWig \
	-o distProx/$subfolder/TEMP_distal.gz --outFileSortedRegions distProx/$subfolder/TSSdistal_sorted.bed --missingDataAsZero -p max --smartLabels

plotHeatmap -m distProx/$subfolder/TEMP_distal.gz \
	--sortRegions no \
	-out distProx/$subfolder/TEMP_distal.pdf --colorList 'white,red'

sed 's/genes/proximal/g' distProx/$subfolder/TSSproximal_sorted.bed > distProx/$subfolder/TSSproximal_sortedNamed.bed
sed 's/genes/distal/g' distProx/$subfolder/TSSdistal_sorted.bed > distProx/$subfolder/TSSdistal_sortedNamed.bed

tail -n +2 "$1" > distProx/$subfolder/TEMP.txt
sh sortProximalDistal.sh distProx/$subfolder/TEMP.txt plot $subfolder
EOF

sbatch ${folder}/PBS/${cName}_makeSortedMatrix.pbs
elif [ "$2" = "plot" ]; then
	allFiles=""
	allColors=""
	OLDIFS=$IFS
	subfolder=$3
	while IFS=$'\t' read -r -a varArray; do
		bigWigFiles=${varArray[0]}  # The control group
		name=${varArray[1]}  # Name of heatmap
		color=${varArray[2]}  # Color hexcode
		IFS=$'\t' read -ra numBigWigs <<< "$bigWigFiles"
		for i in $numBigWigs; do
			allColors+="'white,${color}' "
		done
		allFiles+="$bigWigFiles "
		IFS=$OLDIFS
		cat >${folder}/PBS/${name}_plotSortedHeatmaps'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${name}_plotSortedHeatmaps # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=32GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=12:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/${name}_plotSortedHeatmaps%j.out
#SBATCH --error=${folder}/log/${name}_plotSortedHeatmaps%j.err
################################
# Enter your code to run below #
################################
cd ${folder}

source activate alignment

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R distProx/$subfolder/TSSproximal_sortedNamed.bed distProx/$subfolder/TSSdistal_sortedNamed.bed \
	-S ${bigWigFiles} \
	-o distProx/$subfolder/${name}_TSS.gz --missingDataAsZero -p max --smartLabels

plotHeatmap -m distProx/$subfolder/${name}_TSS.gz \
	--sortRegions no \
	-out distProx/$subfolder/${name}_TSS.pdf --colorList 'white,${color}'
EOF
	sbatch ${folder}/PBS/${name}_plotSortedHeatmaps'.pbs'
	done < ${libraryText}
	allColors=${allColors/%" "/""}
	allFiles=${allFiles/%" "/""}
	cat >${folder}/PBS/all_plotSortedHeatmaps'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=all_plotSortedHeatmaps # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=32GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=12:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/all_plotSortedHeatmaps%j.out
#SBATCH --error=${folder}/log/all_plotSortedHeatmaps%j.err
################################
# Enter your code to run below #
################################
cd ${folder}

source activate alignment

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R distProx/$subfolder/TSSproximal_sortedNamed.bed distProx/$subfolder/TSSdistal_sortedNamed.bed \
	-S ${allFiles} \
	-o distProx/$subfolder/all_TSS.gz --missingDataAsZero -p max --smartLabels

plotHeatmap -m distProx/$subfolder/all_TSS.gz \
	--sortUsingSamples 1 \
	-out distProx/$subfolder/all_TSS.pdf --colorList ${allColors}
EOF
	sbatch ${folder}/PBS/all_plotSortedHeatmaps.pbs
else
	echo ERROR!
	echo Invalid option specified
	echo "\$2 should either have the value 'plot', or no value at all"
	exit 1
fi