#!/bin/bash

# Peakcalling pipeline for CUT&RUN data. Modified from existing pipeline 8/8/2022
# Designed to work in conjunction with cr_alignment.sh, and automatically call HOMER and ChIPseeker motifs
# Requires the peakcalling environment, and the ChIPseeker environment (for later)
# See areas marked "CHANGE FOR FILE EXTENSION"

# Requires the "alignment" environment, and the "qc" environment (for fastqc and multiqc)
# This pipeline takes in three parameters:
#   $1 - FDR cutoff
#   $2 - log cutoff
#   $3 - counts cutoff

# Version Doc: https://docs.google.com/document/d/1HMoJZ5QpscD6h2Gb1v2hWVpKfk1sxz4MzuVZoArSHyw/edit

# TODO: Add EChO

# Default parameter cutoffs
fdr=0.01
lfc=1
cnts=20

suffix1=.sorted.filtered.bam

# If nothing is provided, then use the defaults
if [ -z "$1" ]; then 
  echo Using the default parameters:
  echo "For future reference, this program takes in three parameters:"
  echo "\$1 - FDR cutoff (Default <= 0.01)"
  echo "\$2 - log fold change cutoff (Default >= 1)"
  echo "\$3 - counts cutoff (Default 20)"
elif [ -z "$3" ]; then
  echo "Error. Either pass in all three parameters, or pass in none to use the defaults"
else
  fdr=$1
  lfc=$2
  cnts=$3
fi

# Setup: Create needed folders
mkdir -p PBS
mkdir -p log
mkdir -p mergedBed
mkdir -p mergedHeatmaps

folder=$(pwd)

count=$(find ./mergedBam -mindepth 1 -type f -name "*_sortedMerged.bam" -printf x | wc -c)  # Finds total number of files matching extension. CHANGE FOR FILE EXTENSION

echo "${count} peakcalling to be done..."

# Meta file to know when qc will begin (empty file)
cat >${folder}/'pcMeta.txt' <<EOF
EOF

cd ${folder}/aligned
for file in *IgG*${suffix1}; do
    igGbase=$(basename "$file" "${suffix1}")
    groupprefix=${igGbase%%_IgG*}
    echo $igGbase
    for f in $folder/mergedBam/${groupprefix}*_sortedMerged.bam; do
        base=$(basename "$f" "_sortedMerged.bam")
        if [ "$base" != "$igGbase" ]; then
            cat >${folder}/PBS/${base}_mergePeakcalling'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${base}_mergePC  # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=32GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=2:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/${base}_mergePeakcalling.%j.out
#SBATCH --error=${folder}/log/${base}_mergePeakcalling.%j.err
################################
# Enter your code to run below #
################################

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate peakcalling

cd ${folder}

epic2 \
  --treatment $folder/mergedBam/${base}_sortedMerged.bam \
  --control $folder/aligned/$file \
  --genome hg38 \
  -o mergedBed/${base}.bed

awk '\$9 <= ${fdr} {print \$0}' mergedBed/${base}.bed > mergedBed/${base}_FDR_${fdr}.bed

awk '\$10 >= ${lfc} {print \$0}' mergedBed/${base}_FDR_${fdr}.bed > mergedBed/${base}_log_${lfc}.bed

awk '\$7 > ${cnts} {print \$0}' mergedBed/${base}_log_${lfc}.bed > mergedBed/${base}.bed

rm mergedBed/${base}_FDR_${fdr}.bed
rm mergedBed/${base}_log_${lfc}.bed

conda activate alignment

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R mergedBed/${base}.bed \
	-S mergedBigWig/${base}_merged.bw  \
	-o $folder/mergedHeatmaps/${base}_center.gz --missingDataAsZero -p max

plotHeatmap -m $folder/mergedHeatmaps/${base}_center.gz \
	-out $folder/mergedHeatmaps/${base}_center.pdf --colorList 'white,darkred'

rm $folder/mergedHeatmaps/${base}_center.gz

echo "${base} completed!" >> ${folder}/'pcMeta.txt'

currLine=\$(wc -l < ${folder}/pcMeta.txt)
if ((\$currLine == $count)); then
    source activate base
    sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/homerMotif.sh mergedBed
    sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/ChIPseeker.sh mergedBed
    rm ${folder}/pcMeta.txt
fi
EOF
        echo ${base}
        sbatch ${folder}/PBS/${base}_mergePeakcalling'.pbs'
        fi
    done
done