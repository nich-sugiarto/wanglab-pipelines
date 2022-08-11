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

# TODO: Test pipeline

fdr=0.01
lfc=1
cnts=20
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
mkdir -p epic2
mkdir -p heatmaps

folder=$(cd "$(dirname "$0")";pwd)

count=$(find ./fastq -mindepth 1 -type f -name "*_1.fq.gz" -printf x | wc -c)  # Finds total number of files matching extension. CHANGE FOR FILE EXTENSION
groupCount=$(find ./fastq -mindepth 1 -type f -name "*IgG*_1.fq.gz" -printf x | wc -c)  # Finds total number of control (IgG) files. CHANGE FOR FILE EXTENSION

((count-=groupCount))

echo "${count} peakcalling to be done..."

# Meta file to know when qc will begin (empty file)
cat >${folder}/'meta.txt' <<EOF
EOF

cd fastq
for file in *_IgG_*_1.fq.gz; do
    igGbase=$(basename "$file" "_1.fq.gz")
    igGbase=${igGbase%_CKDL*}
    groupprefix=${igGbase%%_*}

    for f in $folder/fastq/*_1.fq.gz; do
        base=$(basename "$f" "_1.fq.gz")
        prefix=${base%%_*}
        smallBase=${base%_CKDL*}
        if [ "$smallBase" != "$igGbase" ] && [ "$prefix" == "$groupprefix" ]; then
            echo ${smallBase}
            cat >${folder}/PBS/${smallBase}_epic2'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=epic2_${smallBase}  # Name of the job

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
#SBATCH --output=${folder}/log/${smallBase}_epic2.%j.out
#SBATCH --error=${folder}/log/${smallBase}_epic2.%j.err
################################
# Enter your code to run below #
################################

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate peakcalling

cd ${folder}

epic2 \
  --treatment $folder/aligned/${smallBase}.sorted.filtered.bam \
  --control $folder/aligned/${igGbase}.sorted.filtered.bam \
  --genome hg38 \
  -o epic2/${smallBase}.bed

awk '\$9 <= ${fdr} {print \$0}' epic2/${smallBase}.bed > epic2/${smallBase}_FDR_0.01.bed

awk '\$10 >= ${lfc} {print \$0}' epic2/${smallBase}_FDR_0.01.bed > epic2/${smallBase}_log_1.bed

awk '\$7 > ${cnts} {print \$0}' epic2/${smallBase}_log_1.bed > epic2/${smallBase}.bed

rm epic2/${smallBase}_FDR_0.01.bed
rm epic2/${smallBase}_log_1.bed

conda activate alignment

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R epic2/${smallBase}.bed \
	-S normalized_bw/${smallBase}_normalized.bw  \
	-o $folder/heatmaps/${smallBase}_center.gz --missingDataAsZero -p max

plotHeatmap -m $folder/heatmaps/${smallBase}_center.gz \
	-out $folder/heatmaps/${smallBase}_center.pdf --colorList 'white,darkred'

rm $folder/heatmaps/${smallBase}_center.gz

echo "${smallBase} completed!" >> ${folder}/'meta.txt'

currLine=\$(wc -l < ${folder}/meta.txt)
if ((\$currLine == $count)); then
    source activate base
    cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/homerMotif.sh ${folder}
    sh homerMotif.sh epic2
    cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/ChIPseeker.sh ${folder}
    sh ChIPseeker.sh epic2
    rm ${folder}/meta.txt
fi
EOF
        sbatch ${folder}/PBS/${smallBase}_epic2'.pbs'
        fi
    done
done




