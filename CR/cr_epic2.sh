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

# Default parameter cutoffs
fdr=0.01
lfc=1
cnts=20

suffix=.sorted.filtered.bam

# If nothing is provided, then use the defaults
if [ -z "$1" ]; then 
  echo Using the default parameters:
  echo "For future reference, this program takes in three parameters:"
  echo "\$1 - FDR cutoff (Default <= 0.01)"
  echo "\$2 - log fold change cutoff (Default >= 1)"
  echo "\$3 - counts cutoff (Default 20)"
elif [ -z "$3" ]; then
  echo "Error. Either pass in all three parameters, or pass in none to use the defaults"
  echo "Exiting now..."
  exit 1
else
  fdr=$1
  lfc=$2
  cnts=$3
fi

if [ -z "$4" ]; then
    echo "Proceeding with ChIPSeeker and HOMER analysis after aligment and peakcalling..."
    downstream=true
elif [ "$4" == "silence" ]; then
    echo Warning: ChIPSeeker and HOMER will NOT be run after alignment and peakcalling...
    downstream=false
else
    echo "Error! Invalid option specified. Either pass in the word "silence", or do not pass in anything at all"
    echo "Exiting now..."
    exit 1
fi

# Setup: Create needed folders
mkdir -p PBS
mkdir -p log
mkdir -p epic2
mkdir -p heatmaps

folder=$(pwd)

count=$(find ./aligned -mindepth 1 -type f -name "*${suffix}" -printf x | wc -c)  # Finds total number of files matching extension.

echo "${count} peakcalling to be done..."

# Meta file to know when qc will begin (empty file)
cat >${folder}/'meta.txt' <<EOF
EOF

cd aligned
for file in *_IgG${suffix}; do
    igGbase=$(basename "$file" "${suffix}")
    groupprefix=${igGbase%%_IgG}

     cat >${folder}/PBS/${igGbase}_epic2'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=epic2_${igGbase}  # Name of the job

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
#SBATCH --output=${folder}/log/${igGbase}_epic2.%j.out
#SBATCH --error=${folder}/log/${igGbase}_epic2.%j.err
################################
# Enter your code to run below #
################################

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate peakcalling

cd ${folder}

epic2 \
  --treatment $folder/aligned/${igGbase}${suffix} \
  --genome hg38 \
  -o epic2/${igGbase}.bed

awk '\$9 <= ${fdr} {print \$0}' epic2/${igGbase}.bed > epic2/${igGbase}_FDR_${fdr}.bed

awk '\$10 >= ${lfc} {print \$0}' epic2/${igGbase}_FDR_${fdr}.bed > epic2/${igGbase}_log_${lfc}.bed

awk '\$7 > ${cnts} {print \$0}' epic2/${igGbase}_log_${lfc}.bed > epic2/${igGbase}.bed

rm epic2/${igGbase}_FDR_${fdr}.bed
rm epic2/${igGbase}_log_${lfc}.bed

conda activate alignment

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R epic2/${igGbase}.bed \
	-S normalized_bw/${igGbase}_normalized.bw  \
	-o $igGbase/heatmaps/${igGbase}_center.gz --missingDataAsZero -p max

plotHeatmap -m $folder/heatmaps/${igGbase}_center.gz \
	-out $folder/heatmaps/${igGbase}_center.pdf --colorList 'white,darkred'

rm $folder/heatmaps/${igGbase}_center.gz

echo "${igGbase} completed!" >> ${folder}/'meta.txt'

currLine=\$(wc -l < ${folder}/meta.txt)
if ((\$currLine == $count)); then
    source activate base
    if $downstream; then
      sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/homerMotif.sh epic2
      sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/ChIPseeker.sh epic2
    fi
    rm ${folder}/meta.txt
fi
EOF
    for f in $folder/aligned/${groupprefix}*${suffix}; do
        base=$(basename "$f" "${suffix}")
        if [ "$base" != "$igGbase" ]; then
            cat >${folder}/PBS/${base}_epic2'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=epic2_${base}  # Name of the job

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
#SBATCH --output=${folder}/log/${base}_epic2.%j.out
#SBATCH --error=${folder}/log/${base}_epic2.%j.err
################################
# Enter your code to run below #
################################

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate peakcalling

cd ${folder}

epic2 \
  --treatment $folder/aligned/${base}${suffix} \
  --control $folder/aligned/${igGbase}${suffix} \
  --genome hg38 \
  -o epic2/${base}.bed

awk '\$9 <= ${fdr} {print \$0}' epic2/${base}.bed > epic2/${base}_FDR_${fdr}.bed

awk '\$10 >= ${lfc} {print \$0}' epic2/${base}_FDR_${fdr}.bed > epic2/${base}_log_${lfc}.bed

awk '\$7 > ${cnts} {print \$0}' epic2/${base}_log_${lfc}.bed > epic2/${base}.bed

rm epic2/${base}_FDR_${fdr}.bed
rm epic2/${base}_log_${lfc}.bed

conda activate alignment

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R epic2/${base}.bed \
	-S normalized_bw/${base}_normalized.bw  \
	-o $folder/heatmaps/${base}_center.gz --missingDataAsZero -p max

plotHeatmap -m $folder/heatmaps/${base}_center.gz \
	-out $folder/heatmaps/${base}_center.pdf --colorList 'white,darkred'

rm $folder/heatmaps/${base}_center.gz

echo "${base} completed!" >> ${folder}/'meta.txt'

currLine=\$(wc -l < ${folder}/meta.txt)
if ((\$currLine == $count)); then
    source activate base
    if $downstream; then
      sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/homerMotif.sh epic2
      sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/ChIPseeker.sh epic2
    fi
    rm ${folder}/meta.txt
fi
EOF
        sbatch ${folder}/PBS/${base}_epic2'.pbs'
        fi
    done
done




