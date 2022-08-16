#!/bin/bash

# ADD COMMENT DESCRIPTION HERE

# Version Doc: 

# Script automatically generated using the "slouch" command on Tue Aug  9 16:36:50 EDT 2022

folder=$(cd "$(dirname "$0")";pwd)  # Stores current folder as a variable

# Set up required files
mkdir -p PBS
mkdir -p log
mkdir -p EChO_results

cd aligned
for file in *.sorted.filtered.bam; do
	base=$(basename "$file" ".sorted.filtered.bam")
	mkdir -p ${folder}/EChO_results/${base}
	cat >${folder}/PBS/${base}_EChO'.sbatch' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=EChO # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=16GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=5:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/${base}_EChO%j.out
#SBATCH --error=${folder}/log/${base}_EChO%j.err
################################
# Enter your code to run below #
################################

cd ${folder}

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate alignment

# samtools sort -n aligned/${file} -o EChO_results/${base}/${base}.nSorted.bam

# bedtools bamtobed -bedpe -i EChO_results/${base}/${base}.nSorted.bam > ${folder}/EChO_results/${base}/${base}_frags.bed
# cat ${folder}/EChO_results/${base}/${base}_frags.bed | awk -F '\t' -v OFS='\t' ' \$1 == \$4 && ((\$2 - \$6) > -1000  && (\$2 - \$6) < 1000) { print \$1, \$2, \$6 }' > ${folder}/EChO_results/${base}/${base}_frags2.bed

/dartfs-hpc/rc/lab/W/WangX/EChO/EChO_1.0.sh epic2/${base}.bed EChO_results/${base}/${base}_frags2.bed foci EChO_results/${base}/${base}_foci
EOF
sbatch ${folder}/PBS/${base}_EChO'.sbatch'
done