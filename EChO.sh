#!/bin/bash

# Runs EChO analysis to generate a bed file of foci of enriched fragment regions.
# Calls plot_EChO.sh afterwards automatically to plot said regions

# Script automatically generated using the "slouch" command on Tue Aug  9 16:36:50 EDT 2022

folder=$(pwd)  # Stores current folder as a variable

# Set up required files
mkdir -p PBS
mkdir -p log
mkdir -p EChO_results

suffix=.sorted.filtered.bam

count=$(find ./aligned -mindepth 1 -type f -name "*${suffix}" -printf x | wc -c)  # Finds total number of files matching extension.
IgGcount=$(find ./aligned -mindepth 1 -type f -name "*IgG*${suffix}" -printf x | wc -c)  # Finds total number of IgG files matching extension.
((count-=IgGcount))

echo $count files found!
# Meta file to know when qc will begin (empty file)
cat >${folder}/'EChO_meta.txt' <<EOF
EOF

cd aligned

for file in *${suffix}; do
	base=$(basename "$file" "${suffix}")
	# Check to make sure not control
	if [[ ${base} != *"IgG"* ]]; then
		cat >${folder}/PBS/${base}_EChO'.sbatch' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${base}_EChO # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=16GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=24:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/${base}_EChO%j.out
#SBATCH --error=${folder}/log/${base}_EChO%j.err
################################
# Enter your code to run below #
################################

cd ${folder}

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate alignment

# Sorts reads by name
samtools sort -n aligned/${base}.bam -o EChO_results/${base}.nSorted.bam

# Writes alignments to bedpe format
bedtools bamtobed -bedpe -i EChO_results/${base}.nSorted.bam > ${folder}/EChO_results/${base}_frags.bed

rm aligned/${base}.bam
rm EChO_results/${base}.nSorted.bam

# Filters out reads where pairs are on different chromosomes, or are more than 1000 bp apart
cat ${folder}/EChO_results/${base}_frags.bed | \
	awk -F '\t' -v OFS='\t' ' \$1 == \$4 && ((\$2 - \$6) > -1000  && (\$2 - \$6) < 1000) { print \$1, \$2, \$6 }' \
	> ${folder}/EChO_results/${base}_filteredFrags.bed

# Calls EChO to determine foci of enriched regions
/dartfs-hpc/rc/lab/W/WangX/EChO/EChO_1.0.sh epic2/${base}.bed EChO_results/${base}_filteredFrags.bed foci EChO_results/${base}_foci

rm EChO_results/${base}_filteredFrags.bed

echo "${base} completed!" >> ${folder}/'EChO_meta.txt'

currLine=\$(wc -l < ${folder}/EChO_meta.txt.txt)
if ((\$currLine == $count)); then
    cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/plot_EChO.sh ${folder}
    sh plot_EChO.sh
    rm ${folder}/EChO_meta.txt
fi
EOF
		sbatch ${folder}/PBS/${base}_EChO'.sbatch'
	fi
done