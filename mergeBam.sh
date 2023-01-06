#!/bin/bash

# Merges bam files 

# Requires the "alignment" environment, and the "qc" environment (for fastqc and multiqc)

mkdir -p PBS
mkdir -p mergedBam
mkdir -p mergedBigWig

folder=$(cd "$(dirname "$0")";pwd)  # Saves folder as a variable

suffix1=_R1.shifted.nsorted.bam
CR_sfx=_R1.sorted.filtered.bam
count=$(find ./aligned -mindepth 1 -type f -name "*${suffix1}" -printf x | wc -c)  # Finds total number of files matching extension.
echo $count files found!

# Error handling for if there are 0 files
if (($count == 0)); then
	echo Warning! There were no files that were found. Is this CUT\&RUN data?
	echo Will proceed, but will assume that this is CUT\&RUN data... 
	
    suffix1=$CR_sfx
    count=$(find ./aligned -mindepth 1 -type f -name "*${suffix1}" -printf x | wc -c)  # Finds total number of files matching extension.
    echo $count files found!
fi

# Meta file to know when qc will begin (empty file)
cat >${folder}/'mergeMeta.txt' <<EOF
EOF

cd aligned
# Loop over all the fastq files
for file in *${suffix1}; do
	base=$(basename "$file" "${suffix1}")  
	echo ${base}
	cat >${folder}/PBS/${base}'_mergeBam.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${base}_mergeBam # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=32GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=10:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/${base}_mergeBam.%j.out
#SBATCH --error=${folder}/log/${base}_mergeBam.%j.err
################################
# Enter your code to run below #
################################
cd ${folder}
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate alignment

samtools merge -f mergedBam/${base}_merged.bam aligned/${file} aligned/${base}${suffix1/R1/R2} 
samtools index mergedBam/${base}_merged.bam

samtools sort mergedBam/${base}_merged.bam -o mergedBam/${base}_sortedMerged.bam
samtools index mergedBam/${base}_sortedMerged.bam

bamCoverage \
    --bam mergedBam/${base}_sortedMerged.bam \
    -o mergedBigWig/${base}_merged.bw \
    --binSize 1 --normalizeUsing RPKM --ignoreDuplicates \
    --extendReads -of bigwig -p max \
    -bl /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/GRCh38/hg38-blacklist.v2.bed \
    --ignoreForNormalization chrX chrM chrRandom chrUn

rm mergedBam/${base}_merged.bam
rm mergedBam/${base}_merged.bam.bai

echo "${base} completed!" >> ${folder}/'mergeMeta.txt'

currLine=\$(wc -l < ${folder}/mergeMeta.txt)
if ((\$currLine == $count)); then
    source activate base
    if [[ $suffix1 == "_R1.shifted.nsorted.bam" ]]; then
        cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/mergeATACpeaks.sh ${folder}
        sh mergePeakCalling.sh
    elif [[ $suffix1 == "_R1.sorted.filtered.bam" ]]; then
        cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/mergePeakCalling.sh ${folder}
        sh mergePeakCalling.sh
    else
        echo "Something went horribly wrong"
    fi
    rm ${folder}/mergeMeta.txt
fi
EOF
sbatch ${folder}/PBS/${base}'_mergeBam.pbs'
done