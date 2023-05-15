#!/bin/bash

# Alignment pipeline for CUT&RUN data. Modified from existing pipeline 6/27/2022
# Assumes that all fastq files are contained within the same folder, and have suffix _R1_001.fastq.gz.

# Requires the "alignment" environment, and the "qc" environment (for fastqc and multiqc)

# TODO: Keep deduplicated for fastqc?
# TODO: Test pipeline

# Setup: Create needed folders
mkdir -p PBS
mkdir -p deduplicated
mkdir -p trimmed
mkdir -p aligned
mkdir -p normalized_bw
mkdir -p log

if [ -z "$1" ]; then
    echo "Proceeding with ChIPSeeker and HOMER analysis after aligment and peakcalling..."
    downstream=true
elif [ "$1" == "silence" ]; then
    echo Warning: ChIPSeeker and HOMER will NOT be run after alignment and peakcalling...
    downstream=false
else
    echo "Error! Invalid option specified. Either pass in the word "silence", or do not pass in anything at all"
    echo "Exiting now..."
    exit 1
fi

suffix1=_R1_001.fastq.gz

folder=$(pwd)  # Saves folder as a variable

count=$(find ./fastq -mindepth 1 -type f -name "*${suffix1}" -printf x | wc -c)  # Finds total number of files matching extension.
echo $count files found!

# Error handling for if there are 0 files
if (($count == 0)); then
	echo Warning! There were no files that were found.
	echo Check to ensure that the suffix is correct. Otherwise, it may be that each fastq file is contained within its own subfolder.
	echo Would you like to expand all subfolders into this one? [y/n]:
	
	# Prompt user to whether they want to try and expand out the subfolders
	read resp

	if [[ "$resp" == "y" ]]; then
		echo Expanding out...
		cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/move.sh ./
		sh move.sh
		echo Done. Proceeding as normal...
		rm move.sh
	elif [[ "$resp" == "n" ]]; then
		echo Terminating run...
		exit 1
	else
		echo Error! Invalid input. Terminating run...
		exit 1
	fi
fi

# Meta file to know when qc will begin (empty file)
cat >${folder}/'meta.txt' <<EOF
EOF

cd fastq

# Loop over all the fastq files
for file in *${suffix1}; do
    base=$(basename "$file" "${suffix1}")  
    smallBase=${base%_S*}  
    echo ${smallBase}
    cat >${folder}/PBS/${smallBase}'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${smallBase}_alignment  # Name of the job

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
#SBATCH --output=${folder}/log/${smallBase}_alignment.%j.out
#SBATCH --error=${folder}/log/${smallBase}_alignment.%j.err
################################
# Enter your code to run below #
################################
cd ${folder}
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate alignment

clumpify.sh \
    in1=fastq/${file} \
    in2=fastq/${base}${suffix1/R1/R2} \
    out1=deduplicated/${smallBase}_R1_dedup.fastq.gz \
    out2=deduplicated/${smallBase}_R2_dedup.fastq.gz \
    dedupe subs=2  

gunzip deduplicated/${smallBase}_R1_dedup.fastq.gz
gunzip deduplicated/${smallBase}_R2_dedup.fastq.gz

bbduk.sh \
    in1=deduplicated/${smallBase}_R1_dedup.fastq \
    in2=deduplicated/${smallBase}_R2_dedup.fastq \
    out1=trimmed/${smallBase}_R1_trimmed.fastq \
    out2=trimmed/${smallBase}_R2_trimmed.fastq \
    ref=/dartfs-hpc/rc/lab/W/WangX/Nicholas/bbmap/resources/adapters.fa \
    ktrim=r k=21 mink=11 hdist=1 tpe tbo

rm deduplicated/${smallBase}_R1_dedup.fastq*
rm deduplicated/${smallBase}_R2_dedup.fastq*

module load bowtie/2.2.7

bowtie2 -x \
    /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/GRCh38/grch38_1kgmaj \
    -1 trimmed/${smallBase}_R1_trimmed.fastq \
    -2 trimmed/${smallBase}_R2_trimmed.fastq \
    --local --very-sensitive-local --no-unal --no-mixed \
    --no-discordant --phred33 -I 10 -X 700 > aligned/${smallBase}.sam

rm trimmed/${smallBase}_R1_trimmed.fastq
rm trimmed/${smallBase}_R2_trimmed.fastq

samtools view -Sbo aligned/${smallBase}.bam aligned/${smallBase}.sam
samtools sort aligned/${smallBase}.bam -o aligned/${smallBase}.sorted.bam

samtools index aligned/${smallBase}.sorted.bam

bamCoverage \
    --bam aligned/${smallBase}.sorted.bam \
    -o normalized_bw/${smallBase}_normalized.bw \
    --binSize 1 --normalizeUsing RPKM --ignoreDuplicates \
    --extendReads -of bigwig -p max \
    -bl /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/GRCh38/hg38-blacklist.v2.bed \
    --ignoreForNormalization chrX chrM chrRandom chrUn

samtools view -b \
    aligned/${smallBase}.sorted.bam \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
    chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
    > aligned/${smallBase}.chr_filt.bam

bedtools intersect -v -a aligned/${smallBase}.chr_filt.bam -b \
    /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/GRCh38/hg38-blacklist.v2.bed > aligned/${smallBase}.sorted.filtered.bam

samtools index aligned/${smallBase}.sorted.filtered.bam

rm ${folder}/aligned/${smallBase}.chr_filt.bam
rm ${folder}/aligned/${smallBase}.sam
rm ${folder}/aligned/${smallBase}.sorted.bam
rm ${folder}/aligned/${smallBase}.sorted.bam.bai

echo "${smallBase} completed!" >> ${folder}/'meta.txt'

# Checks to see if all the files have been analyzes
# if so, will run fastqc and multiqc on fastq and bam files 
# (requires access to Nick's pipeline folder)

currLine=\$(wc -l < ${folder}/meta.txt)
if ((\$currLine == $count)); then
    source activate base
    rmdir deduplicated/
    rmdir trimmed/
    sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/qc.sh
    if $downstream; then
        sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/CR/cr_epic2.sh
    else
        sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/CR/cr_epic2.sh 0.01 1 20 silence
    fi
    rm ${folder}/meta.txt
fi
EOF
    sbatch ${folder}/PBS/${smallBase}'.pbs'
done