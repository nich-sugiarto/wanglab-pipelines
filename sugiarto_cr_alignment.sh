#!/bin/bash

# Alignment pipeline for CUT&RUN data. Modified from existing pipeline 6/27/2022
# Assumes that all fastq files are contained within the same folder, and have suffix _1.fq.gz.
# See areas marked "CHANGE FOR FILE EXTENSION"
# Stolen from Ruoyun Wang

# Requires the "vanilla" environment, and the "qc" environment (for fastqc and multiqc)

# Version Doc: https://docs.google.com/document/d/1LL0hB4eDcqZ7_UmBX_KwyiYlm5r1fBio-Y4hdBUdZQw/edit

# TODO: Use chromap instead of bowtie2 for faster and more accurate(?) alignment
# TODO: Keep deduplicated for fastqc?
# TODO: Test pipeline
# TODO: Add peakcalling and downstream analysis (ChIPSeeker, Motif Analysis) to run automatically with it

# Setup: Create needed folders
mkdir -p PBS
mkdir -p deduplicated
mkdir -p trimmed
mkdir -p aligned
mkdir -p normalized_bw
mkdir -p log

folder=$(cd "$(dirname "$0")";pwd)  # Saves folder as a variable

count=$(find ./ -mindepth 1 -type f -name "*_1.fq.gz" -printf x | wc -c)  # Finds total number of files matching extension. Needed to know when to start qc. CHANGE FOR FILE EXTENSION

# Meta file to know when qc will begin (empty file)
cat >${folder}/'meta.txt' <<EOF
EOF

# Loop over all the fastq files
# CHANGE FOR FILE EXTENSION
for file in *_1.fq.gz; do
    cd ${folder}/fastq

    base=$(basename "$file" "_1.fq.gz")  # CHANGE FOR FILE EXTENSION
    smallBase=${base%_CKDL*}  # CHANGE FOR FILE EXTENSION

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
#SBATCH --output=${smallBase}_alignment%x.%j.out
#SBATCH --error=${smallBase}_alignment%x.%j.err
################################
# Enter your code to run below #
################################

source activate vanilla

clumpify.sh \
    in1=fastq/${file} \
    in2=fastq/${base}_2.fq.gz \
    out1=deduplicated/${smallBase}_R1_dedup.fastq.gz \
    out2=deduplicated/${smallBase}_R2_dedup.fastq.gz \
    dedupe subs=2  # CHANGE FOR FILE EXTENSION

gunzip deduplicated/${smallBase}_R1_dedup.fastq.gz
gunzip deduplicated/${smallBase}_R2_dedup.fastq.gz

bbduk.sh \
    in1=deduplicated/${smallBase}_R1_dedup.fastq \
    in2=deduplicated/${smallBase}_R2_dedup.fastq \
    out1=trimmed/${smallBase}_R1_trimmed.fastq \
    out2=trimmed/${smallBase}_R2_trimmed.fastq \
    ref=/dartfs-hpc/rc/lab/W/WangX/Nicholas/bbmap/resources/adapters.fa \
    ktrim=r k=21 mink=11 hdist=1 tpe tbo

module load bowtie/2.2.7

bowtie2 -x \
    /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/GRCh38/grch38_1kgmaj \
    -1 trimmed/${smallBase}_R1_trimmed.fastq \
    -2 trimmed/${smallBase}_R2_trimmed.fastq \
    --local --very-sensitive-local --no-unal --no-mixed \
    --no-discordant --phred33 -I 10 -X 700 > aligned/${smallBase}.sam

samtools view -Sbo aligned/${smallBase}.bam aligned/${smallBase}.sam
samtools sort aligned/${smallBase}.bam -o aligned/${smallBase}sorted.bam

samtools index aligned/${smallBase}sorted.bam

bamCoverage \
    --bam aligned/${smallBase}sorted.bam \
    -o normalized_bw/${smallBase}_normalized.bw \
    --binSize 1 --normalizeUsing RPKM --ignoreDuplicates \
    --extendReads -of bigwig -p max \
    -bl /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/GRCh38/hg38-blacklist.v2.bed \
    --ignoreForNormalization chrX chrM chrRandom chrUn

samtools view -b \
    aligned/${smallBase}sorted.bam \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
    chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
    > aligned/${smallBase}_chr_filt.bam

bedtools intersect -v -a aligned/${smallBase}_chr_filt.bam -b \
    /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/GRCh38/hg38-blacklist.v2.bed > aligned/${smallBase}_sorted_filtered.bam

samtools index aligned/${smallBase}_sorted_filtered.bam

rm ${folder}/aligned/${smallBase}.bam
rm ${folder}/aligned/${smallBase}_chr_filt.bam
rm ${folder}/aligned/${smallBase}.sam
rm ${folder}/aligned/${smallBase}sorted.bam
rm ${folder}/aligned/${smallBase}sorted.bam.bai

echo "${smallBase} completed!" >> ${folder}/'meta.txt'

# Checks to see if all the files have been analyzes
# if so, will run fastqc and multiqc on fastq and bam files 
# (requires access to Nick's pipeline folder)

currLine=\$(wc -l < ${folder}/meta.txt)
if ((\$currLine == (($count + 1)))); then
    rm -r deduplicated/
    rm -r trimmed/
    cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/sugiarto_qc.sh ${folder}
    sh sugiarto_qc.sh
    rm ${folder}/meta.txt
fi
EOF
    cd ${folder}/log
    sbatch ${folder}/PBS/${smallBase}'.pbs'
done