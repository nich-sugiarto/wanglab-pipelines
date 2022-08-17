#!/bin/bash

# Fastqc and multiqc data for fastq files. Created 6/27/2022
# Assumes that all fastq files are contained within the same folder

# Requires the "qc" environment

# Vesion control: https://docs.google.com/document/d/1VjdQUg3D2mrZnQEgmvFfvegX_VJpGF-0J2Td-2V7ua4/edit

#TODO: Use/test only on deduplicated reads? 
#FIXME: Are the bam files working?


source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh

folder=$(cd "$(dirname "$0")";pwd)

cat >${folder}/PBS/'fastqc.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=fastqc # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=8GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=24:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/fastqc.%j.out
#SBATCH --error=${folder}/log/fastqc.%j.err
################################
# Enter your code to run below #
################################

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate qc

cd ${folder}

mkdir -p fastqc
mkdir -p fastqc/fastq
mkdir -p fastqc/bam
mkdir -p multiqc

# fastq files
fastqc -o fastqc/fastq fastq/*
multiqc -f -o multiqc -n fastq fastqc/fastq/

# bam files
fastqc -o fastqc/bam aligned/*.bam
multiqc -f -o multiqc -n bam fastqc/bam/
EOF
sbatch ${folder}/PBS/$base'fastqc.pbs'