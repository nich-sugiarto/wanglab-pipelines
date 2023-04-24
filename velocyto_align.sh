#!/bin/bash

# Realigns 10x sample to fit what velocyto needs to do pseudotime

# Script automatically generated using the "slouch" command on Tue Aug  9 11:57:02 EDT 2022

if [ -z "$1" ]; then 
  echo ERROR: TARGET FOLDER WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target folder
  exit 1
fi

target=$1
folder=$(pwd)  # Stores current folder as a variable

# Setup - Generate required folders
mkdir -p PBS
mkdir -p log

cat >${folder}/PBS/velocyto_align'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=velocyto_align # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=128GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=24:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/velocyto_align%j.out
#SBATCH --error=${folder}/log/velocyto_align%j.err
################################
# Enter your code to run below #
################################
source activate velocyto

cd ${folder}

velocyto run10x \
	-m /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/hg38_repeat_rmsk.gtf \
	${folder}/${1} \
	/dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/refdata-gex-GRCh38-2020-A/genes/genes.gtf

EOF

sbatch ${folder}/PBS/velocyto_align'.pbs'
