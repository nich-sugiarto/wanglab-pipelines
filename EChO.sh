#!/bin/bash

# ADD COMMENT DESCRIPTION HERE

# Version Doc: 

# Script automatically generated using the "slouch" command on Tue Aug  9 16:36:50 EDT 2022

folder=$(cd "$(dirname "$0")";pwd)  # Stores current folder as a variable

# Set up required files
mkdir -p PBS
mkdir -p log
cd aligned
for file in *.sorted.filtered.bam; do
	base=$(basename "$file" ".sorted.filtered.bam")
	cat >${folder}/PBS/EChO'.pbs' <<EOF
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
#SBATCH --output=${folder}/log/EChO%j.out
#SBATCH --error=${folder}/log/EChO%j.err
################################
# Enter your code to run below #
################################

EOF
sbatch ${folder}/PBS/EChO'.pbs'
done


