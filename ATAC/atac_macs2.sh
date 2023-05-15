#!/bin/bash

# Alternative to Genrich alignment using macs2 instead

# Script automatically generated using the "slouch" command on Thu Sep  8 16:13:51 EDT 2022

folder=$(pwd)  # Stores current folder as a variable

# Set up required files
mkdir -p PBS
mkdir -p log
mkdir -p macs2
cd aligned

for file in *.shifted.nsorted.bam; do
	base=$(basename "${file}" ".shifted.nsorted.bam")
	cat >${folder}/PBS/${base}_atac_macs2'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${base}_atac_macs2 # Name of the job

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
#SBATCH --output=${folder}/log/${base}_atac_macs2%j.out
#SBATCH --error=${folder}/log/${base}_atac_macs2%j.err
################################
# Enter your code to run below #
################################
cd ${folder}

conda activate alignment

macs2 callpeak -f BAMPE -t ${folder}/aligned/${file} -n ${folder}/macs2/${base}
EOF

	sbatch ${folder}/PBS/${base}_atac_macs2'.pbs'
done