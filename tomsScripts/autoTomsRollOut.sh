#!/bin/bash
folder=$(cd "$(dirname "$0")";pwd)  # Stores current folder as a variable

# Set up required files
mkdir -p PBS
mkdir -p log
mkdir -p finalCounts;

for file in *.fastq; do
	base=$(basename "${file}" "_001.fastq")
	cat >${folder}/finalCounts/${base}_LevenshteinCounts.txt <<EOF1
EOF1
	cat >${folder}/finalCounts/${base}_filt.fastq <<EOF1
EOF1

	cat >${folder}/PBS/${base}_levenshtein'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${base}_levenshtein # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=8GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=192:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/${base}_levenshtein%j.out
#SBATCH --error=${folder}/log/${base}_levenshtein%j.err
################################
# Enter your code to run below #
################################
cd ${folder}

cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/tomsScripts/strictTomsScript.sh ${folder}
sh strictTomsScript.sh ${file} finalCounts/${base}_LevenshteinCounts.txt finalCounts/${base}_filt.fastq
EOF

	sbatch ${folder}/PBS/${base}_levenshtein'.pbs'
done