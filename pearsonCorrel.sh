#!/bin/bash

# Calculates the Pearson correlation between reps to determine consistency

# Setup: Create needed folders
mkdir -p PBS
mkdir -p log
mkdir -p pearsonCorrel

folder=$(pwd)  # Saves folder as a variable

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

cd aligned

# Loop over all the fastq files
for file in *${suffix1}; do
    base=$(basename "$file" "${suffix1}")  
    echo ${base}
    cat >${folder}/PBS/${base}'_pearsonCorrel.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${base}_pearsonCorrel  # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=8GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=10:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/${base}_pearsonCorrel.%j.out
#SBATCH --error=${folder}/log/${base}_pearsonCorrel.%j.err
################################
# Enter your code to run below #
################################
cd ${folder}
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate alignment

multiBamSummary bins \
	--bamfiles aligned/${file} aligned/${base}${suffix1/1/2} \
	 -o pearsonCorrel/${base}.npz \
	 --labels ${base}1 ${base}2 \
	 -bl /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/GRCh38/hg38-blacklist.v2.bed \


plotCorrelation \
    -in pearsonCorrel/${base}.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation of Average Scores Per Transcript" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o pearsonCorrel/${base}_heatmap.pdf

plotCorrelation \
	-in pearsonCorrel/${base}.npz \
	--corMethod pearson --skipZeros \
	--plotTitle "Pearson Correlation of Average Scores Per Transcript" \
	--whatToPlot scatterplot \
	-o pearsonCorrel/${base}_scatter.pdf

EOF
    sbatch ${folder}/PBS/${base}'_pearsonCorrel.pbs'
done