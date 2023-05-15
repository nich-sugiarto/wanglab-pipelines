#!/bin/bash

# Plots the EChO results from EChO.sh. Should automatically trigger

# Script automatically generated using the "slouch" command on Thu Aug 18 16:21:35 EDT 2022

folder=$(pwd)  # Stores current folder as a variable

# Set up required files
mkdir -p PBS
mkdir -p log

# For each sample group (determined by IgG control)
for file in aligned/*_IgG*.sorted.filtered.bam; do
    groupPrefix=${file%%_IgG*.sorted.filtered.bam}  # 
	groupPrefix=${groupPrefix#"aligned/"}
    echo ${groupPrefix}
    cat >${folder}/PBS/${groupPrefix}_EChO'.R' <<EOF
library(ggplot2)

getwd()
df <- data.frame(matrix(ncol = 2, nrow = 0))
EOF
    # For each sample in that group
    for f in ./EChO_results/${groupPrefix}*.EChO.bed; do
        base=$(basename "$f" "_foci.EChO.bed")
        base=${base#"./EchO_results/CR_O422T"}
        echo ${base}
            # Add that file to the overall dataframe
            cat >>${folder}/PBS/${groupPrefix}_EChO'.R' <<EOF

df2 <- read.table("$f", sep = "\t")
df2 <- data.frame(sampleGroup="${base}", fragSize = df2\$V4)
str(df2)

df <- rbind(df, df2)
EOF
    done

    # Create the graph
    cat >>${folder}/PBS/${groupPrefix}_EChO'.R' <<EOF

png("EChO_results/${groupPrefix}.png", width = 1080, height = 720)
    ggplot(df, aes(x=sampleGroup, y=fragSize, fill = sampleGroup)) + geom_violin() + ggtitle("${groupPrefix} Fragment Size Distributions") + geom_boxplot(width = 0.2) + geom_hline(aes(yintercept=147))
dev.off()
EOF

    # Actually run the scripts
    cat >${folder}/PBS/${groupPrefix}_plotEChO'.sbatch' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${groupPrefix}_plotEChO  # Name of the job

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
#SBATCH --output=${folder}/log/${groupPrefix}_plotEChO.%j.out
#SBATCH --error=${folder}/log/${groupPrefix}_plotEChO.%j.err
################################
# Enter your code to run below #
################################
cd ${folder}

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate ATACQC

Rscript ${folder}/PBS/${groupPrefix}_EChO'.R'

rm EChO_results/*${base}*frags*
EOF
    sbatch ${folder}/PBS/${groupPrefix}_plotEChO'.sbatch'
done