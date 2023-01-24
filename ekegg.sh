if [ -z "$1" ]; then 
  echo ERROR: TARGET FOLDER WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target folder
  exit 1
fi

folder=$(cd "$(dirname "$0")";pwd)

for file in $1/*.bed; do
  base=$(basename "$file" ".bed")
  mkdir -p ${folder}/ChIPseeker/${base}
  cat > ${folder}/PBS/${base}_ChIPseeker.r <<EOF

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(annotables)
library(org.Hs.eg.db)
library(DOSE)
library(dplyr)
library(R.utils)

R.utils::setOption("clusterProfiler.download.method","wget")
ekegg <- enrichKEGG(gene = entrezids,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
pdf(file = "$folder/ChIPseeker/${base}/${base}_KEGG.pdf")
dotplot(ekegg,showCategory = 20,font.size=6)
dev.off()

do = enrichDO(entrezids)
pdf(file = "$folder/ChIPseeker/${base}/${base}_DO.pdf")
dotplot(do, showCategory=20,font.size=6)
dev.off()

EOF

  cat >${folder}/PBS/${base}_ChIPseeker'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=ChIPSeeker # Name of the job

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
#SBATCH --output=${folder}/log/%x.%j.out
#SBATCH --error=${folder}/log/%x.%j.err
################################
# Enter your code to run below #
################################
cd ${folder}

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate ChIPseeker

Rscript ${folder}/PBS/${base}_ChIPseeker.r
EOF
  sbatch ${folder}/PBS/${base}_ChIPseeker.pbs
done