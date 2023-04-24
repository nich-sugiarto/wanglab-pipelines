#!/bin/bash

# This pipeline requires 1 positional arguments:
# 	$1 - The Factor (ARID1A, ARID1B, etc.)
folder=$(cd "$(dirname "$0")";pwd)
mkdir -p diffBind_Heatmaps/Clusters
factor=$1
if [ -z "$1" ]; then 
  echo ERROR: NOT ALL PARAMETERS HAVE BEEN SPECIFIED
  echo USAGE:
  echo This pipeline takes in a minimum of 1 positional argument:
  echo '$1 - The Factor (ARID1A, ARID1B, etc.)'
  exit 1
fi

cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/downstreamDiffBind/diffBindHeatmap.R ${folder}/PBS

cat >${folder}/PBS/diffBind_Heatmaps.sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=diffBind_Heatmap
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=7:00:00
#SBATCH -o log/diffBind_Heatmap_%j.txt -e log/diffBind_Heatmap_%j.err.txt
#------- END OF HEADER -------#

source activate R
cd ${folder}
mkdir -p diffBind_Heatmaps
Rscript PBS/diffBindHeatmap.R ${factor}
EOF
sbatch PBS/diffBind_Heatmaps.sbatch
