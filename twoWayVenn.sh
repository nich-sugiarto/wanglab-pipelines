#!/bin/bash

# Sample body for the meta text
# markA	markB

mkdir -p PBS
mkdir -p log
mkdir -p eulerVenn

if [ -z "$1" ]; then 
  echo ERROR: META FILE WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target meta file containing sample information
  echo Meta file shouled be tab-delimited, and formatted as follows:
  echo EXAMPLE:
  echo markA.bed	markB.bed
  echo Note that each should be a RELATIVE path
  exit 1
fi

folder=$(cd "$(dirname "$0")";pwd)  # Save the current folder's location
libraryText=$1  # The passed in file

OLDIFS=$IFS
while IFS=$'\t' read -r -a varArray; do
	markA=${varArray[0]}  # Top row
	markB=${varArray[1]}  # Bottom row
	IFS=$OLDIFS

cA=$(basename "$markA" ".bed")  # Nice name for first mark
cB=$(basename "$markB" ".bed")  # Nice name for second mark

cat >${folder}/PBS/${cA}_${cB}_venn.pbs <<EOF
#!/bin/bash
#SBATCH --job-name=${cA}_${cB}_venn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${cA}_${cB}_venn_%j.txt -e ${folder}/log/${cA}_${cB}_venn_%j.err.txt
cd ${folder}
#------- End of header -------#
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate alignment

# A only
A=$(bedtools intersect \
	-v -a ${markA} \
	-b ${markB} \
	| wc -l)

# B only
B=$(bedtools intersect \
	-v -a ${markB} \
	-b ${markA} \
	| wc -l)

# A&B
AB=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markB} \
	| wc -l)

echo \$A \$B
echo \$AB

cat >${folder}/PBS/${cA}_${cB}_venn.R <<EOF2
set.seed(1) 
library(eulerr)

combo <- c("$cA" = \$A, 
	"$cB" = \$B, 
	"${cA}&${cB}" = \${AB})

pdf(file = "./eulerVenn/${cA}_${cB}_venn.pdf")
plot(euler(combo, shape = "ellipse"), quantities = TRUE)
dev.off()

EOF2

source activate base
conda info --envs
source activate /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/envs/ChIPseeker

Rscript ${folder}/PBS/${cA}_${cB}_venn.R
EOF

sbatch ${folder}/PBS/${cA}_${cB}_venn.pbs
done < ${libraryText}