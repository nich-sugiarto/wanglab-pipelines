#!/bin/bash

# USAGE: 
# This pipeline takes in one positional argument:
# 	$1 - target folder

# Annotates bed files with the nearest gene, distance to said gene (if applicable), functional region, and other relevant information
# Done using Homer's annotatePeaks function

folder=$(cd "$(dirname "$0")";pwd)  # Save current folder as a variable

#  If a name is not provided
if [ -z "$1" ]; then 
  echo ERROR: TARGET FOLDER WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target folder
  exit 1
fi

# Setup - create required folders
mkdir -p PBS
mkdir -p log
mkdir -p geneLists

cd $1

for file in *.bed; do
	base=$(basename "$file" ".bed")

	cat >${folder}/PBS/${base}_pileup'.pbs' <<EOF
#!/bin/bash
#SBATCH --job-name=${base}_alignment
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH -o ${folder}/log/${base}_pileup_%j.txt -e ${folder}/log/${base}_pileup_%j.err.txt
cd ${folder}
################################
# Enter your code to run below #
################################

annotatePeaks.pl ${folder}/$1/${base}.bed hg38 > ${folder}/geneLists/${base}.txt

cat ${folder}/geneLists/${base}.txt | \
	awk -F '\t' -v OFS='\t' ' \$11 != "NA" \
	{ print \$1 "	" \$2 "	" \$3 "	" \$4 "	" \$5 "	" \$6 "	" \$7 "	" \$8 "	" \$9 "	" \$10 "	" \$11 "	" \$16}' \
	> ${folder}/geneLists/${base}_noBlanks.txt
EOF
	sbatch ${folder}/PBS/${base}_pileup'.pbs'
done