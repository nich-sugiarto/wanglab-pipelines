#!/bin/bash

# USAGE: 
# This pipeline takes in one positional argument:
# 	$1 - target folder

# This script looks at all the bed files located in the user's desired folder,
# and runs HOMER's findMotifs program offshore on the discovery HPC. 

# Requires that HOMER is available in you PATH variable

mkdir -p PBS
mkdir -p log

folder=$(cd "$(dirname "$0")";pwd)

# If a name is not provided
if [ -z "$1" ]; then 
  echo ERROR: TARGET FOLDER WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target folder
  exit 1
fi

cd $1

for file in *.bed; do
	base=$(basename "$file" ".bed")
	mkdir -p ${folder}/homer_sizegiven/${base}

	cat >${folder}/PBS/${base}_motif_sizegiven'.pbs' <<EOF
#!/bin/bash -l

#SBATCH --job-name=${base}_motifs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=16:00:00
#SBATCH -o ${folder}/log/${base}_motifs_%j.txt -e ${folder}/log/${base}_motifs_%j.err.txt
# Change to job working directory
cd ${folder}
#---------------------End of header---------------------#

findMotifsGenome.pl $1/${base}.bed hg38 homer_sizegiven/${base}/ -size given
EOF

	sbatch ${folder}/PBS/${base}_motif_sizegiven.pbs
done
