#!/bin/bash

# Creates a five-way venn diagram that is quite frankly unholy. I do not like this script, and yet I cannot find an easier way to write it.
# TODO: Maybe I can write it using a nested for loop?

# Sample body for the meta text
# markA	markB	markC	markD	markE

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
  echo markA.bed	markB.bed	markC.bed	markD.bed	markE.bed
  echo Note that each should be a RELATIVE path
  exit 1
fi

# conda activate alignment

folder=$(pwd)
libraryText=$1

OLDIFS=$IFS
while IFS=$'\t' read -r -a varArray; do
	markA=${varArray[0]}  # Top row
	markB=${varArray[1]}  # Bottom row
	markC=${varArray[2]}  # Left column
	markD=${varArray[3]}  # Right column
	markE=${varArray[4]}  # Right column
	IFS=$OLDIFS
done < ${libraryText}

cA=$(basename "$markA" ".bed")
cB=$(basename "$markB" ".bed")
cC=$(basename "$markC" ".bed")
cD=$(basename "$markD" ".bed")
cE=$(basename "$markE" ".bed")

cat >${folder}/PBS/fiveWayVenn.pbs <<EOF
#!/bin/bash
#SBATCH --job-name=fiveWayVenn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/fiveWayVenn_%j.txt -e ${folder}/log/fiveWayVenn_%j.err.txt
cd ${folder}
#------- End of header -------#
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate alignment

# A only
A=$(bedtools intersect \
	-v -a ${markA} \
	-b ${markB} ${markC} ${markD} ${markE} \
	| wc -l)

# B only
B=$(bedtools intersect \
	-v -a ${markB} \
	-b ${markA} ${markC} ${markD} ${markE} \
	| wc -l)

# C only
C=$(bedtools intersect \
	-v -a ${markC} \
	-b ${markA} ${markB} ${markD} ${markE} \
	| wc -l)

# D only
D=$(bedtools intersect \
	-v -a ${markD} \
	-b ${markA} ${markB} ${markC} ${markE} \
	| wc -l)

# E only
E=$(bedtools intersect \
	-v -a ${markE} \
	-b ${markA} ${markB} ${markC} ${markD} \
	| wc -l)

# A&B
AB=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markB} | \
	bedtools intersect \
		-v -a stdin \
		-b ${markC} ${markD} ${markE} \
		| wc -l)

# A&C
AC=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markC} | \
	bedtools intersect \
		-v -a stdin \
		-b ${markB} ${markD} ${markE} \
		| wc -l)

# A&D
AD=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markD} | \
	bedtools intersect \
		-v -a stdin \
		-b ${markB} ${markC} ${markE} \
		| wc -l)

# A&E
AE=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markE} | \
	bedtools intersect \
		-v -a stdin \
		-b ${markB} ${markC} ${markD} \
		| wc -l)

# B&C
BC=$(bedtools intersect \
	-wa -a ${markB} \
	-b ${markC} | \
	bedtools intersect \
		-v -a stdin \
		-b ${markA} ${markD} ${markE} \
		| wc -l)

# B&D
BD=$(bedtools intersect \
	-wa -a ${markB} \
	-b ${markD} | \
	bedtools intersect \
		-v -a stdin \
		-b ${markA} ${markC} ${markE} \
		| wc -l)

# B&E
BE=$(bedtools intersect \
	-wa -a ${markB} \
	-b ${markE} | \
	bedtools intersect \
		-v -a stdin \
		-b ${markA} ${markC} ${markD} \
		| wc -l)

# C&D
CD=$(bedtools intersect \
	-wa -a ${markC} \
	-b ${markD} | \
	bedtools intersect \
		-v -a stdin \
		-b ${markA} ${markB} ${markE} \
		| wc -l)

# C&E
CE=$(bedtools intersect \
	-wa -a ${markC} \
	-b ${markE} | \
	bedtools intersect \
		-v -a stdin \
		-b ${markA} ${markB} ${markD} \
		| wc -l)

# D&E
DE=$(bedtools intersect \
	-wa -a ${markD} \
	-b ${markE} | \
	bedtools intersect \
		-v -a stdin \
		-b ${markA} ${markB} ${markC} \
		| wc -l)

# A&B&C
ABC=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markB} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markC} | \
		bedtools intersect \
			-v -a stdin \
			-b ${markD} ${markE} \
			| wc -l)

# A&B&D
ABD=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markB} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markD} | \
		bedtools intersect \
			-v -a stdin \
			-b ${markC} ${markE} \
			| wc -l)

# A&B&E
ABE=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markB} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markE} | \
		bedtools intersect \
			-v -a stdin \
			-b ${markC} ${markD} \
			| wc -l)

# A&C&D
ACD=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markC} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markD} | \
		bedtools intersect \
			-v -a stdin \
			-b ${markB} ${markE} \
			| wc -l)

# A&C&E
ACE=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markC} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markE} | \
		bedtools intersect \
			-v -a stdin \
			-b ${markB} ${markD} \
			| wc -l)

# A&D&E
ADE=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markD} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markE} | \
		bedtools intersect \
			-v -a stdin \
			-b ${markB} ${markC} \
			| wc -l)

# B&C&D
BCD=$(bedtools intersect \
	-wa -a ${markB} \
	-b ${markC} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markD} | \
		bedtools intersect \
			-v -a stdin \
			-b ${markA} ${markE} \
			| wc -l)

# B&C&E
BCE=$(bedtools intersect \
	-wa -a ${markB} \
	-b ${markC} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markE} | \
		bedtools intersect \
			-v -a stdin \
			-b ${markA} ${markD} \
			| wc -l)

# B&D&E
BDE=$(bedtools intersect \
	-wa -a ${markB} \
	-b ${markD} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markE} | \
		bedtools intersect \
			-v -a stdin \
			-b ${markA} ${markC} \
			| wc -l)

# C&D&E
CDE=$(bedtools intersect \
	-wa -a ${markC} \
	-b ${markD} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markE} | \
		bedtools intersect \
			-v -a stdin \
			-b ${markA} ${markB} \
			| wc -l)

# A&B&C&D
ABCD=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markB} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markC} | \
		bedtools intersect \
			-wa -a stdin \
			-b ${markD} | \
			bedtools intersect \
				-v -a stdin \
				-b ${markE} \
				| wc -l)

# B&C&D&E
BCDE=$(bedtools intersect \
	-wa -a ${markB} \
	-b ${markC} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markD} | \
		bedtools intersect \
			-wa -a stdin \
			-b ${markE} | \
			bedtools intersect \
				-v -a stdin \
				-b ${markA} \
				| wc -l)

# C&D&E&A
CDEA=$(bedtools intersect \
	-wa -a ${markC} \
	-b ${markD} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markE} | \
		bedtools intersect \
			-wa -a stdin \
			-b ${markA} | \
			bedtools intersect \
				-v -a stdin \
				-b ${markB} \
				| wc -l)

# D&E&A&B
DEAB=$(bedtools intersect \
	-wa -a ${markD} \
	-b ${markE} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markA} | \
		bedtools intersect \
			-wa -a stdin \
			-b ${markB} | \
			bedtools intersect \
				-v -a stdin \
				-b ${markC} \
				| wc -l)

# E&A&B&C
EABC=$(bedtools intersect \
	-wa -a ${markE} \
	-b ${markA} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markB} | \
		bedtools intersect \
			-wa -a stdin \
			-b ${markC} | \
			bedtools intersect \
				-v -a stdin \
				-b ${markD} \
				| wc -l)

# A&B&C&D&E
ABCDE=$(bedtools intersect \
	-wa -a ${markA} \
	-b ${markB} | \
	bedtools intersect \
		-wa -a stdin \
		-b ${markC} | \
		bedtools intersect \
			-wa -a stdin \
			-b ${markD} | \
			bedtools intersect \
				-wa -a stdin \
				-b ${markE} \
				| wc -l)


echo \$A \$B \$C \$D \$E
echo \$AB \$AC \$AD \$AE \$BC \$BD \$BE \$CD \$CE \$DE
echo \$ABC \$ABD \$ABE \$ACD \$ACE \$ADE \$BCD \$BCE \$BCD \$CDE
echo \$ABCD \$BCDE \$CDEA \$DEAB \$EABC
echo \$ABCDE

cat >${folder}/PBS/fiveWayVenn.R <<EOF2
set.seed(1) 
library(eulerr)

combo <- c("$cA" = \$A, 
	"$cB" = \$B, 
	"$cC" = \$C, 
	"$cD" = \$D, 
	"$cE" = \$E, 
	"${cA}&${cB}" = \${AB}, 
	"${cA}&${cC}" = \${AC}, 
	"${cA}&${cD}" = \${AD}, 
	"${cA}&${cE}" = \${AE}, 
	"${cB}&${cC}" = \${BC}, 
	"${cB}&${cD}" = \${BD}, 
	"${cB}&${cE}" = \${BE}, 
	"${cC}&${cD}" = \${CD}, 
	"${cC}&${cE}" = \${CE}, 
	"${cD}&${cE}" = \${DE}, 
	"${cA}&${cB}&${cC}" = \${ABC}, 
	"${cA}&${cB}&${cD}" = \${ABD}, 
	"${cA}&${cB}&${cE}" = \${ABE}, 
	"${cA}&${cC}&${cD}" = \${ACD}, 
	"${cA}&${cC}&${cE}" = \${ACE}, 
	"${cA}&${cD}&${cE}" = \${ADE}, 
	"${cB}&${cC}&${cD}" = \${BCD}, 
	"${cB}&${cC}&${cE}" = \${BCE}, 
	"${cB}&${cD}&${cE}" = \${BDE}, 
	"${cC}&${cD}&${cE}" = \${CDE}, 
	"${cA}&${cB}&${cC}&${cD}" = \${ABCD}, 
	"${cB}&${cC}&${cD}&${cE}" = \${BCDE}, 
	"${cA}&${cC}&${cD}&${cE}" = \${CDEA}, 
	"${cA}&${cB}&${cD}&${cE}" = \${DEAB}, 
	"${cA}&${cB}&${cC}&${cE}" = \${EABC}, 
	"${cA}&${cB}&${cC}&${cD}&${cE}" = \${ABCDE})

pdf(file = "./eulerVenn/venn.pdf")
plot(euler(combo, shape = "ellipse"), quantities = TRUE)
dev.off()

pdf(file = "./eulerVenn/prettyVenn.pdf")
plot(euler(combo, shape = "ellipse"), quantities = TRUE, fills = c("#036603", "#0070c6", "#942AAD", "#017D6F", "#5252b0"))
dev.off()
EOF2

source activate base
conda info --envs
source activate /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/envs/ChIPseeker

Rscript ${folder}/PBS/fiveWayVenn.R
EOF

sbatch ${folder}/PBS/fiveWayVenn.pbs