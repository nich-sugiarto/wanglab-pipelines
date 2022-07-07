#!/bin/bash

# Alignment pipeline for RNA data. Modified from existing pipeline 6/27/2022
# Assumes that all fastq files are contained in a fastq folder, and have suffix _1.fq.gz. See areas marked "CHANGE FOR FILE EXTENSION"
# Stolen and later adapted from Ruoyun Wang

# Requires the "vanilla" and "R" environments, and the "qc" environment (for fastqc and multiqc)

# Version Doc: https://docs.google.com/document/d/1LL0hB4eDcqZ7_UmBX_KwyiYlm5r1fBio-Y4hdBUdZQw/edit

# TODO: Ask/find out what the read length was
# TODO: Check to see if fastqc is recursive (Important!)

mkdir -p counts
mkdir -p clumped
mkdir -p trimmed
mkdir -p aligned
mkdir -p bigwig
mkdir -p results
mkdir -p PBS
mkdir -p log

folder=$(cd "$(dirname "$0")";pwd)

count=$(find ./ -mindepth 1 -type f -name "*1.fastq.gz" -printf x | wc -c)  # Finds total number of files matching extension. Needed to know when to start qc. CHANGE FOR FILE EXTENSION
echo There are $count sets of files

# Meta file to know when qc will begin (empty file)
cat >${folder}/'meta.txt' <<EOF
EOF
cd fastq
for file in *1.fastq.gz  # CHANGE FOR FILE EXTENSION
do

	base=$(basename "$file" "_R1.fastq.gz")  # CHANGE FOR FILE EXTENSION]
	smallBase=${base%_S1*}
	echo ${smallBase}

	cat >${folder}/counts/${smallBase}_RPKM'.R' <<EOF
library(edgeR)
leng <- read.table("${smallBase}_featurecounts_Length.txt",header = TRUE,skip=1)
data <- read.table("${smallBase}_featurecounts_Count.txt",header = TRUE,skip=1)
geneid <- read.table("${smallBase}_featurecounts_Name.txt",header = TRUE,skip=1)

names(data)[names(data) == "aligned.${smallBase}_sorted.bam"] <- "Counts"
rpkm <- rpkm(data,leng)
names(rpkm)[names(rpkm) == "Length"] <- "RPKM"
final<- data.frame(geneid,data,rpkm)
write.table(final,file="${smallBase}_RPKM.csv",row.names = FALSE,quote = FALSE,sep = ",")
EOF

	cat >${folder}/PBS/$smallBase'.sbatch' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${smallBase}_alignment # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=10

#Number of memory
#SBATCH --mem-per-cpu=32GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=10:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
################################
# Enter your code to run below #
################################
cd ${folder}
pwd

source activate vanilla
	
mkdir -p fastqc

# CHANGE FOR FILE EXTENSION
clumpify.sh \
	in1=${folder}/${base}_1.fq.gz \
	in2=${folder}/${base}_2.fq.gz \
	out1=clumped/${smallBase}_R1_clumped.fastq.gz \
	out2=clumped/${smallBase}_R2_clumped.fastq.gz 
	
gunzip clumped/${smallBase}*_clumped.fastq.gz

bbduk.sh \
	in1=clumped/${smallBase}_R1_clumped.fastq \
	in2=clumped/${smallBase}_R2_clumped.fastq \
	out1=trimmed/${smallBase}_R1_trimmed.fastq \
	out2=trimmed/${smallBase}_R2_trimmed.fastq \
	ref=/dartfs-hpc/rc/lab/W/WangX/Nicholas/bbmap/resources/adapters.fa \
	ktrim=r k=21 mink=11 hdist=1 tpe tbo

rm clumped/${smallBase}*

STAR \
	--genomeDir /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/hg38_STAR  \
	--runThreadN 8 \
	--readFilesIn trimmed/${smallBase}_R1_trimmed.fastq trimmed/${smallBase}_R2_trimmed.fastq \
	--outFileNamePrefix aligned/${smallBase} \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes Standard

rm -r aligned/${smallBase}_STARtmp
rm trimmed/${smallBase}*

samtools sort -n -o aligned/${smallBase}_sorted.bam aligned/${smallBase}Aligned.sortedByCoord.out.bam
samtools index aligned/${smallBase}Aligned.sortedByCoord.out.bam

bamCoverage -b aligned/${smallBase}Aligned.sortedByCoord.out.bam -o bigwig/${smallBase}.bw

featureCounts -T 8 -s 0 \
	-g gene_name -a /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/hg38.knownGene.gtf \
	-o counts/${smallBase}_featurecounts.txt \
	-p aligned/${smallBase}_sorted.bam

cp counts/${smallBase}_featurecounts.txt.summary results

cd counts/
	
cut -f 7 ${smallBase}_featurecounts.txt > ${smallBase}_featurecounts_Count.txt
cut -f 1 ${smallBase}_featurecounts.txt > ${smallBase}_featurecounts_Name.txt
cut -f 6 ${smallBase}_featurecounts.txt > ${smallBase}_featurecounts_Length.txt

source activate R
Rscript ${smallBase}_RPKM.R

echo "${smallBase} completed!" >> ${folder}/'meta.txt'

# Checks to see if all the files have been analyzes
# if so, will run fastqc and multiqc on fastq and bam files 
# (requires access to Nick's pipeline folder)

currLine=$(wc -l < ${folder}/meta.txt)
if ((\$currLine == $count));
then
    cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/sugiarto_qc.sh ${folder}
    sh sugiarto_qc_folder.sh
    rm ${folder}/meta.txt
	rmdir clumped/
	rmdir trimmed/
fi
EOF

	cd ${folder}/log

	sbatch ${folder}/PBS/$smallBase'.sbatch'
	cd ${folder}
done