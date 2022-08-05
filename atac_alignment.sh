#!/bin/bash

# Alignment pipeline for ATAC data. Assumes that all fastq files are in the same folder.
# Modified from existing pipeline 7/5/2022
# Assumes that all fastq files are contained within the same folder, and have suffix _1.fq.gz.

# Requires the "vanilla" environment, and the "qc" environment (for fastqc and multiqc)
# Version Doc: https://docs.google.com/document/d/1Jz95yFucvhhb8FvUM8XM37M4In3sG5mHYRaFzZuZiK8/edit

# TODO: Add comments for the generated script itself
# TODO: Test qc triggers correctly

# Taken from Ruoyun Wang

suffix1=_R1_001.fastq.gz
suffix2=_R2_001.fastq.gz

folder=$(cd "$(dirname "$0")";pwd)  # Save current folder as a variable

# Set up - Create required folder
mkdir -p PBS
mkdir -p deduplicated
mkdir -p trimmed
mkdir -p aligned
mkdir -p peaks_called
mkdir -p bigwig
mkdir -p log
mkdir -p heatmap

# CHANGE FOR FILE EXTENSION
count=$(find ./ -mindepth 1 -type f -name "*${suffix1}" -printf x | wc -c)  # Finds total number of files matching extension. Needed to know when to start qc. 
echo $count files found!

# Meta file to know when qc will begin (empty file)
cat >${folder}/'meta.txt' <<EOF
EOF

cd fastq
# Loop over all files, create the alignment script for each
for file in *${suffix1}; do
	base=$(basename "$file" ${suffix1})
	smallBase=${base%_S*}
	echo ${smallBase}
	cat >${folder}/PBS/$smallBase'.pbs' <<EOF
#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${smallBase}_%j.txt -e ${folder}/log/${smallBase}_%j.err.txt
cd ${folder}
################################
# Enter your code to run below #
################################

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate alignment

clumpify.sh \
	in1=./fastq/${base}${suffix1} \
	in2=./fastq/${base}${suffix2} \
	out1=deduplicated/${smallBase}_dedup${suffix1} \
	out2=deduplicated/${smallBase}_dedup${suffix2} \
	dedupe subs=2

gunzip -c deduplicated/${smallBase}_dedup${suffix1} > deduplicated/${smallBase}_dedup_R1.fastq
gunzip -c deduplicated/${smallBase}_dedup${suffix2} > deduplicated/${smallBase}_dedup_R2.fastq

bbduk.sh \
	in1=deduplicated/${smallBase}_dedup_R1.fastq \
	in2=deduplicated/${smallBase}_dedup_R2.fastq \
	out1=trimmed/${smallBase}_R1_trimmed.fastq \
	out2=trimmed/${smallBase}_R2_trimmed.fastq \
	ref=/dartfs-hpc/rc/lab/W/WangX/Nicholas/bbmap/resources/adapters.fa \
	ktrim=r k=21 mink=11 hdist=1 tpe tbo

rm deduplicated/${smallBase}_dedup*

module load bowtie/2.2.7

bowtie2 -p 10 \
	-x /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/GRCh38/grch38_1kgmaj \
	-1 trimmed/${smallBase}_R1_trimmed.fastq -2 trimmed/${smallBase}_R2_trimmed.fastq \
	--very-sensitive -k 10 > aligned/${smallBase}.sam

rm trimmed/${smallBase}_R1_trimmed.fastq trimmed/${smallBase}_R2_trimmed.fastq

samtools view -Sbo aligned/${smallBase}.bam aligned/${smallBase}.sam

samtools view -h aligned/${smallBase}.bam | awk '{if(\$3 != "chrM" && \$3 != "chrUn"){print \$0}}' | samtools view -Shb - > aligned/${smallBase}.filter.bam
samtools sort aligned/${smallBase}.filter.bam -o aligned/${smallBase}.filter.sorted.bam
samtools index aligned/${smallBase}.filter.sorted.bam

alignmentSieve --ATACshift --bam aligned/${smallBase}.filter.sorted.bam -o aligned/${smallBase}.shifted.bam
samtools sort ${folder}/aligned/${smallBase}.shifted.bam -o ${folder}/aligned/${smallBase}.shifted.nsorted.bam
samtools index ${folder}/aligned/${smallBase}.shifted.nsorted.bam

samtools sort -n ${folder}/aligned/${smallBase}.shifted.nsorted.bam -o ${folder}/aligned/${smallBase}.sortQuery.bam
samtools index ${folder}/aligned/${smallBase}.sortQuery.bam

rm aligned/${smallBase}.filter.bam
rm aligned/${smallBase}.bam 
rm aligned/${smallBase}.sam
rm aligned/${smallBase}.shifted.bam

bamCoverage --bam aligned/${smallBase}.shifted.nsorted.bam -o ${folder}/bigwig/${smallBase}_normalized.bw \
	--binSize 1 --normalizeUsing RPKM --ignoreDuplicates --extendReads -of bigwig -p max \
	-bl /dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/GRCh38/hg38-blacklist.v2.bed \
	--ignoreForNormalization chrX chrM chrRandom chrUn

Genrich -t ${folder}/aligned/${smallBase}.sortQuery.bam \
	-o ${folder}/peaks_called/${smallBase}.bed \
	-j  -y  -r  -e chrM  -v -E \
	/dartfs-hpc/rc/lab/W/WangX/Genomes_and_extra/GRCh38/hg38-blacklist.v2.bed

rm aligned/${smallBase}.sortQuery*

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R peaks_called/${smallBase}.bed \
	-S bigwig/${smallBase}_normalized.bw  \
	-o $folder/heatmap/${smallBase}_center.gz --missingDataAsZero -p max

plotHeatmap -m $folder/heatmap/${smallBase}_center.gz \
	-out $folder/heatmap/${smallBase}_center.pdf --colorList 'white,darkred'

rm $folder/heatmap/${smallBase}_center.gz

echo "${smallBase} completed!" >> ${folder}/'meta.txt'

# Checks to see if all the files have been analyzes
# if so, will run fastqc and multiqc on fastq and bam files 
# (requires access to Nick's pipeline folder)

currLine=\$(wc -l < ${folder}/meta.txt)
if ((\$currLine == $count)); then
	source activate base
    rmdir deduplicated/
    rmdir trimmed/
    cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/qc.sh ${folder}
    sh qc.sh
	cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/homerMotif.sh ${folder}
	sh homerMotif.sh peaks_called
	cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/ChIPseeker.sh ${folder}
	sh ChIPseeker.sh peaks_called
    rm ${folder}/meta.txt
	cp /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/atac_qc.sh ${folder}
	sh atac_qc.sh
fi
EOF
	cd ${folder}/log
	sbatch ../PBS/${smallBase}.pbs
	cd ${folder}/fastq
done
