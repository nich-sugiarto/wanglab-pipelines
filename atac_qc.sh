#!/bin/bash

# BROKEN - DO NOT USE

# Quality control specific for ATAC data. Includes library size estimation, fragment size distribution, 
# Promoter/Transcript (PT) body score, Nucleosome Free Regions (NFR) scores, and Transcription Start Site (TSS) scores.

# Code taken from https://bioconductor.org/packages/release/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html#Split_reads

# Version control available at https://docs.google.com/document/d/1lCYRKcYol6praX75HFxb90SaNrrugmx8aanxGQu3Ups/edit

# Script automatically generated using the "slouch" command on Wed Jun 29 15:29:27 EDT 2022

#TODO: Split reads? 
#TODO: Add multicorr plot

folder=$(cd "$(dirname "$0")";pwd)

mkdir -p ATACqc
mkdir -p PBS
mkdir -p log

cd aligned

for file in *.filter.sorted.bam; do
	base=$(basename "$file" ".filter.sorted.bam")
	echo ${base}
	mkdir -p ../ATACqc/${base}
	cat >./${base}_qc.R <<EOF
library(ATACseqQC)
library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # For PT score
library(Rsamtools)

pdf(file="${base}_libComplexity.pdf")
estimateLibComplexity(readsDupFreq("${file}"))
dev.off()

pdf(file="${base}_fragSize.pdf")
fragSize <- fragSizeDist("${file}", "${base}")
dev.off()

seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)

## bamfile tags to be read in
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))

print("Here!")

bamTop100 <- scanBam(BamFile("${file}", yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]]\$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
tags

print("Here!")

seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
which <- as(seqinformation, "GRanges")
bam <- readBamFile("${file}", tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
bam
bam <- shiftGAlignmentsList(bam)

txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
pt <- PTscore(bam, txs)

pdf(file="${base}_PT.pdf")
plot(pt\$log2meanCoverage, pt\$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
dev.off()

nfr <- NFRscore(as(bam, "GAlignments"), txs)
pdf(file="${base}_NFR.pdf")
plot(nfr\$log2meanCoverage, nfr\$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))
dev.off()

print("Here!")
tsse <- TSSEscore(bam, txs)
print("Here!")

pdf(file="${base}_NFR.pdf")
plot(100*(-9:10-.5), tsse\$values, type="b", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")
dev.off()
EOF

	cat >../PBS/${base}'_ATACqc.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${base}_ATACQC # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=16GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=5:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
################################
# Enter your code to run below #
################################
cd ${folder}/aligned

module load R/4.1.2

Rscript ${base}_qc.R

rm ${base}_qc.R
mv ${base}*.pdf ${folder}/ATACqc/${base}
EOF

	cd ../log
	sbatch ../PBS/${base}_ATACqc'.pbs'
	cd ${folder}/aligned
done