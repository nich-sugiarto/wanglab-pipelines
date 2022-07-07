#!/bin/bash

# Written by Nick Sugiarto

# Simple script to move all files ending with .gz in fastq subfolders into fastq
# Assumes that fq files are compressed with .gz and contained in folder fastq/

# Trivial - No version control

cd fastq

# Loop over all folders in fastq
for folder in */; do
	mv ${folder}/*.gz ./
	rm -r ${folder}
done