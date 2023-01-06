#!/bin/bash

mkdir -p combinedCounts

cd finalCounts
for file in *_perfMatch.txt; do
    base=$(basename "${file}" "_perfMatch.txt")
    touch ../combinedCounts/${base}.txt
    cat ${file} > ../combinedCounts/${base}.txt
    cat ${base}_LevenshteinCounts.txt >> ../combinedCounts/${base}.txt
done