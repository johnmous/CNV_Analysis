#!/bin/bash
# Title: A script to calculate the coverage per genomic position per bam file
# Author: Ioannis Moustakas, i.moustakas@uva.nl

# sufix of bam files
sufix=sorted.bam

# list files in current directory
bamFiles=$(ls *$sufix)

for filename in ${bamFiles[@]};
do
  prefix=$(echo $filename | tr "$sufix" "\n")
  prefix=$prefix
  genomeCoverageBed -ibam $filename -d |cut -f 2,3 > coverageAnalysis/$prefix"_Cover.txt"
done