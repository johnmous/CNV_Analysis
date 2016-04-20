#!/bin/bash
# Title: A script to calculate the avarage coverage per gene per bam file, using htseq-count. Meant for CNV analysis
# Author: Ioannis Moustakas, i.moustakas@uva.nl

# sufix of bam files
sufix=sorted.bam

# list files in current directory
bamFiles=$(ls *$sufix)

for filename in ${bamFiles[@]};
do
  prefix=$(echo $filename | tr "$sufix" "\n")
  prefix=$prefix
  samtools view  $filename |cut -f 1-11 \
  |htseq-count -r pos -a 0 -s no -m intersection-nonempty -t gene --idattr locus_tag - Wild_Type_fixName.gff > $prefix"_ReadCountPerGene.txt"
done