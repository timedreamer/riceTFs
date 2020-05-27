#!/usr/bin/bash

module load bedtools/intel/2.27.1

GENE_BED="/scratch/cgsb/coruzzi/jh6577/GenomeFile/rice_MSU7/all_gene.bed"

for file in ../GSE128495/*.txt
do
	sed '1d' ${file}|cut -f 2,3,4|bedtools sort > ${file%.txt}.bed
	bedtools closest -a ${file%.txt}.bed -b $GENE_BED -d |\
awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.txt}_gene_dist.bed
	
done






