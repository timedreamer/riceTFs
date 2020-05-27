#!/usr/bin/bash

module load bedtools/intel/2.27.1

GENE_BED="/scratch/cgsb/coruzzi/jh6577/GenomeFile/rice_MSU7/all_gene.bed"

for file in ../GSE119575/*_summits.bed
do
	sed -e 's/^/Chr/' -e '/^ChrSyng_TIGR/d' ${file}| cut -f 1,2,3|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed

	echo "${file} finish"
	
done

