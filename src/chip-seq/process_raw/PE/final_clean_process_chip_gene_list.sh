#!/usr/bin/bash

# This script is to process the _gene_dist.bed files from `macs_peak_call.sh` 
# to have only unique genes as one column.
# files are saved in #SCRATCH/riceTF/chip-seq/process_raw/final_gene_list
# Run `mv $(find . -name "gene_list_.txt") final_gene_list/` to move 
# final gene list to another folder. 

cd /scratch/jh6577/riceTF/chip-seq/process_raw

for file in $(find . -name '*summits_gene_dist.bed')
do
	cut -f 4 ${file} |sort|uniq > ${file%dist.bed}list.txt
	NUM_TARGET=$(wc -l ${file%dist.bed}list.txt|cut -d " " -f 1)
	echo $NUM_TARGET
	sed '1 i target' ${file%dist.bed}list.txt > ${file%dist.bed}list_$NUM_TARGET.txt
	rm ${file%dist.bed}list.txt
	#echo ${file%dist.bed}list.txt created
done
