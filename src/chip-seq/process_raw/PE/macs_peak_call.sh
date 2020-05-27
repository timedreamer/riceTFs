#!/usr/bin/bash

## Author: Ji Huang
## Date: 2020-05-16

## Loops are unnecessary... I should probably use function...

module load macs2/intel/2.1.1
module load r/gnu/3.5.1
module load bedtools/intel/2.27.1

GENE_BED="/scratch/cgsb/coruzzi/jh6577/GenomeFile/rice_MSU7/all_gene.bed"

cd ../data/04aligned_bamPE/

# IDS
macs2 callpeak -t {ERR2286135,ERR2286136,ERR2286137}.sorted.rmdup.q10.bam -c ERR2286138.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir IDS1_keep  -n IDS1;

for file in IDS1_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r IDS1_keep ../../result/;rm -r IDS1_keep


# NAC2
macs2 callpeak -t SRR6336924.sorted.rmdup.q10.bam -c SRR6336925.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir NAC2_keep -n NAC2;

for file in NAC2_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r NAC2_keep ../../result/;rm -r NAC2_keep


##### TGAP #############
# TGAP untreated rep1
macs2 callpeak -t DRR015049.sorted.rmdup.q10.bam -c DRR015048.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir TGAP_untreated_rep1_keep -n TGAP_untreated_rep1;

for file in TGAP_untreated_rep1_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r TGAP_untreated_rep1_keep ../../result/;rm -r TGAP_untreated_rep1_keep

# TGAP untreated rep2
macs2 callpeak -t DRR015053.sorted.rmdup.q10.bam -c DRR015052.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir TGAP_untreated_rep2_keep -n TGAP_untreated_rep2;

for file in TGAP_untreated_rep2_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r TGAP_untreated_rep2_keep ../../result/;rm -r TGAP_untreated_rep2_keep

# TGAP treated rep1
macs2 callpeak -t DRR015051.sorted.rmdup.q10.bam -c DRR015050.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir TGAP_treated_rep1_keep -n TGAP_treated_rep1;

for file in TGAP_treated_rep1_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r TGAP_treated_rep1_keep ../../result/;rm -r TGAP_treated_rep1_keep

# TGAP treated rep2
macs2 callpeak -t DRR015055.sorted.rmdup.q10.bam -c DRR015054.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir TGAP_treated_rep2_keep -n TGAP_treated_rep2;

for file in TGAP_treated_rep2_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r TGAP_treated_rep2_keep ../../result/;rm -r TGAP_treated_rep2_keep


##### NAC127 #######
# NAC127H
macs2 callpeak -t {SRR10423448,SRR10423449,SRR10423450}.sorted.rmdup.q10.bam -c SRR10423460.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir NAC127H_keep -n NAC127H;

for file in NAC127H_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r NAC127H_keep ../../result/;rm -r NAC127H_keep

# NAC127N
macs2 callpeak -t {SRR10423454,SRR10423455,SRR10423456}.sorted.rmdup.q10.bam -c SRR10423462.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir NAC127N_keep -n NAC127N;

for file in NAC127N_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r NAC127N_keep ../../result/;rm -r NAC127N_keep

##### NAC129 #####
# NAC129H
macs2 callpeak -t {SRR10423451,SRR10423452,SRR10423453}.sorted.rmdup.q10.bam -c SRR10423461.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir NAC129H_keep -n NAC129H;

for file in NAC129H_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r NAC129H_keep ../../result/;rm -r NAC129H_keep

# NAC129N
macs2 callpeak -t {SRR10423457,SRR10423458,SRR10423459}.sorted.rmdup.q10.bam -c SRR10423463.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir NAC129N_keep -n NAC129N;

for file in NAC129N_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r NAC129N_keep ../../result/;rm -r NAC129N_keep


##### SNAC1 #########

#SNAC1-OE Normal
macs2 callpeak -t {SRR8746736,SRR8746737,SRR8746738}.sorted.rmdup.q10.bam -c {SRR8746739,SRR8746740,SRR8746741}.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir SNAC1-OE_keep -n SNAC1-OE;

for file in SNAC1-OE_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r SNAC1-OE_keep ../../result/;rm -r SNAC1-OE_keep

#ZH11 SNAC1 Normal
macs2 callpeak -t {SRR8746742,SRR8746743,SRR8746744}.sorted.rmdup.q10.bam -c {SRR8746745,SRR8746746,SRR8746747}.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir SNAC1-ZH11-Normal_keep -n SNAC1-ZH11-Normal;

for file in SNAC1-ZH11-Normal_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r SNAC1-ZH11-Normal_keep ../../result/;rm -r SNAC1-ZH11-Normal_keep

#ZH11 SNAC1 Drought
macs2 callpeak -t {SRR8746748,SRR8746749,SRR8746750}.sorted.rmdup.q10.bam -c {SRR8746751,SRR8746752,SRR8746753}.sorted.rmdup.q10.bam -f BAMPE -g 2E08 --keep-dup all --outdir SNAC1-ZH11-Drought_keep -n SNAC1-ZH11-Drought;

for file in SNAC1-ZH11-Drought_keep/*_summits.bed
do
	cut -f 1,2,3 ${file}|bedtools sort > ${file%.bed}.n.bed;
	bedtools closest -a ${file%.bed}.n.bed -b $GENE_BED -d |\
		awk -F "\t" 'BEGIN {OFS=FS} {if ($9<2001) {print $1,$2,$3,$7,$9}}' > ${file%.bed}_gene_dist.bed
	rm ${file%.bed}.n.bed
	echo "${file} finish"
done

cp -r SNAC1-ZH11-Drought_keep ../../result/;rm -r SNAC1-ZH11-Drought_keep