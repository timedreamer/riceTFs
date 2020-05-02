# The snakemake pipeline for analyzing public SRA datasets from
# Pair-end (PE) RNA-Seq.

# Author: Ji Huang
# Date: 2018-08-08
# Last modified: 2020-04-30

# Read SRA IDs from samples.csv
with open('samples.csv', 'r') as t:
    SAMPLES = [line.rstrip('\n') for line in t]

SAMPLES=list(filter(None, SAMPLES))

configfile: 'config.yaml'

SRADIR = config['SRADIR']

STARREF = config['STARREF']

FCOUNTREF = config['FCOUNTREF']

FCRESULT_STAR = config['FCRESULT_STAR']
FCSUMMARY_STAR = config['FCSUMMARY_STAR']

TRIM_OUTDIR = config['TRIM_OUTDIR']


# Rules -----------------------------------------------------------------------------------
localrules: all

rule all:
    input:
        # fastp result
        expand('../data/trim_FASTQ/{sample}_fastp.html', sample=SAMPLES),

        FCRESULT_STAR,
        FCSUMMARY_STAR,
        

rule fastq_dump_SRA:
    output:
        temp("../data/sra/{sample}_1.fastq"),
        temp("../data/sra/{sample}_2.fastq"),
    threads: 4,
    resources: io_weigh = 10,
    priority: 10,

    shell:"""
    fasterq-dump --outdir {SRADIR} --split-files --threads {threads} {wildcards.sample}
"""

rule QCADAPTERS:
    input:
        r1 = "../data/sra/{sample}_1.fastq",
        r2 = "../data/sra/{sample}_2.fastq",
    output:
        o1 = temp("../data/trim_FASTQ/{sample}_1_trimed.fq.gz"),
        o2 = temp("../data/trim_FASTQ/{sample}_2_trimed.fq.gz"),
        html = "../data/trim_FASTQ/{sample}_fastp.html",
        json = "../data/trim_FASTQ/{sample}_fastp.json"
    threads: 2,
    resources: io_weigh = 5,
    priority: 50,

    shell:
       'fastp -q 20 -l {config[LENGTH_CUT]} --thread {threads} -y -a AGATCGGAAGAGC \
       --adapter_sequence_r2 AGATCGGAAGAGC -t 1 -i {input.r1} -I {input.r2} \
       -o {output.o1} -O {output.o2} -h {output.html} -j {output.json}'

rule STAR_ALN:
    input:
        r1 = "../data/trim_FASTQ/{sample}_1_trimed.fq.gz",
        r2 = "../data/trim_FASTQ/{sample}_2_trimed.fq.gz",

    output:
        bam = temp("../data/star_output/{sample}Aligned.out.bam"),
        sj = "../data/star_output/{sample}SJ.out.tab",
        log = "../data/star_output/{sample}Log.final.out",
    threads: 8

    shell:
        "STAR --genomeDir {STARREF} --readFilesCommand zcat \
        --runThreadN {threads} --readFilesIn {input.r1} {input.r2} \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within --twopassMode None --outFilterMultimapNmax 10 \
        --outReadsUnmapped None --outSAMstrandField intronMotif \
        --outFileNamePrefix '../data/star_output/{wildcards.sample}'"

rule FCSTAR:
    input:
        expand("../data/star_output/{sample}Aligned.out.bam",sample=SAMPLES)
    output:
        FCRESULT_STAR,
        FCSUMMARY_STAR,
    threads: 20,

    shell:
        "featureCounts -T {threads} -a {FCOUNTREF} -o {output[0]} {input} 2> {output[1]}"

