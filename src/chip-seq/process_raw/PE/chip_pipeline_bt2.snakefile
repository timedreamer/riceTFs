# Read sample informations samples.csv
from csv import DictReader

TREATMENT,CONTROL=[],[]
with open("samples.csv") as f:
    for row in  DictReader(f):
        TREATMENT.append(row["treatment"])
        CONTROL.append(row["control"])

# Combine all SRA IDs.
SAMPLES = set(TREATMENT + CONTROL)

# Read in configuration files.
configfile: 'config_bt2.yaml'

SRADIR = config['SRADIR']

BT2INDEX = config['BT2INDEX']

REFGENOME = config['REFGENOME']

LENGTH_CUT = config['LENGTH_CUT']

MAPQ_FILTER = config['MAPQ_FILTER']

SEQ_END = config['SEQ_END']



# I need to define "SE" or "PE", to get the dump files. 
if config["SEQ_END"] == "SE":
   DUMP = "../data/02raw_fastq{SEQ_END}/{sample}.fastq.gz"
   FORMAT = "BAM"

elif config["SEQ_END"] == "PE":
   DUMP = ["../data/02raw_fastq{SEQ_END}/{sample}_1.fastq.gz",
   "../data/02raw_fastq{SEQ_END}/{sample}_2.fastq.gz"]
   FORMAT = "BAMPE"
else:
    sys.exit("SEQ_END only accepts PE or SE!")

# Rules -----------------------------------------------------------------------------------
localrules: all

rule all:
    input:
        # QC file
        expand("../data/03trim_fastq{SEQ_END}/{sample}_fastp.html", SEQ_END=SEQ_END, sample=SAMPLES),

        # # BT2 alignment
        #expand('../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.bam', SEQ_END=SEQ_END, sample=SAMPLES),
        # # BT2 MAPQ filter alignment
        expand('../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.q{MAPQ_FILTER}.bam', SEQ_END=SEQ_END, sample=SAMPLES, MAPQ_FILTER=MAPQ_FILTER),
 
        # # BT2 alignment index
        #expand('../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.bam.bai', SEQ_END=SEQ_END, sample=SAMPLES),
        expand('../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.q{MAPQ_FILTER}.bam.bai', SEQ_END=SEQ_END, sample=SAMPLES, MAPQ_FILTER=MAPQ_FILTER),

        # # BT2 alignment flagstat
        #expand('../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.bam.flagstat', SEQ_END=SEQ_END, sample=SAMPLES),
        #expand('../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.q{MAPQ_FILTER}.bam.flagstat', SEQ_END=SEQ_END, sample=SAMPLES, MAPQ_FILTER=MAPQ_FILTER),

        # Bigwig files
        #expand('../data/05bigwig{SEQ_END}/{sample}_q{MAPQ_FILTER}.bw', sample=SAMPLES, MAPQ_FILTER=MAPQ_FILTER, SEQ_END=SEQ_END),

        # Peak files
        #expand("../result/macs2_peaks/{treatment}_vs_{control}", zip, treatment=TREATMENT, control=CONTROL),

        # deepTools_matrix and three plots.
        #"../result/deepTools/computeMatrix.gz"



rule DOWNLOAD_SRA:
    output:
        temp("../data/01sra{SEQ_END}/{sample}.sra")
    threads: 2

    shell: """
        prefetch --output-directory "../data/01sra{SEQ_END}" {wildcards.sample}
    """

rule FASTQ_DUMP:
    input: 
        sra = "../data/01sra{SEQ_END}/{sample}.sra"
    output:
        temp(DUMP),

    threads: 4,

    run:
        if config['SEQ_END'] == "PE":
            shell(
                r"""
                /home/jh6577/miniconda3/envs/py37/bin/parallel-fastq-dump -s {input} --outdir "../data/02raw_fastq{SEQ_END}/" --gzip --threads {threads} --split-files
                """)
        else:
            shell(
                r"""
                /home/jh6577/miniconda3/envs/py37/bin/parallel-fastq-dump -s {input} --outdir "../data/02raw_fastq{SEQ_END}/" --gzip --threads {threads}
                """)

# PE and SE differs in "QC_ADAPTERS" and "BT2_ALIGN" two steps.
rule QC_ADAPTERS_PE:
    input:
        r1 = "../data/02raw_fastqPE/{sample}_1.fastq.gz",
        r2 = "../data/02raw_fastqPE/{sample}_2.fastq.gz",
    output:
        o1 = temp("../data/03trim_fastqPE/{sample}_1_trimed.fq.gz"),
        o2 = temp("../data/03trim_fastqPE/{sample}_2_trimed.fq.gz"),
        html = "../data/03trim_fastqPE/{sample}_fastp.html",
        json = "../data/03trim_fastqPE/{sample}_fastp.json"
    threads: 2

    shell:
       'fastp -q 20 -l {LENGTH_CUT} --thread {threads} -y \
       --adapter_sequence=AGATCGGAAGAGC --adapter_sequence_r2=AGATCGGAAGAGC \
       -t 1 -f 3 -i {input.r1} -I {input.r2} -o {output.o1} -O {output.o2} -h {output.html} -j {output.json}'

# BT2 alignment
rule BT2_ALIGN_PE:
    input:  
        seq1 = "../data/03trim_fastqPE/{sample}_1_trimed.fq.gz",
        seq2 = "../data/03trim_fastqPE/{sample}_2_trimed.fq.gz",
    output: temp("../data/04aligned_bamPE/{sample}.sorted.rmdup.bam"),
    threads: 12,
    log: 
        bt2 = "log/{sample}_bt2.txt",
        rmdup = "log/{sample}_rmdup.txt"
    message: "aligning {input}: {threads} threads",

    shell:"""
        bowtie2 -p 10 -1 {input.seq1} -2 {input.seq2} -x {BT2INDEX} 2>{log.bt2}|samblaster --removeDups 2>{log.rmdup}|samtools view -Sb -F 4 - |samtools sort -m 4G -@ 2 -T {output}.tmp -o {output}
        """



# rule QC_ADAPTERS_SE:
#     input:
#         r1 = "../data/02raw_fastqSE/{sample}.fastq.gz",
#     output:
#         o1 = temp("../data/03trim_fastqSE/{sample}_trimed.fq.gz"),
#         html = "../data/03trim_fastqSE/{sample}_fastp.html",
#         json = "../data/03trim_fastqSE/{sample}_fastp.json"
#     threads: 2

#     shell:
#        'fastp -q 20 -l {LENGTH_CUT} --thread {threads} -y -a AGATCGGAAGAGC \
#        -t 1 -f 3 -i {input.r1} -o {output.o1} -h {output.html} -j {output.json}'

## BT2 alignment
# rule BT2_ALIGN_SE:
#     input:  
#         o1 = "../data/03trim_fastq{SEQ_END}/{sample}_trimed.fq.gz",
#     output: "../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.bam",
#     threads: 12,
#     log: 
#         bt2 = "log/{sample}_bt2.txt",
#         rmdup = "log/{sample}_rmdup.txt"
#     message: "aligning {input}: {threads} threads",

#     shell:"""
#         bowtie2 -p 10 -U {input} -x {BT2INDEX} 2>{log.bt2}|samblaster --removeDups 2>{log.rmdup}|samtools view -Sb -F 4 - |samtools sort -m 4G -@ 2 -T {output}.tmp -o {output}
#         """

rule MAPQ_FILTER:
    input: "../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.bam",
    output: "../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.q{MAPQ_FILTER}.bam"
    threads: 1,
    message: "BAM filer on MAPQ = {MAPQ_FILTER}",

    shell:"""
        samtools view -b -q {MAPQ_FILTER} {input} > {output}
    """

## samtools index
rule BAM1_INDEX_FLAG:
    input:  
        "../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.bam"
    output: 
        "../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.bam.bai",
        "../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.bam.flagstat"
    threads: 1
    message: "index and flagstat bam {input}: {threads} threads"
    shell:"""
        samtools index {input}& samtools flagstat {input} > {output[1]};
        """

## samtools index
rule BAM2_INDEX_FLAG:
    input:  
        "../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.q{MAPQ_FILTER}.bam",
    output: 
        "../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.q{MAPQ_FILTER}.bam.bai",
        "../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.q{MAPQ_FILTER}.bam.flagstat"
    threads: 1
    message: "index and flagstat bam {input}: {threads} threads"
    shell:"""
        samtools index {input}& samtools flagstat {input} > {output[1]};
        """

## make_bigwigs
# rule MAKE_BIGWIGS:
#     input : 
#         "../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.q{MAPQ_FILTER}.bam",
#         "../data/04aligned_bam{SEQ_END}/{sample}.sorted.rmdup.q{MAPQ_FILTER}.bam.bai",
#     output: "../data/05bigwig{SEQ_END}/{sample}_q{MAPQ_FILTER}.bw"
#     threads: 4
#     message: "making bigwig for {input}"
#     shell:
#         """
#         bamCoverage -b {input[0]} --normalizeUsing RPKM --binSize 30 --smoothLength 120 -p {threads} --extendReads 200 -o {output}
#         """

# rule CALL_PEAKS:
#     input: 
#         control=control_bam, 
#         treatment=treatment_bam,

#     output: 
#         peak_dir = "../result/macs2_peaks/{treatment}_vs_{control}"
#     params:
#         name = "{treatment}_vs_{control}",

#     message: "call_peaks macs2"
#     shell:
#         """
#         module purge
#         module load macs2/intel/2.1.1
#         module load r/gnu/3.5.1
#         module load bedtools/intel/2.27.1
#         macs2 callpeak -t {input.treatment} \
#             -c {input.control} \
#              -f {FORMAT} -g {config[SIZE]} \
#             --outdir {output} -n {params.name}
#         cd {output} && Rscript {params.name}_model.r;
#         bedtools closest -a {GENE_BED} -b {params.name}_summits.bed -D ref |awk '{{ if (($11 < {RANGE}) && ($11 > -{RANGE})) {{ print }} }}'> {params.name}_intersected_genes.bed
#         """
#         # Have to use double braces to escape '{}' http://bit.ly/2ZHol0i Snakemake FAQ

# rule DEEP_TOOLS:
#     input:
#         expand("../data/05bigwig{SEQ_END}/{sample}_q{MAPQ_FILTER}.bw", SEQ_END=SEQ_END, sample=SAMPLES, MAPQ_FILTER=MAPQ_FILTER)
#     output:
#         "../result/deepTools/computeMatrix.gz"
#     threads: 8
#     message: "calculating deepTools_matrix for all bigwig files"
#     shell:
#         "computeMatrix scale-regions -S {input} -R {GENE_BED} -b 3000 -a 3000 --regionBodyLength 5000 -p 8 -o {output} &&\
#         plotProfile -m {output} -out ../result/deepTools/deepTools_profile.png --plotTitle 'profile' &&\
#         plotHeatmap -m {output} -out ../result/deepTools/deepTools_heatmap.png &&\
#         cd ../data/04align* &&\
#         plotFingerprint -b *q10.bam --skipZeros --region 1 --numberOfSamples 50000 -T 'Fingerprints' --plotFile ../../result/deepTools/deepTools_fingerprint.png"

# #|awk '{{ if (($11 < {RANGE}) && ($11 > -{RANGE})) {{ print }} }}'