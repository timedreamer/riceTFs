# The Snakemake pipeline for ChIP-Seq analysis

Author: Ji Huang

Date: 2019-09-03
Update: 2020-05-17

This repo is the snakemake analysis pipeline for the publicly available ChIP-Seq dataset (on NCBI-SRA). It contains most of the common ChIP-Seq analysis steps: 

![pipe](https://i.imgur.com/i8An0Rm.png)

Make sure `sample.csv` contains two columns: `treatment` and `control`. Each row is a comparison for peak calling. You can use same control for multiple peak calls. 

### Update: 2020-05-17
For the riceTFs projects, I stopped at the *q10_filter* step, and call `macs` manually. The script is the `macs_peak_call.sh`.


Files included:

1. `chip_pipeline_bt2.snakefile`: The main snakemake files contains all the rules.
2. `cluster_bt2.json`: Json file for the cluter configuration. The number of `thread` is defined in the `chip_pipeline_bt2.snakefile` file.
3. `config_bt2.yaml`: Contains most of the parameters that need to change.
4. `samples.csv`: CSV file for the sample you want to process.
5. `bt2_snake_slurm`: script to send out individual jobs.

To run the pipeline, modify files accordingly and run `./bt2_snake_slurm` in a computing node.

---
TO-DO:

1. Combine multiple SRA/Fastq files into one.
2. Better ways to handle SE and PE end files?


