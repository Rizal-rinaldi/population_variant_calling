Population Mapping and Variant Calling Version 2
This repository contains a Snakemake pipeline for integrated population mapping and variant calling. It is designed to handle raw FASTQ files, perform quality control, align reads to a reference genome, and call variants using a parallelized FreeBayes approach.

Author: Rizal (based on the framework by Carolina Pita Barros and WUR pipelines)

Mode: Local & Cluster execution

Date: 2025

Table of contents
Introduction

Workflow

Setup and Configuration

Execution

Output Structure

Introduction
This pipeline automates the transition from raw sequencing reads to a filtered population VCF file. It is optimized for multi-sample datasets and includes automated directory setup and logging.

Workflow
The pipeline implements the following steps:

Trimming: fastp for adapter removal and quality filtering (sliding window).

Indexing: Automated indexing of the reference genome with samtools and bwa-mem2.

Mapping: bwa-mem2 for alignment, piped through samblaster for duplicate marking.

Processing: samtools for merging library lanes, coordinate sorting, and indexing.

Quality Control: Qualimap BamQC for mapping statistics and coverage reports.

Variant Calling: Parallelized FreeBayes across genomic regions, followed by vcffilter (QUAL > 20).

Setup and Configuration
1. Prerequisites
Ensure you have the following tools installed (e.g., via Conda):

Snakemake

bwa-mem2

samtools & samblaster

fastp

qualimap

freebayes

vcflib & htslib

2. Configuration
Edit the config.yaml file to match your environment. This file decouples your data paths from the main logic.

YAML
OUTDIR: "results_bosjava_zoo/"
ASSEMBLY: "/path/to/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa"
GFF_FILE: "/path/to/Bos_taurus.ARS-UCD2.0.115.gff3"
PREFIX: "Var_calling_bosjava"
NUM_CHRS: 30
PATHS_WITH_FILES:
  group1: "/path/to/raw_fastq_data"
Execution
Dry run
Always check the workflow logic before execution:

Bash
snakemake -n
Local Run
Run the pipeline on a local machine using a specific number of cores:

Bash
snakemake --cores 24
SLURM Run
To run on an HPC cluster using SLURM, you can use the following command:

Bash
snakemake --cores 1 \
--jobs 20 \
--cluster "sbatch -p {partition} -c {threads} --mem=32G -o logs_slurm/%x_%j.out -e logs_slurm/%x_%j.err"
Output Structure
The pipeline organizes results into the following directories within your OUTDIR:

trimmed_reads/: Cleaned FASTQ files and fastp HTML/JSON reports.

processed_reads/: Final sorted and indexed BAM files (*.sorted.bam).

mapping_stats/: Qualimap reports and a consolidated sample_quality_summary.tsv.

results/variant_calling/: The final population VCF ({PREFIX}.vcf.gz) and its index.

logs_slurm/: Standard output and error logs from cluster jobs.