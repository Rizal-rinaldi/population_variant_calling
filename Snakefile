configfile: "config.yaml"
import os
from snakemake.utils import makedirs

#############################################
# Integrated Population Mapping + Variant Calling
# Author : Rizal (based on Carolina Pita Barros, WUR)
# Mode   : Local real-run
# Date   : 2025
#############################################

# --- Setup ---
if "OUTDIR" in config:
    print("\n[INFO] Saving to", config["OUTDIR"], "\n")
    workdir: config["OUTDIR"]

makedirs("logs_slurm")

ASSEMBLY   = config["ASSEMBLY"]
GFF_FILE   = config["GFF_FILE"]
PREFIX     = config["PREFIX"]
NUM_CHRS   = config["NUM_CHRS"]
PATHS      = config["PATHS_WITH_FILES"]

# --- Step 1: Identify FASTQ files ---
def define_to_map_together(given_path):
    TO_MAP_TOGETHER, FULLPATHS, SAMPLES = [], [], []
    possible_ext = [".fastq", ".fq.gz", ".fastq.gz", ".fq"]
    for root, dirs, files in os.walk(given_path):
        for name in files:
            for ext in possible_ext:
                if name.endswith(ext):
                    fullpath = os.path.join(root, name)
                    sample = fullpath.rsplit("/", 2)[1]
                    name_sub = name.rsplit("_", 2)[0]
                    if name_sub not in TO_MAP_TOGETHER:
                        TO_MAP_TOGETHER.append(name_sub)
                    FULLPATHS.append(fullpath)
                    if sample not in SAMPLES:
                        SAMPLES.append(sample)
    return TO_MAP_TOGETHER, FULLPATHS, SAMPLES

MAP, ALL_SAMPLES = {}, []
for file in PATHS.values():
    maptogether, fullpaths, samples = define_to_map_together(file)
    for var in maptogether:
        MAP[var] = [f for f in fullpaths if var in f]
    ALL_SAMPLES += samples


# --- Step 2: Workflow summary ---
rule all:
    input:
        # QC & Mapping
        expand("mapping_stats/qualimap/{sample_merged}/genome_results.txt", sample_merged=ALL_SAMPLES),
        "mapping_stats/sample_quality_summary.tsv",
	f"results/variant_calling/{PREFIX}.vcf.gz"
        # Variant calling & annotation
        #f"results/final_VCF/{PREFIX}.gatk_annot.stats",  
        #f"results/PCA/{PREFIX}_PCA.pdf",
	#f"results/hwe/{PREFIX}.hwe.tsv"


# ============================================================
# =============== MAPPING STAGE (BWA + SAMTOOLS) =============
# ============================================================
rule samtools_faidx:
    input:
        ASSEMBLY
    output:
        f"{ASSEMBLY}.fai"
    message:
        "Indexing reference FASTA with samtools faidx"
    shell:
        "samtools faidx {input}"
rule bwa_index:
    input:
        ref = ASSEMBLY,
        fai = rules.samtools_faidx.output
    output:
        multiext(ASSEMBLY, ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
    message:
        "Building BWA index"
    shell:
        "bwa-mem2 index {input.ref}"


def create_input_names(wildcards):
    return MAP[wildcards.sample]
# ============================================================
# =============== FASTP TRIMMING STAGE (NEW) =================
# ============================================================

rule fastp_trim:
    input:
        reads = create_input_names
    output:
        r1 = temp("trimmed_reads/{sample}_1.trimmed.fq.gz"),
        r2 = temp("trimmed_reads/{sample}_2.trimmed.fq.gz"),
        html = "trimmed_reads/{sample}_fastp.html",
        json = "trimmed_reads/{sample}_fastp.json"
    threads: 12
    message:
        "Running fastp on {wildcards.sample}"
    shell:
        """
        mkdir -p trimmed_reads
        fastp \
            -i {input.reads[0]} \
            -I {input.reads[1]} \
            -o {output.r1} \
            -O {output.r2} \
            -w {threads} \
            --html {output.html} \
            --json {output.json} \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20 \
            --length_required 50 \
            --detect_adapter_for_pe
        """


rule bwa_map:
    input:
        assembly = ASSEMBLY,
        idx = rules.bwa_index.output,
        r1 = rules.fastp_trim.output.r1,
        r2 = rules.fastp_trim.output.r2
    output:
        temp("mapped_reads/{sample}.bam")
    threads: 12
    message:
        "Mapping {wildcards.sample} with bwa-mem2"
    shell:
        r"""
        mkdir -p mapped_reads
        sample="{wildcards.sample}"
        echo "[Mapping] $sample"

        bwa-mem2 mem -t {threads} \
            -R "@RG\\tID:${{sample}}\\tSM:${{sample}}\\tPL:ILLUMINA" \
            {input.assembly} {input.r1} {input.r2} 2>/dev/null \
        | samblaster -r \
        | samtools view -b -o {output}
        """


def create_names_to_merge(wildcards):
    input_files = []
    mapped_dir = expand("mapped_reads/{sample}.bam", sample=MAP.keys())
    for file in mapped_dir:
        file_nodir = os.path.basename(file)
        if file_nodir.startswith(wildcards.sample_merged):
            input_files.append(file)
    return input_files

rule merge_mapped:
    input:
        create_names_to_merge
    output:
        temp("merged_reads/{sample_merged}.bam")
    run:
        shell("mkdir -p merged_reads")
        if len(input) > 1:
            shell("samtools merge -@ 12 {output} {input}")
        else:
            shell("mv {input} {output}")
 	
rule samtools_sort:
    input:
        "merged_reads/{sample_merged}.bam"
    output:
        "processed_reads/{sample_merged}.sorted.bam"
    shell:
        """
        mkdir -p processed_reads
        samtools sort -m 2G -@ 12 -O bam {input} -o {output}
        """

rule samtools_index:
    input:
        "processed_reads/{sample_merged}.sorted.bam"
    output:
        "processed_reads/{sample_merged}.sorted.bam.bai"
    shell:
        "samtools index -@ 12 {input}"

rule qualimap_report:
    input:
        bam = "processed_reads/{sample_merged}.sorted.bam",
        bai = "processed_reads/{sample_merged}.sorted.bam.bai"
    output:
        "mapping_stats/qualimap/{sample_merged}/genome_results.txt"
    params:
        outdir = "mapping_stats/qualimap/{sample_merged}/"
    shell:
        """
        mkdir -p mapping_stats/qualimap/{wildcards.sample_merged}
        unset DISPLAY
        qualimap bamqc -bam {input.bam} --java-mem-size=14G -nt 12 -outformat PDF -outdir {params.outdir}
        """

rule qualimap_summary:
    input:
        expand("mapping_stats/qualimap/{sample}/genome_results.txt", sample=ALL_SAMPLES)
    output:
        "mapping_stats/sample_quality_summary.tsv"
    params:
        script = os.path.join(workflow.basedir, "scripts/create_qualimap_summary.sh")
    shell:
        r"""
        set -euo pipefail
        bash {params.script}
        """

# ============================================================
# ============== VARIANT CALLING + ANNOTATION ================
# ============================================================

rule create_bam_list:
    input:
        bam = expand("processed_reads/{sample_merged}.sorted.bam", sample_merged=ALL_SAMPLES),
	bai = expand("processed_reads/{sample_merged}.sorted.bam.bai", sample_merged=ALL_SAMPLES)
    output:
        "bam_list.txt"
    shell:
        """
        realpath processed_reads/*.sorted.bam > {output}
        """

rule var_calling_freebayes:
    input:
        ref = ASSEMBLY,
        bam = rules.create_bam_list.output
    output:	
        vcf = "results/variant_calling/{prefix}.vcf.gz",
        idx = "results/variant_calling/{prefix}.vcf.gz.tbi"
    params:
        chunksize = 100000,
        scripts_dir = os.path.join(workflow.basedir, "scripts")
    shell:
        """
        mkdir -p results/variant_calling
        {params.scripts_dir}/freebayes-parallel.sh <({params.scripts_dir}/fasta_generate_regions.py {input.ref}.fai {params.chunksize}) 24 \
          -f {input.ref} -L {input.bam} \
          --use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 \
          --haplotype-length 0 --ploidy 2 --min-alternate-count 2 \
        | vcffilter -f 'QUAL > 20' | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
