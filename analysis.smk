import pandas as pd
import glob
import os
import random 

# Include plotting rules from separate file
include: "plotting.smk"

file_path = "/home/AD/wdecoster/inquiSTR_paper/1000G_cohort.tsv"
# Load a TSV file into a DataFrame, containing URLs of 1000 Genomes cram files
df = pd.read_csv(file_path, sep='\t', usecols=['sample', 'hg38_path', 'source', 'Superpopulation code'])
df = df[df["source"] == "Noyvert/Schloissnig"] # these are cram files

# PCA cohort: 100 mixed samples, 20 per superpopulation
# we can't take a random sample, because that would make the result of the workflow change every iteration, triggering a rerun
samples_per_superpopulation = 20
selected_samples = pd.concat([
    df[df["Superpopulation code"] == "AFR"].head(samples_per_superpopulation),
    df[df["Superpopulation code"] == "AMR"].head(samples_per_superpopulation),
    df[df["Superpopulation code"] == "EAS"].head(samples_per_superpopulation),
    df[df["Superpopulation code"] == "EUR"].head(samples_per_superpopulation),
    df[df["Superpopulation code"] == "SAS"].head(samples_per_superpopulation),
]).drop_duplicates(subset="sample")

# Relatedness cohort: 100 EUR samples, must include HG01512 (father of HG01514) and realigned HG01514
# The first 20 EUR overlap with the PCA cohort and don't need to be genotyped twice
relatedness_samples = pd.concat([
    df[df["sample"] == "HG01512"],
    df[df["Superpopulation code"] == "EUR"],
]).drop_duplicates(subset="sample").head(100)
# HG01514 is added via the realign rule and included in the relatedness combine

# All unique samples needing genotyping (union of PCA and relatedness cohorts)
# HG01514 has special handling via the realign rule
all_genotype_samples = pd.concat([selected_samples, relatedness_samples]).drop_duplicates(subset="sample")

# Build lookup dicts for per-sample polymorphic genotyping
POLYMORPHIC_SAMPLE_URLS = {row['sample']: row['hg38_path'] for _, row in all_genotype_samples.iterrows()}
POLYMORPHIC_SAMPLES = list(POLYMORPHIC_SAMPLE_URLS.keys())

# Sample lists for separate combine steps
PCA_SAMPLES = list(selected_samples['sample'])
RELATEDNESS_SAMPLES = list(relatedness_samples['sample'])


reference = "/home/AD/wdecoster/database/GRCh38.fa"
inquiSTR = "/home/AD/wdecoster/repositories/inquiSTR/target/x86_64-unknown-linux-musl/release/inquiSTR"
TRGT = "/home/AD/wdecoster/bin/trgt",
LongTR = "/home/AD/wdecoster/anaconda3/envs/longtr/bin/LongTR"
STRAGLR = "/home/AD/wdecoster/repositories/straglr/straglr.py" # requires environment with straglr dependencies

# Benchmark parameters
TECHNOLOGIES = ["ont", "pacbio"]
THREAD_COUNTS = list(range(1, 13))  # 1 to 12 threads
REPLICATES = list(range(1,6))  # 5 replicates
MAX_LOCUS = 10000  # limit to loci shorter than 10kb for genotyping

# Randomize order of thread counts for benchmarking
THREAD_ORDER = THREAD_COUNTS * len(REPLICATES)
random.shuffle(THREAD_ORDER)

# if the "tool_versions.txt" file exists, remove it to ensure fresh capture of tool versions
if os.path.exists("tool_versions.txt"):
    os.remove("tool_versions.txt")

# pathogenic expansions to genotype are in /home/AD/wdecoster/inquiSTR_paper/puretarget-data and have to be globbed to a list and genotyped with inquiSTR --pathogenic
puretarget_files = [os.path.basename(f).replace(".bam", "") for f in glob.glob("/home/AD/wdecoster/inquiSTR_paper/puretarget-data/*.bam")]

rule all:
    input:
        "benchmarking/results.tsv",
        "benchmarking/runtime_plot.html",
        "benchmarking/memory_plot.html",
        #"genotyping/adotto_combined_selected_samples.tsv",
        # PacBio TRGT
        expand("tool_comparison/pacbio-trgt-adotto_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-trgt-adotto_rep{replicate}.time", replicate=REPLICATES),
        # PacBio inquiSTR
        expand("tool_comparison/pacbio-inquistr-adotto_rep{replicate}.inq.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        # PacBio inquiSTR filter catalog (TRGT format)
        expand("tool_comparison/adotto-variable-catalog-pacbio_rep{replicate}.bed", replicate=REPLICATES),
        expand("tool_comparison/filter-inquistr-pacbio_rep{replicate}.time", replicate=REPLICATES),
        # PacBio inquiSTR + filter catalog (LongTR format)
        expand("tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.inq.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.time", replicate=REPLICATES),
        expand("tool_comparison/adotto-variable-catalog-pacbio-longtr_rep{replicate}.bed", replicate=REPLICATES),
        expand("tool_comparison/filter-inquistr-pacbio-longtr_rep{replicate}.time", replicate=REPLICATES),
        # PacBio TRGT filtered
        expand("tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        # PacBio LongTR
        expand("tool_comparison/pacbio-longtr-adotto_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-longtr-adotto_rep{replicate}.time", replicate=REPLICATES),
        # PacBio LongTR filtered
        expand("tool_comparison/pacbio-longtr-adotto-filtered_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-longtr-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        # ONT inquiSTR
        expand("tool_comparison/ont-inquistr-adotto_rep{replicate}.inq.gz", replicate=REPLICATES),
        expand("tool_comparison/ont-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        # ONT inquiSTR filter catalog (TRGT format, unused by LongTR)
        expand("tool_comparison/adotto-variable-catalog-ont_rep{replicate}.bed", replicate=REPLICATES),
        expand("tool_comparison/filter-inquistr-ont_rep{replicate}.time", replicate=REPLICATES),
        # ONT inquiSTR + filter catalog (LongTR format)
        expand("tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.inq.gz", replicate=REPLICATES),
        expand("tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.time", replicate=REPLICATES),
        expand("tool_comparison/adotto-variable-catalog-ont-longtr_rep{replicate}.bed", replicate=REPLICATES),
        expand("tool_comparison/filter-inquistr-ont-longtr_rep{replicate}.time", replicate=REPLICATES),
        # ONT LongTR
        expand("tool_comparison/ont-longtr-adotto_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/ont-longtr-adotto_rep{replicate}.time", replicate=REPLICATES),
        # ONT LongTR filtered
        expand("tool_comparison/ont-longtr-adotto-filtered_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/ont-longtr-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        # inquiSTR convert (rep1 only - representative conversion)
        "tool_comparison/pacbio-trgt-adotto_rep1.inq.gz",
        "tool_comparison/pacbio-longtr-adotto_rep1.inq.gz",
        "tool_comparison/ont-longtr-adotto_rep1.inq.gz",
        # inquiSTR benchmark: compare inquiSTR genotypes against TRGT/LongTR genotypes
        "tool_comparison/benchmark_inquistr_vs_trgt_pacbio.tsv",
        "tool_comparison/benchmark_inquistr_vs_trgt_pacbio.html",
        "tool_comparison/benchmark_inquistr_vs_trgt_pacbio_discrepancies.tsv",
        "tool_comparison/benchmark_inquistr_vs_longtr_pacbio.tsv",
        "tool_comparison/benchmark_inquistr_vs_longtr_pacbio.html",
        "tool_comparison/benchmark_inquistr_vs_longtr_pacbio_discrepancies.tsv",
        "tool_comparison/benchmark_inquistr_vs_longtr_ont.tsv",
        "tool_comparison/benchmark_inquistr_vs_longtr_ont.html",
        "tool_comparison/benchmark_inquistr_vs_longtr_ont_discrepancies.tsv",
        # inquiSTR --require-spanning vs TRGT benchmark
        "tool_comparison/pacbio-inquistr-adotto-requirespanning.inq.gz",
        "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio.tsv",
        "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio.html",
        "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio_discrepancies.tsv",
        expand("benchmarking/accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        "tool_comparison/results.tsv",
        "tool_comparison/runtime_plot.html",
        "tool_comparison/memory_plot.html",
        "tool_comparison/straglr_chr21_catalog.bed",
        expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-inquistr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES),
        "tool_comparison/straglr_chr21_results.tsv",
        "tool_comparison/straglr_chr21_runtime_plot.html",
        "tool_comparison/straglr_chr21_memory_plot.html",
        expand("tool_comparison/straglr_chr21_accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        expand("tool_comparison/straglr_chr21_accuracy_{technology}.html", technology=TECHNOLOGIES),
        expand("tool_comparison/inquistr_chr21_accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        expand("tool_comparison/inquistr_chr21_accuracy_{technology}.html", technology=TECHNOLOGIES),
        "tool_versions.txt",
        "straglr_version.txt",
        "inquiSTR_version.txt"

rule straglr_chr21:
    """Target rule for the chr21 STRaglr vs inquiSTR runtime comparison subset."""
    input:
        "tool_comparison/straglr_chr21_catalog.bed",
        expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.tsv", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.bed", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-inquistr-chr21_rep{replicate}.inq.gz", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-inquistr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES),
        expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.inq.gz", technology=TECHNOLOGIES, replicate=["1"]),
        "tool_comparison/straglr_chr21_results.tsv",
        "tool_comparison/straglr_chr21_runtime_plot.html",
        "tool_comparison/straglr_chr21_memory_plot.html",
        expand("tool_comparison/straglr_chr21_accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        expand("tool_comparison/straglr_chr21_accuracy_{technology}.html", technology=TECHNOLOGIES),
        expand("tool_comparison/inquistr_chr21_accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        expand("tool_comparison/inquistr_chr21_accuracy_{technology}.html", technology=TECHNOLOGIES)

rule versions:
    input:
        "tool_versions.txt",
        "straglr_version.txt"

rule straglr_version:
    """Capture STRaglr version using the STRaglr conda environment."""
    input:
        STRAGLR
    output:
        "straglr_version.txt"
    params:
        straglr = STRAGLR
    conda:
        "envs/straglr.yml"
    log:
        "logs/straglr_version.log"
    shell:
        """
        {params.straglr} --version > {output} 2> {log} || {params.straglr} -h > {output} 2>> {log}
        """

rule inquiSTR_version:
    """Capture inquiSTR version to trigger reruns when version changes"""
    input:
        inquiSTR
    output:
        "inquiSTR_version.txt"
    params:
        inquiSTR = inquiSTR
    log:
        "logs/inquiSTR_version.log"
    shell:
        """
        {params.inquiSTR} --version > {output} 2>&1
        """

rule capture_tool_versions:
    input:
        straglr_version = "straglr_version.txt"
    output:
        "tool_versions.txt"
    params:
        inquiSTR = inquiSTR,
        TRGT = TRGT,
        LongTR = LongTR
    log:
        "logs/capture_tool_versions.log"
    shell:
        """
        {{
            echo "Tool Versions Summary"
            echo "====================="
            echo ""
            echo "inquiSTR:"
            {params.inquiSTR} --version 2>&1 || echo "Version command not available"
            echo ""
            echo "TRGT:"
            {params.TRGT} --version 2>&1 || echo "Version command not available"
            echo ""
            echo "LongTR:"
            {params.LongTR} --version 2>&1 || echo "Version command not available"
            echo ""
            echo "STRaglr:"
            cat {input.straglr_version}
            echo ""
            echo "Generated on: $(date)"
        }} > {output} 2> {log}
        """

# using polymorphic repeats from illumina https://zenodo.org/records/8329210/files/polymorphic_repeats.hg38.bed?download=1
rule polymorphic:
    input:
        "polymorphic/repeats_relate.tsv",
        "polymorphic/repeats_pca.html",
        "polymorphic/repeats_pca_colored.html",

rule create_manifest_selected_samples:
    output:
        manifest = "selected_samples_manifest.tsv"
    log:
        "logs/create_manifest.log"
    run:
        with open(output.manifest, 'w') as f:
            f.write("bam_path\tsample_name\n")
            for idx, row in selected_samples.iterrows():
                f.write(f"{row['hg38_path']}\t{row['sample']}\n")

rule genotype_polymorphic_sample:
    input:
        bed = "polymorphic_repeats.hg38.bed",
    output:
        "polymorphic/individual/{sample}.tsv"
    params:
        url = lambda wildcards: POLYMORPHIC_SAMPLE_URLS[wildcards.sample],
        reference = "/home/AD/wdecoster/database/1KG_ONT_VIENNA_hg38.fa",
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS,
        cram = lambda wildcards: f"polymorphic/tmp/{wildcards.sample}.cram"
    threads: 4
    log:
        "logs/genotype_polymorphic_{sample}.log"
    shell:
        """
        mkdir -p polymorphic/tmp 2> {log}
        wget -q -O {params.cram} {params.url} 2>> {log}
        wget -q -O {params.cram}.crai {params.url}.crai 2>> {log}
        {params.inquiSTR} call {params.cram} \
            --region-file {input.bed} \
            --reference {params.reference} \
            --threads {threads} \
            --unphased \
            --max-locus {params.max_locus} > {output} 2>> {log}
        rm -f {params.cram} {params.cram}.crai 2>> {log}
        """

rule realign_HG01514:
    input:
        "/home/AD/wdecoster/study322-ONT_genomes/all_files/rr_HG01514/LCYT/209418/v7.3.11/cram/rr_HG01514§LCYT.cram"
    output:
        cram = "polymorphic/tmp/HG01514_realigned.cram",
        crai = "polymorphic/tmp/HG01514_realigned.cram.crai"
    params:
        reference = "/home/AD/wdecoster/database/1KG_ONT_VIENNA_hg38.fa",
    threads:
        20
    conda:
        "envs/minimap2.yml" # also includes samtools
    log:
        "logs/realign_HG01514.log"
    shell:
        """
        mkdir -p polymorphic/tmp
        samtools fastq -@ {threads} {input} 2>> {log} \
        | minimap2 -ax map-ont -t {threads} {params.reference} - 2>> {log} \
        | samtools sort --write-index -o {output.cram} - 2>> {log}
        """

rule genotype_HG01514:
    # run inquiSTR call on the realigned HG01514
    input:
        bed = "polymorphic_repeats.hg38.bed",
        cram = "polymorphic/tmp/HG01514_realigned.cram",
        crai = "polymorphic/tmp/HG01514_realigned.cram.crai"
    output:
        "polymorphic/individual/HG01514.tsv"
    params:
        reference = "/home/AD/wdecoster/database/1KG_ONT_VIENNA_hg38.fa",
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS
    threads: 4
    log:
        "logs/genotype_HG01514.log"
    shell:
        """
        {params.inquiSTR} call {input.cram} \
            --region-file {input.bed} \
            --reference {params.reference} \
            --threads {threads} \
            --unphased \
            --max-locus {params.max_locus} > {output} 2> {log}
        """


rule combine_pca:
    input:
        expand("polymorphic/individual/{sample}.tsv", sample=PCA_SAMPLES),
    output:
        "polymorphic/repeats_pca_combined.tsv"
    log:
        "logs/combine_pca.log"
    threads:
        4
    params:
        inquiSTR = inquiSTR,
    shell:
        """
        {params.inquiSTR} combine --threads {threads} {input} > {output} 2> {log}
        """

rule combine_relatedness:
    input:
        samples = expand("polymorphic/individual/{sample}.tsv", sample=RELATEDNESS_SAMPLES),
        hg01514 = "polymorphic/individual/HG01514.tsv"
    output:
        "polymorphic/repeats_relatedness_combined.tsv"
    log:
        "logs/combine_relatedness.log"
    threads:
        4
    params:
        inquiSTR = inquiSTR,
    shell:
        """
        {params.inquiSTR} combine --threads {threads} {input.samples} {input.hg01514} > {output} 2> {log}
        """



rule polymorphic_relate:
    input:
        combined = "polymorphic/repeats_relatedness_combined.tsv",
    output:
        "polymorphic/repeats_relate.tsv"
    threads: 16
    params:
        inquiSTR = inquiSTR,
        min_spacing = 100000, # also the default, but set explicitly for transparency
        tolerance = 1, # also the default, but set explicitly for transparency
    log:
        "logs/polymorphic_relate.log"
    shell:
        """
        {params.inquiSTR} relate \
            --output {output} \
            --threads {threads} \
            --min-spacing {params.min_spacing} \
            --tolerance {params.tolerance} \
            {input.combined} \
            > {log} 2>&1
        """


rule polymorphic_pca:
    input:
        combined = "polymorphic/repeats_pca_combined.tsv",
    output:
        plot = "polymorphic/repeats_pca.html",
        scores = "polymorphic/repeats_pca_scores.tsv"
    threads: 16
    params:
        inquiSTR = inquiSTR
    log:
        "logs/polymorphic_pca.log"
    shell:
        """
        {params.inquiSTR} pca \
            --output {output.plot} \
            --threads {threads} \
            --scores {output.scores} \
            {input.combined} \
            > {log} 2>&1
        """


rule genotype_adotto:
    input:
        manifest = "selected_samples_manifest.tsv",
        bed = "adotto.bed.gz",
        version = "inquiSTR_version.txt"
    output:
        "genotyping/adotto_combined_selected_samples.tsv"
    params:
        reference = reference,
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS # limit to loci shorter than 10kb for genotyping
    threads: 16
    log:
        "logs/genotype_adotto.log"
    shell:
        """
        {params.inquiSTR} batch {input.manifest} \
            --region-file {input.bed} \
            --output {output} \
            --noextend \
            --threads {threads} \
            --parallel-samples 4 \
            --max-locus {params.max_locus} \
            --reference {params.reference} \
            --resume \
            --save-individual genotyping/adotto_individual \
            > {log} 2>&1
        """


# Benchmark rules
rule benchmark_call:
    input:
        cram = "{technology}.cram", # for ont, this refers to the `ont_downsampled.cram` file
        bed = "adotto.bed.gz",
        version = "inquiSTR_version.txt"
    output:
        result = "benchmarking/data/{technology}_threads{threads}_rep{replicate}.tsv",
        timing = "benchmarking/data/{technology}_threads{threads}_rep{replicate}.time"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS # limit to loci shorter than 10kb for genotyping
    threads: lambda wildcards: int(wildcards.threads)
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    log:
        "logs/benchmark_{technology}_threads{threads}_rep{replicate}.log"
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.cram} \
            --region-file {input.bed} \
            --threads {threads} \
            --noextend \
            --max-locus {params.max_locus} \
            --reference {params.reference} \
            > {output.result} 2> {log}
        """


rule aggregate_benchmark_results:
    input:
        expand("benchmarking/data/{technology}_threads{threads}_rep{replicate}.time",
               technology=TECHNOLOGIES,
               threads=THREAD_COUNTS,
               replicate=REPLICATES)
    output:
        "benchmarking/results.tsv"
    run:
        import re
        results = []
        
        for timing_file in input:
            # Parse filename to get metadata
            match = re.search(r'benchmarking/data/(\w+)_threads(\d+)_rep(\d+)\.time', timing_file)
            if match:
                technology, threads, replicate = match.groups()
                
                # Parse timing output from /usr/bin/time -v
                elapsed_time = None
                max_memory_kb = None
                with open(timing_file, 'r') as f:
                    for line in f:
                        if 'Elapsed (wall clock) time' in line:
                            # Format is "Elapsed (wall clock) time (h:mm:ss or m:ss): 6:15.40"
                            # Split at "): " to get the time value
                            time_str = line.split('): ', 1)[1].strip()
                            # Convert to seconds
                            parts = time_str.split(':')
                            if len(parts) == 3:  # h:mm:ss
                                h, m, s = parts
                                elapsed_time = int(h) * 3600 + int(m) * 60 + float(s)
                            elif len(parts) == 2:  # mm:ss
                                m, s = parts
                                elapsed_time = int(m) * 60 + float(s)
                            elif len(parts) == 1:  # just ss
                                elapsed_time = float(parts[0])
                        elif 'Maximum resident set size' in line:
                            # Format is "Maximum resident set size (kbytes): 123456"
                            max_memory_kb = int(line.split(':')[1].strip())
                
                if elapsed_time is not None and max_memory_kb is not None:
                    results.append({
                        'technology': technology,
                        'threads': int(threads),
                        'replicate': int(replicate),
                        'elapsed_seconds': elapsed_time,
                        'max_memory_gb': max_memory_kb / (1024 * 1024)  # Convert KB to GB
                    })
        
        # Write results
        import pandas as pd
        df = pd.DataFrame(results)
        df = df.sort_values(['technology', 'threads', 'replicate'])
        df.to_csv(output[0], sep='\t', index=False)


rule download_adotto:
    output:
        catalog = "adotto_TRGT.bed.gz",
    log:
        "logs/download_adotto.log"
    params:
        url = "https://zenodo.org/records/8329210/files/adotto_repeats.hg38.bed.gz",
    shell:
        """
        wget -O {output.catalog} {params.url} &> {log}
        """

rule decompress_catalog:
    input:
        catalog = "adotto_TRGT.bed.gz"
    output:
        catalog = "adotto_LongTR.bed"
    shell:
        "zcat {input.catalog} | awk 'BEGIN{{OFS=\"\\t\"}} {{match($4, /MOTIFS=([^;]+)/, m); $4=m[1]; print}}' > {output.catalog}"

rule TRGT_adotto:
    input:
        catalog = "adotto_TRGT.bed.gz",
        pacbio = "pacbio.cram",
    output:
        vcf = "tool_comparison/pacbio-trgt-adotto_rep{replicate}.vcf.gz",
        timing = "tool_comparison/pacbio-trgt-adotto_rep{replicate}.time"
    log:
        "logs/TRGT_adotto_rep{replicate}.log"
    params:
        TRGT = "/home/AD/wdecoster/bin/trgt",
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.TRGT} genotype --genome {params.reference} \
            --repeats {input.catalog} \
            --reads {input.pacbio} \
            --threads {threads} \
            --output-prefix tool_comparison/pacbio-trgt-adotto_rep{wildcards.replicate} &> {log}
        """

rule inquiSTR_adotto:
    input:
        catalog = "adotto_TRGT.bed.gz",
        pacbio = "pacbio.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/pacbio-inquistr-adotto_rep{replicate}.inq.gz",
        timing = "tool_comparison/pacbio-inquistr-adotto_rep{replicate}.time"
    log:
        "logs/inquiSTR_adotto_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS # limit to loci shorter than 10kb for genotyping
    threads:
        4
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.pacbio} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend 2> {log} | gzip > {output.inq} 2>> {log}
        """

rule filter_inquiSTR_pacbio:
    input:
        pacbio_inq = "tool_comparison/pacbio-inquistr-adotto_rep{replicate}.inq.gz",
    output:
        catalog = "tool_comparison/adotto-variable-catalog-pacbio_rep{replicate}.bed",
        timing = "tool_comparison/filter-inquistr-pacbio_rep{replicate}.time"
    log:
        "logs/filter_inquistr_pacbio_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} filter {input.pacbio_inq} --minchange 20 2> {log} | cut -f1-4 | grep -v '^#' > {output.catalog} 2>> {log}"""

rule inquiSTR_adotto_longtr:
    """Run inquiSTR call on PacBio using the LongTR-format catalog, so the filtered
    output catalog has LongTR-compatible 4th-field motifs."""
    input:
        catalog = "adotto_LongTR.bed",
        pacbio = "pacbio.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.inq.gz",
        timing = "tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.time"
    log:
        "logs/inquiSTR_adotto_longtr_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    threads:
        4
    resources:
        benchmark_slot=1
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.pacbio} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend 2> {log} | gzip > {output.inq} 2>> {log}
        """


rule filter_inquiSTR_pacbio_longtr:
    """Filter the LongTR-catalog inquiSTR PacBio calls to variable loci,
    producing a LongTR-compatible catalog (motif-only 4th field)."""
    input:
        pacbio_inq = "tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.inq.gz",
    output:
        catalog = "tool_comparison/adotto-variable-catalog-pacbio-longtr_rep{replicate}.bed",
        timing = "tool_comparison/filter-inquistr-pacbio-longtr_rep{replicate}.time"
    log:
        "logs/filter_inquistr_pacbio_longtr_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} filter {input.pacbio_inq} --minchange 20 2> {log} | cut -f1-4 | grep -v '^#' | grep -v '^chromosome' > {output.catalog} 2>> {log}"""


rule TRGT_adotto_filtered:
    input:
        catalog = "tool_comparison/adotto-variable-catalog-pacbio_rep{replicate}.bed",
        pacbio = "pacbio.cram",
    output:
        vcf = "tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.vcf.gz",
        timing = "tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.time"
    log:
        "logs/TRGT_adotto_filtered_rep{replicate}.log"
    params:
        TRGT = "/home/AD/wdecoster/bin/trgt",
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.TRGT} genotype --genome {params.reference} \
            --repeats {input.catalog} \
            --reads {input.pacbio} \
            --threads {threads} \
            --output-prefix tool_comparison/pacbio-trgt-adotto-filtered_rep{wildcards.replicate} &> {log}
        """

rule LongTR_pacbio:
    input:
        catalog = "adotto_LongTR.bed",
        pacbio = "pacbio.cram",
    output:
        vcf = "tool_comparison/pacbio-longtr-adotto_rep{replicate}.vcf.gz",
        timing = "tool_comparison/pacbio-longtr-adotto_rep{replicate}.time"
    log:
        "logs/LongTR_pacbio_rep{replicate}.log"
    params:
        LongTR = LongTR,
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.LongTR} \
            --bams {input.pacbio} \
            --fasta {params.reference} \
            --regions {input.catalog} \
            --tr-vcf {output.vcf} &> {log}
        """

rule LongTR_pacbio_filtered:
    input:
        catalog = "tool_comparison/adotto-variable-catalog-pacbio-longtr_rep{replicate}.bed",
        pacbio = "pacbio.cram",
    output:
        vcf = "tool_comparison/pacbio-longtr-adotto-filtered_rep{replicate}.vcf.gz",
        timing = "tool_comparison/pacbio-longtr-adotto-filtered_rep{replicate}.time"
    log:
        "logs/LongTR_pacbio_filtered_rep{replicate}.log"
    params:
        LongTR = LongTR,
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.LongTR} \
            --bams {input.pacbio} \
            --fasta {params.reference} \
            --regions {input.catalog} \
            --tr-vcf {output.vcf} &> {log}
        """

rule inquiSTR_ont:
    input:
        catalog = "adotto_TRGT.bed.gz",
        ont = "ont.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/ont-inquistr-adotto_rep{replicate}.inq.gz",
        timing = "tool_comparison/ont-inquistr-adotto_rep{replicate}.time"
    log:
        "logs/inquiSTR_ont_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    threads:
        4
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.ont} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend 2> {log} | gzip > {output.inq} 2>> {log}
        """

rule filter_inquiSTR_ont:
    input:
        ont_inq = "tool_comparison/ont-inquistr-adotto_rep{replicate}.inq.gz",
    output:
        catalog = "tool_comparison/adotto-variable-catalog-ont_rep{replicate}.bed",
        timing = "tool_comparison/filter-inquistr-ont_rep{replicate}.time"
    log:
        "logs/filter_inquistr_ont_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} filter {input.ont_inq} --minchange 20 2> {log} | cut -f1-4 > {output.catalog} 2>> {log}
        """

rule inquiSTR_ont_longtr:
    """Run inquiSTR call on ONT using the LongTR-format catalog, so the filtered
    output catalog has LongTR-compatible 4th-field motifs."""
    input:
        catalog = "adotto_LongTR.bed",
        ont = "ont.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.inq.gz",
        timing = "tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.time"
    log:
        "logs/inquiSTR_ont_longtr_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    threads:
        4
    resources:
        benchmark_slot=1
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.ont} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend 2> {log} | gzip > {output.inq} 2>> {log}
        """


rule filter_inquiSTR_ont_longtr:
    """Filter the LongTR-catalog inquiSTR ONT calls to variable loci,
    producing a LongTR-compatible catalog (motif-only 4th field)."""
    input:
        ont_inq = "tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.inq.gz"
    output:
        catalog = "tool_comparison/adotto-variable-catalog-ont-longtr_rep{replicate}.bed",
        timing = "tool_comparison/filter-inquistr-ont-longtr_rep{replicate}.time"
    log:
        "logs/filter_inquistr_ont_longtr_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        /usr/bin/time -v -o {output.timing} \\
        {params.inquiSTR} filter {input.ont_inq} --minchange 20 2> {log} | cut -f1-4 | grep -v '^#' | grep -v '^chromosome' > {output.catalog} 2>> {log}
        """

rule LongTR_ont:
    input:
        catalog = "adotto_LongTR.bed",
        ont = "ont.cram",
    output:
        vcf = "tool_comparison/ont-longtr-adotto_rep{replicate}.vcf.gz",
        timing = "tool_comparison/ont-longtr-adotto_rep{replicate}.time"
    log:
        "logs/LongTR_ont_rep{replicate}.log"
    params:
        LongTR = LongTR,
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.LongTR} \
            --bams {input.ont} \
            --fasta {params.reference} \
            --regions {input.catalog} \
            --tr-vcf {output.vcf} &> {log}
        """

rule LongTR_ont_filtered:
    input:
        catalog = "tool_comparison/adotto-variable-catalog-ont-longtr_rep{replicate}.bed",
        ont = "ont.cram",
    output:
        vcf = "tool_comparison/ont-longtr-adotto-filtered_rep{replicate}.vcf.gz",
        timing = "tool_comparison/ont-longtr-adotto-filtered_rep{replicate}.time"
    log:
        "logs/LongTR_ont_filtered_rep{replicate}.log"
    params:
        LongTR = LongTR,
        reference = reference
    resources:
        benchmark_slot=1  # Ensure only one benchmark runs at a time
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.LongTR} \
            --bams {input.ont} \
            --fasta {params.reference} \
            --regions {input.catalog} \
            --tr-vcf {output.vcf} &> {log}
        """

rule convert_trgt_pacbio:
    """
    TRGT coordinates should not be corrected between VCF and inquiSTR call format
    """
    input:
        vcf = "tool_comparison/pacbio-trgt-adotto_rep{replicate}.vcf.gz"
    output:
        inq = "tool_comparison/pacbio-trgt-adotto_rep{replicate}.inq.gz"
    log:
        "logs/convert_trgt_pacbio_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        {params.inquiSTR} convert --off-by-one {input.vcf} 2> {log} | gzip > {output.inq}
        """

rule convert_longtr_pacbio:
    input:
        vcf = "tool_comparison/pacbio-longtr-adotto_rep{replicate}.vcf.gz",
    output:
        inq = "tool_comparison/pacbio-longtr-adotto_rep{replicate}.inq.gz"
    log:
        "logs/convert_longtr_pacbio_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        {params.inquiSTR} convert --off-by-one {input.vcf} 2> {log} | gzip > {output.inq}
        """

rule convert_longtr_ont:
    input:
        vcf = "tool_comparison/ont-longtr-adotto_rep{replicate}.vcf.gz",
    output:
        inq = "tool_comparison/ont-longtr-adotto_rep{replicate}.inq.gz"
    log:
        "logs/convert_longtr_ont_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        {params.inquiSTR} convert --off-by-one {input.vcf} 2> {log} | gzip > {output.inq}
        """

rule aggregate_tool_comparison:
    input:
        # PacBio timing files
        pacbio_inquistr_time = expand("tool_comparison/pacbio-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        pacbio_trgt_time = expand("tool_comparison/pacbio-trgt-adotto_rep{replicate}.time", replicate=REPLICATES),
        pacbio_longtr_time = expand("tool_comparison/pacbio-longtr-adotto_rep{replicate}.time", replicate=REPLICATES),
        pacbio_filter_time = expand("tool_comparison/filter-inquistr-pacbio_rep{replicate}.time", replicate=REPLICATES),
        pacbio_trgt_filtered_time = expand("tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        pacbio_longtr_filtered_time = expand("tool_comparison/pacbio-longtr-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        pacbio_longtr_inquistr_time = expand("tool_comparison/pacbio-inquistr-adotto-longtr_rep{replicate}.time", replicate=REPLICATES),
        pacbio_longtr_filter_time = expand("tool_comparison/filter-inquistr-pacbio-longtr_rep{replicate}.time", replicate=REPLICATES),
        # ONT timing files
        ont_inquistr_time = expand("tool_comparison/ont-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        ont_longtr_time = expand("tool_comparison/ont-longtr-adotto_rep{replicate}.time", replicate=REPLICATES),
        ont_filter_time = expand("tool_comparison/filter-inquistr-ont_rep{replicate}.time", replicate=REPLICATES),
        ont_longtr_filtered_time = expand("tool_comparison/ont-longtr-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        ont_longtr_inquistr_time = expand("tool_comparison/ont-inquistr-adotto-longtr_rep{replicate}.time", replicate=REPLICATES),
        ont_longtr_filter_time = expand("tool_comparison/filter-inquistr-ont-longtr_rep{replicate}.time", replicate=REPLICATES)
    output:
        "tool_comparison/results.tsv"
    run:
        import re
        import pandas as pd

        def parse_timing_file(filepath):
            """Parse /usr/bin/time -v output file."""
            elapsed_time = None
            max_memory_kb = None
            with open(filepath, 'r') as f:
                for line in f:
                    if 'Elapsed (wall clock) time' in line:
                        time_str = line.split('): ', 1)[1].strip()
                        parts = time_str.split(':')
                        if len(parts) == 3:  # h:mm:ss
                            h, m, s = parts
                            elapsed_time = int(h) * 3600 + int(m) * 60 + float(s)
                        elif len(parts) == 2:  # mm:ss
                            m, s = parts
                            elapsed_time = int(m) * 60 + float(s)
                        elif len(parts) == 1:  # just ss
                            elapsed_time = float(parts[0])
                    elif 'Maximum resident set size' in line:
                        max_memory_kb = int(line.split(':')[1].strip())
            return elapsed_time, max_memory_kb

        def get_replicate(filepath):
            match = re.search(r'rep(\d+)', filepath)
            return int(match.group(1)) if match else None

        def add_single_tool(file_list, tool, technology, results):
            """Parse individual tool timing files and append to results."""
            for filepath in file_list:
                replicate = get_replicate(filepath)
                elapsed, memory = parse_timing_file(filepath)
                if elapsed and memory and replicate:
                    results.append({
                        'tool': tool,
                        'technology': technology,
                        'replicate': replicate,
                        'elapsed_seconds': elapsed,
                        'max_memory_gb': memory / (1024 * 1024)
                    })

        def add_combined_tool(inq_files, filter_files, tool_files, tool_label, technology, results):
            """Sum timing across three steps (inquiSTR + filter + tool) and append to results."""
            for inq_file, filt_file, tool_file in zip(inq_files, filter_files, tool_files):
                replicate = get_replicate(inq_file)
                inq_elapsed, inq_memory = parse_timing_file(inq_file)
                filt_elapsed, filt_memory = parse_timing_file(filt_file)
                tool_elapsed, tool_memory = parse_timing_file(tool_file)
                if all([inq_elapsed, filt_elapsed, tool_elapsed, replicate]):
                    results.append({
                        'tool': tool_label,
                        'technology': technology,
                        'replicate': replicate,
                        'elapsed_seconds': inq_elapsed + filt_elapsed + tool_elapsed,
                        'max_memory_gb': max(inq_memory, filt_memory, tool_memory) / (1024 * 1024)
                    })

        results = []

        # PacBio single tools
        add_single_tool(input.pacbio_inquistr_time, 'inquiSTR', 'pacbio', results)
        add_single_tool(input.pacbio_trgt_time, 'TRGT', 'pacbio', results)
        add_single_tool(input.pacbio_longtr_time, 'LongTR', 'pacbio', results)

        # PacBio combined pipelines
        add_combined_tool(
            input.pacbio_inquistr_time, input.pacbio_filter_time, input.pacbio_trgt_filtered_time,
            'inquiSTR+TRGT', 'pacbio', results
        )
        add_combined_tool(
            input.pacbio_longtr_inquistr_time, input.pacbio_longtr_filter_time, input.pacbio_longtr_filtered_time,
            'inquiSTR+LongTR', 'pacbio', results
        )

        # ONT single tools
        add_single_tool(input.ont_inquistr_time, 'inquiSTR', 'ont', results)
        add_single_tool(input.ont_longtr_time, 'LongTR', 'ont', results)

        # ONT combined pipeline
        add_combined_tool(
            input.ont_longtr_inquistr_time, input.ont_longtr_filter_time, input.ont_longtr_filtered_time,
            'inquiSTR+LongTR', 'ont', results
        )

        df = pd.DataFrame(results)
        df = df.sort_values(['technology', 'tool', 'replicate'])
        df.to_csv(output[0], sep='\t', index=False)

rule benchmark_inquistr_vs_trgt_pacbio:
    """
    Compare inquiSTR PacBio genotypes (--test) against TRGT PacBio genotypes converted to
    inquiSTR format (--truth). Uses rep1 as a representative call for accuracy comparison.
    """
    input:
        test = "tool_comparison/pacbio-inquistr-adotto_rep1.inq.gz",
        truth = "tool_comparison/pacbio-trgt-adotto_rep1.inq.gz",
    output:
        txt = "tool_comparison/benchmark_inquistr_vs_trgt_pacbio.tsv",
        plot = "tool_comparison/benchmark_inquistr_vs_trgt_pacbio.html",
        diff_out = "tool_comparison/benchmark_inquistr_vs_trgt_pacbio_discrepancies.tsv"
    params:
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS,
        tolerance = 3
    log:
        "logs/benchmark_inquistr_vs_trgt_pacbio.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --test {input.test} \
            --truth {input.truth} \
            --plot {output.plot} \
            --diff-out {output.diff_out} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            > {output.txt} 2> {log}
        """


rule inquiSTR_adotto_requirespanning:
    """Run inquiSTR call on PacBio with --require-spanning to restrict genotypes to
    loci fully covered by spanning reads, for direct comparison against TRGT."""
    input:
        catalog = "adotto_TRGT.bed.gz",
        pacbio = "pacbio.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/pacbio-inquistr-adotto-requirespanning.inq.gz"
    log:
        "logs/inquiSTR_adotto_requirespanning.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    threads:
        4
    shell:
        """
        {params.inquiSTR} call {input.pacbio} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend \
            --require-spanning 2> {log} | gzip > {output.inq} 2>> {log}
        """


rule benchmark_inquistr_requirespanning_vs_trgt_pacbio:
    """Compare inquiSTR PacBio --require-spanning genotypes against TRGT genotypes.
    This isolates spanning-read-only calls to assess how much accuracy improves
    when soft-clipped genotypes are excluded."""
    input:
        test = "tool_comparison/pacbio-inquistr-adotto-requirespanning.inq.gz",
        truth = "tool_comparison/pacbio-trgt-adotto_rep1.inq.gz",
    output:
        txt = "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio.tsv",
        plot = "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio.html",
        diff_out = "tool_comparison/benchmark_inquistr_requirespanning_vs_trgt_pacbio_discrepancies.tsv"
    params:
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS,
        tolerance = 3
    log:
        "logs/benchmark_inquistr_requirespanning_vs_trgt_pacbio.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --test {input.test} \
            --truth {input.truth} \
            --plot {output.plot} \
            --diff-out {output.diff_out} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            > {output.txt} 2> {log}
        """


rule benchmark_inquistr_vs_longtr_pacbio:
    """
    Compare inquiSTR PacBio genotypes (--test) against LongTR PacBio genotypes converted to
    inquiSTR format (--truth). Uses rep1 as a representative call for accuracy comparison.
    """
    input:
        test = "tool_comparison/pacbio-inquistr-adotto_rep1.inq.gz",
        truth = "tool_comparison/pacbio-longtr-adotto_rep1.inq.gz",
    output:
        txt = "tool_comparison/benchmark_inquistr_vs_longtr_pacbio.tsv",
        plot = "tool_comparison/benchmark_inquistr_vs_longtr_pacbio.html",
        diff_out = "tool_comparison/benchmark_inquistr_vs_longtr_pacbio_discrepancies.tsv"
    params:
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS
    log:
        "logs/benchmark_inquistr_vs_longtr_pacbio.log"
    shell:
        """
        {params.inquiSTR} benchmark \\
            --test {input.test} \\
            --truth {input.truth} \\
            --plot {output.plot} \\
            --diff-out {output.diff_out} \\
            --max-locus {params.max_locus} \\
            > {output.txt} 2> {log}
        """


rule benchmark_inquistr_vs_longtr_ont:
    """
    Compare inquiSTR ONT genotypes (--test) against LongTR ONT genotypes converted to
    inquiSTR format (--truth). Uses rep1 as a representative call for accuracy comparison.
    """
    input:
        test = "tool_comparison/ont-inquistr-adotto_rep1.inq.gz",
        truth = "tool_comparison/ont-longtr-adotto_rep1.inq.gz",
    output:
        txt = "tool_comparison/benchmark_inquistr_vs_longtr_ont.tsv",
        plot = "tool_comparison/benchmark_inquistr_vs_longtr_ont.html",
        diff_out = "tool_comparison/benchmark_inquistr_vs_longtr_ont_discrepancies.tsv"
    params:
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS
    log:
        "logs/benchmark_inquistr_vs_longtr_ont.log"
    shell:
        """
        {params.inquiSTR} benchmark \\
            --test {input.test} \\
            --truth {input.truth} \\
            --plot {output.plot} \\
            --diff-out {output.diff_out} \\
            --max-locus {params.max_locus} \\
            > {output.txt} 2> {log}
        """


rule inquiSTR_accuracy:
    """
    execute the inquiSTR benchmark analysis (named here accuracy because we already have a benchmark rule which is on time and memory)
    this command takes the inquiSTR call output and compares it to the adotto truth genotypes (in bed format)
    we will do this fo both ont and pacbio technologies, for a specific thread count and replicate (as that shouldn't matter for accuracy)
    """
    input:
        truth_bed = "/home/AD/wdecoster/optimize_inquiSTR/adotto/HG002_GRCh38_TandemRepeats_v1.0.bed.gz",
        genotypes = "benchmarking/data/{technology}_threads4_rep1.tsv",
        version = "inquiSTR_version.txt"
    output:
        txt = "benchmarking/accuracy_{technology}.tsv",
        plot = "benchmarking/accuracy_{technology}.html"
    params:
        inquiSTR = inquiSTR,
        tolerance = 3,
        max_locus = MAX_LOCUS # limit to same length as used in genotyping
    log:
        "logs/inquiSTR_accuracy_{technology}.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --truth {input.truth_bed} \
            --mode MAX \
            --tier1 \
            --plot {output.plot} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            --test {input.genotypes} \
            > {output.txt} 2> {log}
        """

rule pathogenic:
    input:
        inquistr = expand("puretarget-calls/{sample}.inq", sample=puretarget_files),
        combined = "puretarget-calls/combined.tsv",
        heatmap = "puretarget-calls/heatmap.html",

rule combine_puretarget:
    input:
        expand("puretarget-calls/{sample}.inq", sample=puretarget_files),
    output:
        "puretarget-calls/combined.tsv"
    log:
        "logs/combine_puretarget.log"
    threads: 4
    params:
        inquiSTR = inquiSTR,
    shell:
        """
        {params.inquiSTR} combine --threads {threads} {input} > {output} 2> {log}
        """

rule genotype_puretarget:
    input:
        bam = "/home/AD/wdecoster/inquiSTR_paper/puretarget-data/{sample}.bam",
    output:
        inq = "puretarget-calls/{sample}.inq",
    log:
        "logs/genotype_puretarget_{sample}.log"
    params:
        reference = reference,
        inquiSTR = inquiSTR,
    shell:
        """
        {params.inquiSTR} call --preset pathogenic --imbalance 0.1 --unphased {input.bam} > {output.inq} 2> {log}
        """


rule create_chr21_catalog_for_straglr:
    """Create a minimal STR catalog containing only chr21 loci."""
    input:
        catalog = "adotto_TRGT.bed.gz"
    output:
        catalog = "tool_comparison/straglr_chr21_catalog.bed"
    log:
        "logs/create_chr21_catalog_for_straglr.log"
    shell:
        """
        zcat {input.catalog} | awk 'BEGIN{{OFS="\t"}} $1=="chr21" {{
            n=split($4, fields, ";");
            motif="";
            for (i=1; i<=n; i++) {{
                if (fields[i] ~ /^STRUC=/) {{
                    motif=fields[i];
                    sub(/^STRUC=\\(/, "", motif);
                    sub(/\\)n.*$/, "", motif);
                    break;
                }}
            }}
            if (motif != "" && length(motif) >= 2 && length(motif) <= 50) print $1, $2, $3, motif;
        }}' > {output.catalog} 2> {log}
        """


rule run_straglr_chr21:
    """Run STRaglr on a chr21-only catalog for both PacBio and ONT."""
    input:
        catalog = "tool_comparison/straglr_chr21_catalog.bed",
        cram = "{technology}.cram"
    output:
        tsv = "tool_comparison/{technology}-straglr-chr21_rep{replicate}.tsv",
        bed = "tool_comparison/{technology}-straglr-chr21_rep{replicate}.bed",
        timing = "tool_comparison/{technology}-straglr-chr21_rep{replicate}.time"
    log:
        "logs/straglr_chr21_{technology}_rep{replicate}.log"
    params:
        straglr = STRAGLR,
        reference = reference
    resources:
        benchmark_slot=1
    threads:
        4
    conda:
        "envs/straglr.yml"
    shell:
        """
        prefix=$(echo {output.tsv} | sed 's/\\.tsv$//')
        /usr/bin/time -v -o {output.timing} \
        {params.straglr} {input.cram} {params.reference} $prefix \
            --loci {input.catalog} \
            --genotype_in_size \
            --nprocs {threads} &> {log}
        """


rule run_inquistr_chr21:
    """Run inquiSTR on the same chr21-only catalog for both PacBio and ONT."""
    input:
        catalog = "tool_comparison/straglr_chr21_catalog.bed",
        cram = "{technology}.cram",
        version = "inquiSTR_version.txt"
    output:
        inq = "tool_comparison/{technology}-inquistr-chr21_rep{replicate}.inq.gz",
        timing = "tool_comparison/{technology}-inquistr-chr21_rep{replicate}.time"
    log:
        "logs/inquistr_chr21_{technology}_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference,
        max_locus = MAX_LOCUS
    resources:
        benchmark_slot=1
    threads:
        4
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.cram} \
            --region-file {input.catalog} \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --noextend 2> {log} | gzip > {output.inq} 2>> {log}
        """


rule aggregate_straglr_chr21_results:
    input:
        straglr_time = expand("tool_comparison/{technology}-straglr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES),
        inquistr_time = expand("tool_comparison/{technology}-inquistr-chr21_rep{replicate}.time", technology=TECHNOLOGIES, replicate=REPLICATES)
    output:
        "tool_comparison/straglr_chr21_results.tsv"
    run:
        import re
        import pandas as pd

        def parse_timing_file(filepath):
            elapsed_time = None
            max_memory_kb = None
            with open(filepath, 'r') as f:
                for line in f:
                    if 'Elapsed (wall clock) time' in line:
                        time_str = line.split('): ', 1)[1].strip()
                        parts = time_str.split(':')
                        if len(parts) == 3:
                            h, m, s = parts
                            elapsed_time = int(h) * 3600 + int(m) * 60 + float(s)
                        elif len(parts) == 2:
                            m, s = parts
                            elapsed_time = int(m) * 60 + float(s)
                        elif len(parts) == 1:
                            elapsed_time = float(parts[0])
                    elif 'Maximum resident set size' in line:
                        max_memory_kb = int(line.split(':')[1].strip())
            return elapsed_time, max_memory_kb

        def parse_metadata(filepath):
            match = re.search(r'tool_comparison/(ont|pacbio)-(straglr|inquistr)-chr21_rep(\d+)\.time$', filepath)
            if not match:
                return None, None, None
            technology, tool, replicate = match.groups()
            tool_name = 'Straglr' if tool == 'straglr' else 'inquiSTR'
            return technology, tool_name, int(replicate)

        results = []
        for timing_file in input:
            technology, tool_name, replicate = parse_metadata(str(timing_file))
            if technology is None:
                continue
            elapsed, memory = parse_timing_file(timing_file)
            if elapsed is None or memory is None:
                continue
            results.append({
                'technology': technology,
                'tool': tool_name,
                'replicate': replicate,
                'elapsed_seconds': elapsed,
                'max_memory_gb': memory / (1024 * 1024)
            })

        df = pd.DataFrame(results)
        df = df.sort_values(['technology', 'tool', 'replicate'])
        df.to_csv(output[0], sep='\t', index=False)

rule convert_straglr_to_inquistr_format:
    """Convert STRaglr BED output to inquiSTR format for benchmarking.
    Removes header lines (starting with '#') and keeps columns 1-5 and 8.
    Adds inquiSTR-compatible header.
    """
    input:
        bed = "tool_comparison/{technology}-straglr-chr21_rep{replicate}.bed"
    output:
        inq = "tool_comparison/{technology}-straglr-chr21_rep{replicate}.inq.gz"
    log:
        "logs/convert_straglr_to_inquistr_{technology}_rep{replicate}.log"
    shell:
        """
           awk 'BEGIN {{OFS="\t"; print "chromosome", "begin", "end", "info", "straglr_H1", "straglr_H2"}} \
               $0 !~ /^#/ {{print $1, $2, $3, $4, $5, $8}}' {input.bed} \
             | gzip > {output.inq} 2> {log}
        """


rule benchmark_straglr_chr21_accuracy:
    """Benchmark STRaglr rep1 genotypes (rep1 for both ont and pacbio) against adotto truth."""
    input:
        test_inq = "tool_comparison/{technology}-straglr-chr21_rep1.inq.gz",
        truth_bed = "/home/AD/wdecoster/optimize_inquiSTR/adotto/HG002_GRCh38_TandemRepeats_v1.0.bed.gz",
        version = "inquiSTR_version.txt"
    output:
        txt = "tool_comparison/straglr_chr21_accuracy_{technology}.tsv",
        plot = "tool_comparison/straglr_chr21_accuracy_{technology}.html"
    params:
        inquiSTR = inquiSTR,
        tolerance = 3,
        max_locus = MAX_LOCUS
    log:
        "logs/benchmark_straglr_chr21_accuracy_{technology}.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --test {input.test_inq} \
            --truth {input.truth_bed} \
            --plot {output.plot} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            > {output.txt} 2> {log}
        """


rule benchmark_inquistr_chr21_accuracy:
    """Benchmark inquiSTR rep1 genotypes (rep1 for both ont and pacbio) against adotto truth."""
    input:
        test_inq = "tool_comparison/{technology}-inquistr-chr21_rep1.inq.gz",
        truth_bed = "/home/AD/wdecoster/optimize_inquiSTR/adotto/HG002_GRCh38_TandemRepeats_v1.0.bed.gz",
        version = "inquiSTR_version.txt"
    output:
        txt = "tool_comparison/inquistr_chr21_accuracy_{technology}.tsv",
        plot = "tool_comparison/inquistr_chr21_accuracy_{technology}.html"
    params:
        inquiSTR = inquiSTR,
        tolerance = 3,
        max_locus = MAX_LOCUS
    log:
        "logs/benchmark_inquistr_chr21_accuracy_{technology}.log"
    shell:
        """
        {params.inquiSTR} benchmark \
            --test {input.test_inq} \
            --truth {input.truth_bed} \
            --plot {output.plot} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            > {output.txt} 2> {log}
        """
