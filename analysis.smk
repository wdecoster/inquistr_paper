import pandas as pd

# Include plotting rules from separate file
include: "plotting.smk"

file_path = "/home/AD/wdecoster/inquiSTR_paper/1000G_cohort.tsv"
# Load a TSV file into a DataFrame, containing URLs of 1000 Genomes cram files
df = pd.read_csv(file_path, sep='\t', usecols=['sample', 'hg38_path', 'source'])
df = df[df["source"] == "Gustafson"]

# get a list of 100 samples
selected_samples = df[:100]


reference = "/home/AD/wdecoster/database/GRCh38.fa"
inquiSTR = "/home/AD/wdecoster/repositories/inquiSTR/target/x86_64-unknown-linux-musl/release/inquiSTR"

# Benchmark parameters
TECHNOLOGIES = ["ont", "pacbio"]
THREAD_COUNTS = list(range(1, 13))  # 1 to 12 threads
REPLICATES = list(range(1,6))  # 5 replicates
MAX_LOCUS = 10000  # limit to loci shorter than 10kb for genotyping

# Randomize order of thread counts for benchmarking
import random
THREAD_ORDER = THREAD_COUNTS * len(REPLICATES)
random.shuffle(THREAD_ORDER)

# if the "tool_versions.txt" file exists, remove it to ensure fresh capture of tool versions
import os
if os.path.exists("tool_versions.txt"):
    os.remove("tool_versions.txt")

rule all:
    input:
        "benchmarking/results.tsv",
        "benchmarking/runtime_plot.html",
        "benchmarking/memory_plot.html",
        #"genotyping/adotto_combined_selected_samples.tsv",
        expand("tool_comparison/pacbio-trgt-adotto_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-trgt-adotto_rep{replicate}.time", replicate=REPLICATES),
        expand("tool_comparison/pacbio-inquistr-adotto_rep{replicate}.inq.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        expand("tool_comparison/adotto-variable-catalog_rep{replicate}.bed.gz", replicate=REPLICATES),
        expand("tool_comparison/filter-inquiSTR-adotto_rep{replicate}.time", replicate=REPLICATES),
        expand("tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        expand("benchmarking/accuracy_{technology}.tsv", technology=TECHNOLOGIES),
        "tool_comparison/results.tsv",
        "tool_comparison/runtime_plot.html",
        "tool_comparison/memory_plot.html",
        "tool_versions.txt",
        "inquiSTR_version.txt"

rule versions:
    input:
        "tool_versions.txt"

rule inquiSTR_version:
    """Capture inquiSTR version to trigger reruns when version changes"""
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

# using polymorphic repeats from illumina https://zenodo.org/records/8329210/files/polymorphic_repeats.hg38.bed?download=1
rule polymorphic:
    input:
        "polymorphic/repeats_relate.tsv",
        "polymorphic/repeats_pca.html",

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

rule create_polymorphic_manifest:
    output:
        manifest = "polymorphic_manifest.tsv"
    log:
        "logs/create_polymorphic_manifest.log"
    run:
        with open(output.manifest, 'w') as f:
            f.write("bam_path\tsample_name\n")
            for idx, row in df.iterrows():
                f.write(f"{row['hg38_path']}\t{row['sample']}\n")


rule genotype_polymorphic:
    input:
        manifest = "polymorphic_manifest.tsv",
        bed = "polymorphic_repeats.hg38.bed",
        version = "inquiSTR_version.txt"
    output:
        "polymorphic/repeats_combined.tsv"
    params:
        reference = reference,
        inquiSTR = inquiSTR,
        max_locus = MAX_LOCUS # limit to loci shorter than 10kb for genotyping
    threads: 4
    log:
        "logs/genotype_polymorphic.log"
    shell:
        """
        {params.inquiSTR} batch {input.manifest} \
            --output {output} \
            --region-file {input.bed} \
            --unphased \
            --threads {threads} \
            --reference {params.reference} \
            --max-locus {params.max_locus} \
            --resume \
            --keep-going \
            > {log} 2>&1
        """


rule polymorphic_relate:
    input:
        "polymorphic/repeats_combined.tsv",
        version = "inquiSTR_version.txt"
    output:
        "polymorphic/repeats_relate.tsv"
    threads: 16
    params:
        inquiSTR = inquiSTR
    log:
        "logs/polymorphic_relate.log"
    shell:
        """
        {params.inquiSTR} relate \
            --output {output} \
            --threads {threads} \
            {input} \
            > {log} 2>&1
        """


rule polymorphic_pca:
    input:
        "polymorphic/repeats_combined.tsv",
        version = "inquiSTR_version.txt"
    output:
        "polymorphic/repeats_pca.html"
    threads: 16
    params:
        inquiSTR = inquiSTR
    log:
        "logs/polymorphic_pca.log"
    shell:
        """
        {params.inquiSTR} pca \
            --output {output} \
            --threads {threads} \
            {input} \
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
            --unphased \
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
            --unphased 2> {log} | gzip > {output.inq} 2>> {log}
        """

rule filter_inquiSTR_adotto:
    input:
        pacbio_inq = "tool_comparison/pacbio-inquistr-adotto_rep{replicate}.inq.gz",
        version = "inquiSTR_version.txt"
    output:
        catalog = "tool_comparison/adotto-variable-catalog_rep{replicate}.bed.gz",
        timing = "tool_comparison/filter-inquiSTR-adotto_rep{replicate}.time"
    log:
        "logs/filter_inquiSTR_adotto_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} filter {input.pacbio_inq} --minchange 20 2> {log} | cut -f1-4 | gzip > {output.catalog} 2>> {log}"""

rule TRGT_adotto_filtered:
    input:
        catalog = "tool_comparison/adotto-variable-catalog_rep{replicate}.bed.gz",
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

rule aggregate_tool_comparison:
    input:
        inquistr_time = expand("tool_comparison/pacbio-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        trgt_time = expand("tool_comparison/pacbio-trgt-adotto_rep{replicate}.time", replicate=REPLICATES),
        filter_time = expand("tool_comparison/filter-inquiSTR-adotto_rep{replicate}.time", replicate=REPLICATES),
        trgt_filtered_time = expand("tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.time", replicate=REPLICATES)
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
        
        results = []
        
        # Parse inquiSTR timing for each replicate
        for filepath in input.inquistr_time:
            match = re.search(r'rep(\d+)', filepath)
            replicate = int(match.group(1)) if match else None
            elapsed, memory = parse_timing_file(filepath)
            if elapsed and memory and replicate:
                results.append({
                    'tool': 'inquiSTR',
                    'replicate': replicate,
                    'elapsed_seconds': elapsed,
                    'max_memory_gb': memory / (1024 * 1024)
                })
        
        # Parse TRGT timing for each replicate
        for filepath in input.trgt_time:
            match = re.search(r'rep(\d+)', filepath)
            replicate = int(match.group(1)) if match else None
            elapsed, memory = parse_timing_file(filepath)
            if elapsed and memory and replicate:
                results.append({
                    'tool': 'TRGT',
                    'replicate': replicate,
                    'elapsed_seconds': elapsed,
                    'max_memory_gb': memory / (1024 * 1024)
                })
        
        # Parse inquiSTR + filter + TRGT filtered timing (sum of all three) for each replicate
        for i, (inq_file, filt_file, trgt_file) in enumerate(zip(input.inquistr_time, input.filter_time, input.trgt_filtered_time)):
            match = re.search(r'rep(\d+)', inq_file)
            replicate = int(match.group(1)) if match else None
            
            inq_elapsed, inq_memory = parse_timing_file(inq_file)
            filt_elapsed, filt_memory = parse_timing_file(filt_file)
            trgt_filt_elapsed, trgt_filt_memory = parse_timing_file(trgt_file)
            
            if all([inq_elapsed, filt_elapsed, trgt_filt_elapsed, replicate]):
                total_elapsed = inq_elapsed + filt_elapsed + trgt_filt_elapsed
                # For memory, use the maximum of the three steps
                max_memory = max(inq_memory, filt_memory, trgt_filt_memory)
                results.append({
                    'tool': 'inquiSTR+TRGT',
                    'replicate': replicate,
                    'elapsed_seconds': total_elapsed,
                    'max_memory_gb': max_memory / (1024 * 1024)
                })
        
        # Write results
        df = pd.DataFrame(results)
        df = df.sort_values(['tool', 'replicate'])
        df.to_csv(output[0], sep='\t', index=False)
        
        if all([inq_elapsed, filt_elapsed, trgt_filt_elapsed]):
            total_elapsed = inq_elapsed + filt_elapsed + trgt_filt_elapsed
            # For memory, use the maximum of the three steps
            max_memory = max(inq_memory, filt_memory, trgt_filt_memory)
            results.append({
                'tool': 'inquiSTR+TRGT',
                'elapsed_seconds': total_elapsed,
                'max_memory_gb': max_memory / (1024 * 1024)
            })
        
        # Write results
        df = pd.DataFrame(results)
        df.to_csv(output[0], sep='\t', index=False)

rule capture_tool_versions:
    output:
        "tool_versions.txt"
    params:
        inquiSTR = inquiSTR,
        TRGT = "/home/AD/wdecoster/bin/trgt"
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
            echo "Generated on: $(date)"
        }} > {output} 2> {log}
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
            --bed {input.truth_bed} \
            --mode MAX \
            --tier1 \
            --plot {output.plot} \
            --tolerance {params.tolerance} \
            --max-locus {params.max_locus} \
            {input.genotypes} \
            > {output.txt} 2> {log}
        """
