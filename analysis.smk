import pandas as pd

file_path = "/home/AD/wdecoster/inquiSTR_paper/1000G_cohort.tsv"
# Load a TSV file into a DataFrame, containing URLs of 1000 Genomes cram files
df = pd.read_csv(file_path, sep='\t', usecols=['sample', 'hg38_path', 'source'])
df = df[df["source"] == "Gustafson"]
reference = "/home/AD/wdecoster/database/GRCh38s.fa"
inquiSTR = "/home/AD/wdecoster/repositories/inquiSTR/target/x86_64-unknown-linux-musl/release/inquiSTR"

# Benchmark parameters
TECHNOLOGIES = ["ont", "pacbio"]
THREAD_COUNTS = list(range(1, 13))  # 1 to 12 threads
REPLICATES = [1, 2, 3]

# Randomize order of thread counts for benchmarking
import random
THREAD_ORDER = THREAD_COUNTS * len(REPLICATES)
random.shuffle(THREAD_ORDER)

rule all:
    input:
        "codis_relate.tsv",
        "codis_pca.html",
        "benchmark_results.tsv"

rule create_codis_manifest:
    output:
        manifest = "codis_manifest.tsv"
    log:
        "logs/create_codis_manifest.log"
    run:
        with open(output.manifest, 'w') as f:
            f.write("bam_path\tsample_name\n")
            for idx, row in df.iterrows():
                f.write(f"{row['hg38_path']}\t{row['sample']}\n")


rule genotype_codis:
    input:
        "codis_manifest.tsv"
    output:
        "codis_combined.tsv"
    params:
        reference = reference,
        inquiSTR = inquiSTR
    threads: 8
    log:
        "logs/genotype_codis.log"
    shell:
        """
        {params.inquiSTR} batch {input} \
            --output {output} \
            --preset codis \
            --unphased \
            --threads {threads} \
            --reference {params.reference} \
            --resume \
            > {log} 2>&1
        """


rule codis_relate:
    input:
        "codis_combined.tsv"
    output:
        "codis_relate.tsv"
    threads: 16
    params:
        inquiSTR = inquiSTR
    log:
        "logs/codis_relate.log"
    shell:
        """
        {params.inquiSTR} relate \
            --combined {input} \
            --output {output} \
            --threads {threads} \
            > {log} 2>&1
        """


rule codis_pca:
    input:
        "codis_combined.tsv"
    output:
        "codis_pca.html"
    threads: 16
    params:
        inquiSTR = inquiSTR
    log:
        "logs/codis_pca.log"
    shell:
        """
        {params.inquiSTR} pca \
            --combined {input} \
            --output {output} \
            --threads {threads} \
            > {log} 2>&1
        """


# Benchmark rules
rule benchmark_call:
    input:
        cram = "{technology}.cram",
        bed = "adotto.bed.gz"
    output:
        result = "benchmarks/{technology}_threads{threads}_rep{replicate}.tsv",
        timing = "benchmarks/{technology}_threads{threads}_rep{replicate}.time"
    params:
        inquiSTR = inquiSTR,
        reference = reference
    threads: lambda wildcards: int(wildcards.threads)
    log:
        "logs/benchmark_{technology}_threads{threads}_rep{replicate}.log"
    shell:
        """
        /usr/bin/time -v -o {output.timing} \
        {params.inquiSTR} call {input.cram} \
            --region-file {input.bed} \
            --threads {threads} \
            --reference {params.reference} \
            > {output.result} 2> {log}
        """


rule aggregate_benchmark_results:
    input:
        expand("benchmarks/{technology}_threads{threads}_rep{replicate}.time",
               technology=TECHNOLOGIES,
               threads=THREAD_COUNTS,
               replicate=REPLICATES)
    output:
        "benchmark_results.tsv"
    run:
        import re
        results = []
        
        for timing_file in input:
            # Parse filename to get metadata
            match = re.search(r'benchmarks/(\w+)_threads(\d+)_rep(\d+)\.time', timing_file)
            if match:
                technology, threads, replicate = match.groups()
                
                # Parse timing output from /usr/bin/time -v
                elapsed_time = None
                with open(timing_file, 'r') as f:
                    for line in f:
                        if 'Elapsed (wall clock) time' in line:
                            # Format is "h:mm:ss.ss" or "m:ss.ss"
                            time_str = line.split(':', 1)[1].strip()
                            # Convert to seconds
                            parts = time_str.split(':')
                            if len(parts) == 3:  # h:mm:ss
                                h, m, s = parts
                                elapsed_time = int(h) * 3600 + int(m) * 60 + float(s)
                            elif len(parts) == 2:  # mm:ss
                                m, s = parts
                                elapsed_time = int(m) * 60 + float(s)
                            break
                
                if elapsed_time is not None:
                    results.append({
                        'technology': technology,
                        'threads': int(threads),
                        'replicate': int(replicate),
                        'elapsed_seconds': elapsed_time
                    })
        
        # Write results
        import pandas as pd
        df = pd.DataFrame(results)
        df = df.sort_values(['technology', 'threads', 'replicate'])
        df.to_csv(output[0], sep='\t', index=False)
