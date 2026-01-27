import pandas as pd

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
REPLICATES = [1, 2, 3]

# Randomize order of thread counts for benchmarking
import random
THREAD_ORDER = THREAD_COUNTS * len(REPLICATES)
random.shuffle(THREAD_ORDER)

rule all:
    input:
        "benchmarking/results.tsv",
        "benchmarking/runtime_plot.html",
        "benchmarking/memory_plot.html",
        "genotyping/adotto_combined_selected_samples.tsv",
        expand("tool_comparison/pacbio-trgt-adotto_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-trgt-adotto_rep{replicate}.time", replicate=REPLICATES),
        expand("tool_comparison/pacbio-inquistr-adotto_rep{replicate}.inq.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-inquistr-adotto_rep{replicate}.time", replicate=REPLICATES),
        expand("tool_comparison/adotto-variable-catalog_rep{replicate}.bed.gz", replicate=REPLICATES),
        expand("tool_comparison/filter-inquiSTR-adotto_rep{replicate}.time", replicate=REPLICATES),
        expand("tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.vcf.gz", replicate=REPLICATES),
        expand("tool_comparison/pacbio-trgt-adotto-filtered_rep{replicate}.time", replicate=REPLICATES),
        "tool_comparison/results.tsv",
        "tool_comparison/runtime_plot.html",
        "tool_comparison/memory_plot.html",
        "tool_versions.txt"

rule list_versions:
    input:
        "tool_versions.txt"

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
        bed = "polymorphic_repeats.hg38.bed"
    output:
        "polymorphic/repeats_combined.tsv"
    params:
        reference = reference,
        inquiSTR = inquiSTR
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
            --resume \
            --keep-going \
            > {log} 2>&1
        """


rule polymorphic_relate:
    input:
        "polymorphic/repeats_combined.tsv"
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
        "polymorphic/repeats_combined.tsv"
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
        "selected_samples_manifest.tsv"
    output:
        "genotyping/adotto_combined_selected_samples.tsv"
    params:
        reference = reference,
        inquiSTR = inquiSTR
    threads: 16
    log:
        "logs/genotype_adotto.log"
    shell:
        """
        {params.inquiSTR} batch {input} \
            --output {output} \
            --preset adotto \
            --unphased \
            --threads {threads} \
            --parallel-samples 4 \
            --reference {params.reference} \
            --resume \
            --save-individual genotyping/adotto_individual \
            > {log} 2>&1
        """


# Benchmark rules
rule benchmark_call:
    input:
        cram = "{technology}.cram",
        bed = "adotto.bed.gz"
    output:
        result = "benchmarking/data/{technology}_threads{threads}_rep{replicate}.tsv",
        timing = "benchmarking/data/{technology}_threads{threads}_rep{replicate}.time"
    params:
        inquiSTR = inquiSTR,
        reference = reference
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

rule plot_benchmark_results:
    input:
        "benchmarking/results.tsv"
    output:
        "benchmarking/runtime_plot.html"
    run:
        import pandas as pd
        import plotly.graph_objects as go
        
        # Read the data
        df = pd.read_csv(input[0], sep='\t')
        
        # Convert seconds to minutes
        df['elapsed_minutes'] = df['elapsed_seconds'] / 60
        
        # Calculate mean per technology and thread count
        mean_df = df.groupby(['technology', 'threads'])['elapsed_minutes'].mean().reset_index()
        
        # Color mapping
        colors = {
            'pacbio': 'purple',
            'ont': 'blue'
        }
        
        # Create the plot
        fig = go.Figure()
        
        # Add individual points and lines for each technology
        for tech in df['technology'].unique():
            tech_data = df[df['technology'] == tech]
            tech_mean = mean_df[mean_df['technology'] == tech]
            
            # Add scatter points for individual measurements
            fig.add_trace(go.Scatter(
                x=tech_data['threads'],
                y=tech_data['elapsed_minutes'],
                mode='markers',
                name=f'{tech} (replicates)',
                marker=dict(
                    color=colors[tech],
                    size=8,
                    opacity=0.5
                ),
                showlegend=False
            ))
            
            # Add line for mean values
            fig.add_trace(go.Scatter(
                x=tech_mean['threads'],
                y=tech_mean['elapsed_minutes'],
                mode='lines+markers',
                name=f'{tech}',
                line=dict(
                    color=colors[tech],
                    width=3
                ),
                marker=dict(
                    color=colors[tech],
                    size=10
                ),
                showlegend=True
            ))
        
        # Update layout
        fig.update_layout(
            title='Benchmark Results: Runtime vs Thread Count',
            xaxis_title='Number of Threads',
            yaxis_title='Elapsed Time (minutes)',
            plot_bgcolor='white',
            hovermode='closest',
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99
            )
        )
        
        # Update axes
        fig.update_xaxes(
            showgrid=True,
            gridcolor='lightgray',
            showline=True,
            linewidth=2,
            linecolor='black'
        )
        fig.update_yaxes(
            showgrid=True,
            gridcolor='lightgray',
            showline=True,
            linewidth=2,
            linecolor='black'
        )
        
        # Save the plot
        fig.write_html(output[0])

rule plot_memory_usage:
    input:
        "benchmarking/results.tsv"
    output:
        "benchmarking/memory_plot.html"
    run:
        import pandas as pd
        import plotly.graph_objects as go
        
        # Read the data
        df = pd.read_csv(input[0], sep='\t')
        
        # Calculate mean per technology and thread count
        mean_df = df.groupby(['technology', 'threads'])['max_memory_gb'].mean().reset_index()
        
        # Color mapping
        colors = {
            'pacbio': 'purple',
            'ont': 'blue'
        }
        
        # Create the plot
        fig = go.Figure()
        
        # Add individual points and lines for each technology
        for tech in df['technology'].unique():
            tech_data = df[df['technology'] == tech]
            tech_mean = mean_df[mean_df['technology'] == tech]
            
            # Add scatter points for individual measurements
            fig.add_trace(go.Scatter(
                x=tech_data['threads'],
                y=tech_data['max_memory_gb'],
                mode='markers',
                name=f'{tech} (replicates)',
                marker=dict(
                    color=colors[tech],
                    size=8,
                    opacity=0.5
                ),
                showlegend=False
            ))
            
            # Add line for mean values
            fig.add_trace(go.Scatter(
                x=tech_mean['threads'],
                y=tech_mean['max_memory_gb'],
                mode='lines+markers',
                name=f'{tech}',
                line=dict(
                    color=colors[tech],
                    width=3
                ),
                marker=dict(
                    color=colors[tech],
                    size=10
                ),
                showlegend=True
            ))
        
        # Update layout
        fig.update_layout(
            title='Benchmark Results: Memory Usage vs Thread Count',
            xaxis_title='Number of Threads',
            yaxis_title='Maximum Memory Usage (GB)',
            plot_bgcolor='white',
            hovermode='closest',
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99
            )
        )
        
        # Update axes
        fig.update_xaxes(
            showgrid=True,
            gridcolor='lightgray',
            showline=True,
            linewidth=2,
            linecolor='black'
        )
        fig.update_yaxes(
            showgrid=True,
            gridcolor='lightgray',
            showline=True,
            linewidth=2,
            linecolor='black'
        )
        
        # Save the plot
        fig.write_html(output[0])

rule download_adotto:
    output:
        catalog = "adotto_TRGT.bed.gz",
    log:
        "logs/download_adotto.log"
    params:
        url = "https://zenodo.org/records/13987414/files/adotto_TRregions_v1.2.1.bed.gz",
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
    output:
        inq = "tool_comparison/pacbio-inquistr-adotto_rep{replicate}.inq.gz",
        timing = "tool_comparison/pacbio-inquistr-adotto_rep{replicate}.time"
    log:
        "logs/inquiSTR_adotto_rep{replicate}.log"
    params:
        inquiSTR = inquiSTR,
        reference = reference
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
            --unphased | gzip > {output.inq} 2> {log}
        """

rule filter_inquiSTR_adotto:
    input:
        "tool_comparison/pacbio-inquistr-adotto_rep{replicate}.inq.gz"
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
        {params.inquiSTR} filter {input} --minchange 20 | cut -f1-4 | gzip > {output.catalog} 2> {log}"""

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

rule plot_tool_comparison_time:
    input:
        "tool_comparison/results.tsv"
    output:
        "tool_comparison/runtime_plot.html"
    run:
        import pandas as pd
        import plotly.graph_objects as go
        
        # Read the data
        df = pd.read_csv(input[0], sep='\t')
        
        # Convert seconds to minutes
        df['elapsed_minutes'] = df['elapsed_seconds'] / 60
        
        # Calculate mean per tool
        mean_df = df.groupby('tool')['elapsed_minutes'].mean().reset_index()
        mean_df = mean_df.sort_values('tool')  # Ensure consistent ordering
        
        # Color mapping
        colors = {
            'inquiSTR': '#1f77b4',
            'TRGT': '#ff7f0e',
            'inquiSTR+TRGT': '#2ca02c'
        }
        
        # Create bar plot with individual points
        fig = go.Figure()
        
        # Add individual replicate points as scatter
        for tool in df['tool'].unique():
            tool_data = df[df['tool'] == tool]
            fig.add_trace(go.Scatter(
                x=[tool] * len(tool_data),
                y=tool_data['elapsed_minutes'],
                mode='markers',
                name=f'{tool} (replicates)',
                marker=dict(
                    color=colors[tool],
                    size=10,
                    opacity=0.6,
                    line=dict(width=1, color='white')
                ),
                showlegend=False
            ))
        
        # Add mean bars
        fig.add_trace(go.Bar(
            x=mean_df['tool'],
            y=mean_df['elapsed_minutes'],
            marker_color=[colors[tool] for tool in mean_df['tool']],
            text=[f"{val:.2f} min" for val in mean_df['elapsed_minutes']],
            textposition='outside',
            name='Mean',
            opacity=0.7
        ))
        
        # Update layout
        fig.update_layout(
            title='Tool Comparison: Runtime (mean with individual replicates)',
            xaxis_title='Tool',
            yaxis_title='Elapsed Time (minutes)',
            plot_bgcolor='white',
            showlegend=False
        )
        
        # Update axes
        fig.update_xaxes(
            showgrid=False,
            showline=True,
            linewidth=2,
            linecolor='black'
        )
        fig.update_yaxes(
            showgrid=True,
            gridcolor='lightgray',
            showline=True,
            linewidth=2,
            linecolor='black'
        )
        
        # Save the plot
        fig.write_html(output[0])

rule plot_tool_comparison_memory:
    input:
        "tool_comparison/results.tsv"
    output:
        "tool_comparison/memory_plot.html"
    run:
        import pandas as pd
        import plotly.graph_objects as go
        
        # Read the data
        df = pd.read_csv(input[0], sep='\t')
        
        # Calculate mean per tool
        mean_df = df.groupby('tool')['max_memory_gb'].mean().reset_index()
        mean_df = mean_df.sort_values('tool')  # Ensure consistent ordering
        
        # Color mapping
        colors = {
            'inquiSTR': '#1f77b4',
            'TRGT': '#ff7f0e',
            'inquiSTR+TRGT': '#2ca02c'
        }
        
        # Create bar plot with individual points
        fig = go.Figure()
        
        # Add individual replicate points as scatter
        for tool in df['tool'].unique():
            tool_data = df[df['tool'] == tool]
            fig.add_trace(go.Scatter(
                x=[tool] * len(tool_data),
                y=tool_data['max_memory_gb'],
                mode='markers',
                name=f'{tool} (replicates)',
                marker=dict(
                    color=colors[tool],
                    size=10,
                    opacity=0.6,
                    line=dict(width=1, color='white')
                ),
                showlegend=False
            ))
        
        # Add mean bars
        fig.add_trace(go.Bar(
            x=mean_df['tool'],
            y=mean_df['max_memory_gb'],
            marker_color=[colors[tool] for tool in mean_df['tool']],
            text=[f"{val:.2f} GB" for val in mean_df['max_memory_gb']],
            textposition='outside',
            name='Mean',
            opacity=0.7
        ))
        
        # Update layout
        fig.update_layout(
            title='Tool Comparison: Memory Usage (mean with individual replicates)',
            xaxis_title='Tool',
            yaxis_title='Maximum Memory Usage (GB)',
            plot_bgcolor='white',
            showlegend=False
        )
        
        # Update axes
        fig.update_xaxes(
            showgrid=False,
            showline=True,
            linewidth=2,
            linecolor='black'
        )
        fig.update_yaxes(
            showgrid=True,
            gridcolor='lightgray',
            showline=True,
            linewidth=2,
            linecolor='black'
        )
        
        # Save the plot
        fig.write_html(output[0])

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
