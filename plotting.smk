"""
Plotting rules for inquiSTR paper analysis
"""

file_path = "/home/AD/wdecoster/inquiSTR_paper/1000G_cohort.tsv"


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
                    size=12,
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
                    size=12
                ),
                showlegend=True
            ))
        
        # Update layout
        fig.update_layout(
            title='Runtime vs. threads',
            xaxis_title='Number of threads',
            yaxis_title='Run time (minutes)',
            plot_bgcolor='white',
            hovermode='closest',
            font=dict(size=30),
            title_font_size=42,
            width=700,
            height=900,
            hoverlabel=dict(font_size=24),
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99,
                font=dict(size=30),
            )
        )
        
        # Update axes
        fig.update_xaxes(
            showgrid=True,
            gridcolor='lightgray',
            showline=True,
            linewidth=2,
            linecolor='black',
            title_font_size=34,
            tickfont_size=28,
        )
        fig.update_yaxes(
            showgrid=True,
            gridcolor='lightgray',
            showline=True,
            linewidth=2,
            linecolor='black',
            title_font_size=34,
            tickfont_size=28,
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
                    size=12,
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
                    width=5
                ),
                marker=dict(
                    color=colors[tech],
                    size=12
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
            font=dict(size=30),
            title_font_size=42,
            hoverlabel=dict(font_size=24),
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99,
                font=dict(size=30),
            )
        )
        
        # Update axes
        fig.update_xaxes(
            showgrid=True,
            gridcolor='lightgray',
            showline=True,
            linewidth=2,
            linecolor='black',
            title_font_size=34,
            tickfont_size=28,
        )
        fig.update_yaxes(
            showgrid=True,
            gridcolor='lightgray',
            showline=True,
            linewidth=2,
            linecolor='black',
            title_font_size=34,
            tickfont_size=28,
        )
        
        # Save the plot
        fig.write_html(output[0])

rule plot_tool_comparison_time:
    input:
        "tool_comparison/results.tsv"
    output:
        "tool_comparison/runtime_plot.html"
    run:
        import pandas as pd
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        # Read the data
        df = pd.read_csv(input[0], sep='\t')
        df['elapsed_minutes'] = df['elapsed_seconds'] / 60

        # Tools per subplot
        subplots = [
            ('PacBio', 'pacbio', ['inquiSTR', 'TRGT', 'inquiSTR+TRGT']),
            ('ONT', 'ont', ['inquiSTR', 'LongTR', 'inquiSTR+LongTR']),
        ]

        tool_colors = {
            'inquiSTR':       '#1f77b4',
            'TRGT':           '#ff7f0e',
            'inquiSTR+TRGT':  '#d62728',
            'LongTR':         '#2ca02c',
            'inquiSTR+LongTR':'#9467bd',
        }

        fig = make_subplots(rows=1, cols=2, subplot_titles=[s[0] for s in subplots],
                            horizontal_spacing=0.15)

        tool_seen = set()
        for col_i, (label, tech, tools) in enumerate(subplots, start=1):
            tech_df = df[df['technology'] == tech]
            present_tools = [t for t in tools if t in tech_df['tool'].values]
            n = len(present_tools)
            group_width = 0.8
            bw = group_width / n

            for j, tool in enumerate(present_tools):
                offset = (j - (n - 1) / 2) * bw
                x = offset
                color = tool_colors[tool]
                tool_df = tech_df[tech_df['tool'] == tool]
                mean_val = tool_df['elapsed_minutes'].mean()

                fig.add_trace(go.Bar(
                    x=[x], y=[mean_val],
                    name=tool, marker_color=color, opacity=0.8,
                    width=bw * 0.9,
                    text=[f"{mean_val:.1f} min"],
                    textposition='outside',
                    textfont=dict(size=16),
                    legendgroup=tool,
                    showlegend=(tool not in tool_seen),
                ), row=1, col=col_i)
                tool_seen.add(tool)

                fig.add_trace(go.Scatter(
                    x=[x] * len(tool_df),
                    y=tool_df['elapsed_minutes'],
                    mode='markers',
                    marker=dict(color=color, size=12, opacity=0.6,
                                line=dict(width=1, color='white')),
                    legendgroup=tool, showlegend=False,
                ), row=1, col=col_i)

        fig.update_layout(
            title='Tool Comparison: Runtime',
            title_x=0.5,
            barmode='overlay',
            plot_bgcolor='white',
            font=dict(size=20),
            title_font_size=24,
            width=600,
            height=600,
            showlegend=False,
        )
        # Make subplot titles larger
        for annotation in fig['layout']['annotations']:
            annotation['font'] = dict(size=22)

        for col_i, (label, tech, tools) in enumerate(subplots, start=1):
            present_tools = [t for t in tools if t in df[df['technology'] == tech]['tool'].values]
            n = len(present_tools)
            bw = 0.8 / n
            tick_positions = [(j - (n - 1) / 2) * bw for j in range(n)]
            xaxis = f"xaxis{col_i}" if col_i > 1 else "xaxis"
            yaxis = f"yaxis{col_i}" if col_i > 1 else "yaxis"
            fig.update_layout(**{
                xaxis: dict(
                    tickmode='array', tickvals=tick_positions, ticktext=present_tools,
                    tickfont=dict(size=16), showgrid=False, showline=True,
                    linewidth=2, linecolor='black',
                ),
                yaxis: dict(
                    title='Run time (minutes)', title_font_size=20, tickfont_size=18,
                    showgrid=True, gridcolor='lightgray', showline=True,
                    linewidth=2, linecolor='black', zeroline=False,
                ),
            })

        fig.write_html(output[0])

rule plot_tool_comparison_memory:
    input:
        "tool_comparison/results.tsv"
    output:
        "tool_comparison/memory_plot.html"
    run:
        import pandas as pd
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        # Read the data
        df = pd.read_csv(input[0], sep='\t')

        # Tools per subplot
        subplots = [
            ('PacBio', 'pacbio', ['inquiSTR', 'TRGT', 'inquiSTR+TRGT']),
            ('ONT', 'ont', ['inquiSTR', 'LongTR', 'inquiSTR+LongTR']),
        ]

        tool_colors = {
            'inquiSTR':       '#1f77b4',
            'TRGT':           '#ff7f0e',
            'inquiSTR+TRGT':  '#d62728',
            'LongTR':         '#2ca02c',
            'inquiSTR+LongTR':'#9467bd',
        }

        fig = make_subplots(rows=1, cols=2, subplot_titles=[s[0] for s in subplots],
                            horizontal_spacing=0.15)

        tool_seen = set()
        for col_i, (label, tech, tools) in enumerate(subplots, start=1):
            tech_df = df[df['technology'] == tech]
            present_tools = [t for t in tools if t in tech_df['tool'].values]
            n = len(present_tools)
            group_width = 0.8
            bw = group_width / n

            for j, tool in enumerate(present_tools):
                offset = (j - (n - 1) / 2) * bw
                x = offset
                color = tool_colors[tool]
                tool_df = tech_df[tech_df['tool'] == tool]
                mean_val = tool_df['max_memory_gb'].mean()

                fig.add_trace(go.Bar(
                    x=[x], y=[mean_val],
                    name=tool, marker_color=color, opacity=0.8,
                    width=bw * 0.9,
                    text=[f"{mean_val:.1f} GB"],
                    textposition='outside',
                    textfont=dict(size=16),
                    legendgroup=tool,
                    showlegend=False,
                ), row=1, col=col_i)
                tool_seen.add(tool)

                fig.add_trace(go.Scatter(
                    x=[x] * len(tool_df),
                    y=tool_df['max_memory_gb'],
                    mode='markers',
                    marker=dict(color=color, size=12, opacity=0.6,
                                line=dict(width=1, color='white')),
                    legendgroup=tool, showlegend=False,
                ), row=1, col=col_i)

        fig.update_layout(
            title='Tool Comparison: Memory Usage',
            title_x=0.5,
            barmode='overlay',
            plot_bgcolor='white',
            font=dict(size=20),
            title_font_size=24,
            width=1000,
            height=600,
            showlegend=False,
        )
        # Make subplot titles larger
        for annotation in fig['layout']['annotations']:
            annotation['font'] = dict(size=22)

        for col_i, (label, tech, tools) in enumerate(subplots, start=1):
            present_tools = [t for t in tools if t in df[df['technology'] == tech]['tool'].values]
            n = len(present_tools)
            bw = 0.8 / n
            tick_positions = [(j - (n - 1) / 2) * bw for j in range(n)]
            xaxis = f"xaxis{col_i}" if col_i > 1 else "xaxis"
            yaxis = f"yaxis{col_i}" if col_i > 1 else "yaxis"
            fig.update_layout(**{
                xaxis: dict(
                    tickmode='array', tickvals=tick_positions, ticktext=present_tools,
                    tickfont=dict(size=16), showgrid=False, showline=True,
                    linewidth=2, linecolor='black',
                ),
                yaxis: dict(
                    title='Maximum Memory Usage (GB)', title_font_size=20, tickfont_size=18,
                    showgrid=True, gridcolor='lightgray', showline=True,
                    linewidth=2, linecolor='black', zeroline=False,
                ),
            })

        fig.write_html(output[0])


rule plot_stratified_accuracy:
    """Accuracy against the adotto/GIAB HG002 truth versus truth allele length.

    Three rows over a shared set of length bins. The first two are the same measurement under a
    fixed and an allele-scaled tolerance; they coincide up to 300 bp, so divergence beyond that
    is the fixed 3 bp window tightening in relative terms rather than the genotypes degrading.
    The third is the call rate, which separates dropout from miscalling.

    The per-locus error distributions (`stratified_errors.tsv.gz`, with the p90/p99 columns in
    the summary table) are deliberately not plotted here: they answer a different question than
    the one this figure is for. The data is retained so the panel can be added back if needed.
    """
    input:
        "tool_comparison/stratified_accuracy.tsv"
    output:
        "tool_comparison/stratified_accuracy.html"
    run:
        import pandas as pd
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        df = pd.read_csv(input[0], sep='\t')

        # inquiSTR-spanning is computed in the tables but not plotted: it tracks inquiSTR so
        # closely that it only hides the inquiSTR line. Compare the two in the TSV instead.
        subplots = [
            ('PacBio', 'pacbio', ['inquiSTR', 'TRGT', 'LongTR']),
            ('ONT', 'ont', ['inquiSTR', 'LongTR']),
        ]

        tool_colors = {
            'inquiSTR':          '#1f77b4',
            'inquiSTR-spanning': '#9467bd',
            'TRGT':              '#ff7f0e',
            'LongTR':            '#2ca02c',
        }
        tool_labels = {'inquiSTR-spanning': 'inquiSTR (spanning reads only)'}

        # (row title, column in the table, factor to reach a percentage)
        rows = [
            ('Within 3 bp (%)', 'within_tolerance_percent', 1),
            ('Within 3 bp or 1% (%)', 'within_scaled_tolerance_percent', 1),
            ('Loci genotyped (%)', 'call_rate', 100),
        ]

        # Bin order as written by the workflow, taken from the table itself
        bins = list(dict.fromkeys(df['length_bin']))

        fig = make_subplots(
            rows=len(rows), cols=2,
            subplot_titles=[s[0] for s in subplots] + [''] * (2 * (len(rows) - 1)),
            shared_xaxes=True,
            vertical_spacing=0.055,
            horizontal_spacing=0.12,
        )

        tool_seen = set()
        for col_i, (label, tech, tools) in enumerate(subplots, start=1):
            tech_df = df[df['technology'] == tech]
            for row_i, (row_title, column, scale) in enumerate(rows, start=1):
                for tool in tools:
                    tool_df = tech_df[tech_df['tool'] == tool].set_index('length_bin')
                    if tool_df.empty:
                        continue
                    tool_df = tool_df.reindex(bins)

                    fig.add_trace(go.Scatter(
                        x=bins, y=tool_df[column] * scale,
                        mode='lines+markers',
                        name=tool_labels.get(tool, tool),
                        line=dict(color=tool_colors[tool], width=3),
                        marker=dict(color=tool_colors[tool], size=10,
                                    line=dict(width=1, color='white')),
                        legendgroup=tool,
                        showlegend=(tool not in tool_seen),
                    ), row=row_i, col=col_i)
                    tool_seen.add(tool)

        # Annotate each x tick with the number of catalog loci in that bin, which is the
        # denominator every curve is computed over and is identical across tools
        tick_text = {}
        for label, tech, tools in subplots:
            counts = (df[df['technology'] == tech]
                      .groupby('length_bin')['n_targeted'].max())
            tick_text[tech] = [f"{b}<br>n={int(counts.get(b, 0)):,}" for b in bins]

        fig.update_layout(
            title='Accuracy versus truth allele length (adotto/GIAB HG002)',
            title_x=0.5,
            plot_bgcolor='white',
            font=dict(size=18),
            title_font_size=24,
            width=1200,
            height=1100,
            margin=dict(t=120, b=220, l=110),
            legend=dict(orientation='h', yanchor='top', y=-0.14, xanchor='center', x=0.5),
        )
        for annotation in fig['layout']['annotations']:
            annotation['font'] = dict(size=22)

        # Panel letters, added after the loop above so they keep their own font. Each letter
        # labels a whole row, since a row is one measurement shown for both technologies.
        for row_i, letter in zip(range(1, len(rows) + 1), 'ABC'):
            subplot = fig.get_subplot(row=row_i, col=1)
            fig.add_annotation(
                text=f'<b>{letter}</b>',
                xref='paper', yref='paper',
                x=subplot.xaxis.domain[0] - 0.085,
                y=subplot.yaxis.domain[1] + 0.008,
                xanchor='left', yanchor='bottom',
                showarrow=False,
                font=dict(size=26, color='black'),
            )

        axis_style = dict(showline=True, linewidth=2, linecolor='black')
        for col_i, (label, tech, tools) in enumerate(subplots, start=1):
            for row_i, (row_title, column, scale) in enumerate(rows, start=1):
                is_last = row_i == len(rows)
                fig.update_xaxes(
                    tickmode='array',
                    tickvals=bins,
                    ticktext=tick_text[tech] if is_last else ['' for _ in bins],
                    title_text='Truth allele length (bp)' if is_last else None,
                    title_font_size=20, tickfont_size=14, showgrid=False,
                    row=row_i, col=col_i, **axis_style,
                )
                fig.update_yaxes(
                    title_text=row_title if col_i == 1 else None,
                    title_font_size=16, tickfont_size=16, range=[0, 105],
                    showgrid=True, gridcolor='lightgray', zeroline=False,
                    row=row_i, col=col_i, **axis_style,
                )

        # The browser's "download plot as png" button exports at the on-screen size by
        # default, which is too coarse for print. scale=3 gives 3600x3300 px.
        fig.write_html(output[0], config={
            'toImageButtonOptions': {
                'format': 'png',
                'filename': 'stratified_accuracy',
                'scale': 3,
            },
        })


rule plot_medaka_benchmark:
    """Runtime and peak memory for the three ONT callers on the full adotto catalog.

    Kept separate from the headline runtime/memory figures while medaka tandem is still being
    evaluated; merge into plot_tool_comparison_time / _memory only if it joins the main
    comparison.
    """
    input:
        "tool_comparison/medaka_results.tsv"
    output:
        "tool_comparison/medaka_benchmark_plot.html"
    run:
        import pandas as pd
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        df = pd.read_csv(input[0], sep='\t')
        df['elapsed_minutes'] = df['elapsed_seconds'] / 60

        tool_colors = {
            'inquiSTR':      '#1f77b4',
            'LongTR':        '#2ca02c',
            'medaka tandem': '#8c564b',
        }
        tools = [t for t in ['inquiSTR', 'LongTR', 'medaka tandem'] if t in df['tool'].values]
        panels = [
            ('Runtime', 'elapsed_minutes', 'Run time (minutes)', '{:.1f} min'),
            ('Peak memory', 'max_memory_gb', 'Peak memory (GB)', '{:.1f} GB'),
        ]

        fig = make_subplots(rows=1, cols=2, subplot_titles=[p[0] for p in panels],
                            horizontal_spacing=0.18)

        for col_i, (title, column, y_title, fmt) in enumerate(panels, start=1):
            for j, tool in enumerate(tools):
                tool_df = df[df['tool'] == tool]
                mean_val = tool_df[column].mean()

                fig.add_trace(go.Bar(
                    x=[j], y=[mean_val],
                    marker_color=tool_colors[tool], opacity=0.8, width=0.6,
                    text=[fmt.format(mean_val)], textposition='outside', textfont=dict(size=16),
                    showlegend=False,
                ), row=1, col=col_i)

                # individual replicates on top of the mean, as in the other benchmark figures
                fig.add_trace(go.Scatter(
                    x=[j] * len(tool_df), y=tool_df[column],
                    mode='markers',
                    marker=dict(color=tool_colors[tool], size=11, opacity=0.6,
                                line=dict(width=1, color='white')),
                    showlegend=False,
                ), row=1, col=col_i)

        fig.update_layout(
            title='ONT tool comparison: medaka tandem',
            title_x=0.5,
            plot_bgcolor='white',
            font=dict(size=18),
            title_font_size=24,
            width=1000,
            height=600,
            margin=dict(t=110, b=120, l=90),
        )
        for annotation in fig['layout']['annotations']:
            annotation['font'] = dict(size=20)

        axis_style = dict(showline=True, linewidth=2, linecolor='black')
        for col_i, (title, column, y_title, fmt) in enumerate(panels, start=1):
            fig.update_xaxes(
                tickmode='array', tickvals=list(range(len(tools))), ticktext=tools,
                tickfont_size=16, showgrid=False, row=1, col=col_i, **axis_style,
            )
            fig.update_yaxes(
                title_text=y_title, title_font_size=18, tickfont_size=16,
                showgrid=True, gridcolor='lightgray', zeroline=False,
                row=1, col=col_i, **axis_style,
            )

        fig.write_html(output[0], config={
            'toImageButtonOptions': {
                'format': 'png',
                'filename': 'medaka_benchmark',
                'scale': 3,
            },
        })


rule plot_coverage_accuracy:
    """Concordance against the adotto/GIAB HG002 truth as a function of ONT sequencing depth.

    Two panels: the agreement metrics, and the number of loci actually genotyped. The second
    matters because at low depth a locus can simply fail the read-support threshold, and
    concordance is only computed over loci that produced a genotype.
    """
    input:
        "coverage/accuracy_by_coverage.tsv"
    output:
        "coverage/accuracy_by_coverage.html"
    run:
        import pandas as pd
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        df = pd.read_csv(input[0], sep='\t').sort_values('coverage')

        metrics = [
            ('Exact match', 'exact_percent', '#1f77b4'),
            ('Within 1 bp', 'within_1bp_percent', '#ff7f0e'),
            ('Within 3 bp', 'within_3bp_percent', '#2ca02c'),
        ]

        fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.14,
                            subplot_titles=['Concordance with truth', 'Loci genotyped'])

        for label, column, color in metrics:
            fig.add_trace(go.Scatter(
                x=df['coverage'], y=df[column],
                mode='lines+markers', name=label,
                line=dict(color=color, width=3),
                marker=dict(color=color, size=10, line=dict(width=1, color='white')),
            ), row=1, col=1)

        fig.add_trace(go.Scatter(
            x=df['coverage'], y=df['loci_assessed'],
            mode='lines+markers', name='Loci assessed',
            line=dict(color='#9467bd', width=3),
            marker=dict(color='#9467bd', size=10, line=dict(width=1, color='white')),
            showlegend=False,
        ), row=1, col=2)

        fig.update_layout(
            title='Genotype concordance versus ONT coverage (HG002)',
            title_x=0.5,
            plot_bgcolor='white',
            font=dict(size=18),
            title_font_size=24,
            width=1200,
            height=600,
            margin=dict(t=110, b=140, l=90),
            legend=dict(orientation='h', yanchor='top', y=-0.18, xanchor='center', x=0.5),
        )
        for annotation in fig['layout']['annotations']:
            annotation['font'] = dict(size=20)

        # Panel letters, added after the loop above so they keep their own font
        for col_i, letter in zip((1, 2), 'AB'):
            subplot = fig.get_subplot(row=1, col=col_i)
            fig.add_annotation(
                text=f'<b>{letter}</b>',
                xref='paper', yref='paper',
                x=subplot.xaxis.domain[0] - 0.075,
                y=subplot.yaxis.domain[1] + 0.06,
                xanchor='left', yanchor='bottom',
                showarrow=False,
                font=dict(size=26, color='black'),
            )

        axis_style = dict(showline=True, linewidth=2, linecolor='black')
        for col_i, y_title in enumerate(['Genotypes matching truth (%)', 'Loci assessed'], start=1):
            fig.update_xaxes(
                title_text='Coverage (x)', title_font_size=20, tickfont_size=16,
                tickmode='array', tickvals=list(df['coverage']),
                showgrid=False, row=1, col=col_i, **axis_style,
            )
            fig.update_yaxes(
                title_text=y_title, title_font_size=18, tickfont_size=16,
                showgrid=True, gridcolor='lightgray', zeroline=False,
                row=1, col=col_i, **axis_style,
            )
        fig.update_yaxes(range=[0, 105], row=1, col=1)

        fig.write_html(output[0], config={
            'toImageButtonOptions': {
                'format': 'png',
                'filename': 'accuracy_by_coverage',
                'scale': 3,
            },
        })


rule plot_straglr_chr21_runtime:
    input:
        "tool_comparison/straglr_chr21_results.tsv"
    output:
        "tool_comparison/straglr_chr21_runtime_plot.html"
    run:
        import pandas as pd
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        df = pd.read_csv(input[0], sep='\t')
        df['tool'] = df['tool'].replace({'STRaglr': 'Straglr'})
        df['elapsed_minutes'] = df['elapsed_seconds'] / 60

        subplots = [
            ('PacBio', 'pacbio'),
            ('ONT', 'ont'),
        ]

        tool_order = ['inquiSTR', 'Straglr']
        tool_colors = {
            'inquiSTR': '#1f77b4',
            'Straglr': '#d62728',
        }

        fig = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=[s[0] for s in subplots],
            horizontal_spacing=0.12,
        )

        for col_i, (label, tech) in enumerate(subplots, start=1):
            tech_df = df[df['technology'] == tech]
            present_tools = [t for t in tool_order if t in tech_df['tool'].values]
            n = len(present_tools)
            if n == 0:
                continue

            group_width = 0.8
            bw = group_width / n

            for j, tool in enumerate(present_tools):
                offset = (j - (n - 1) / 2) * bw
                x = offset
                color = tool_colors[tool]
                tool_df = tech_df[tech_df['tool'] == tool]
                mean_val = tool_df['elapsed_minutes'].mean()

                fig.add_trace(go.Bar(
                    x=[x], y=[mean_val],
                    name=tool,
                    marker_color=color,
                    opacity=0.85,
                    width=bw * 0.9,
                    text=[f"{mean_val:.2f} min"],
                    textposition='outside',
                    textfont=dict(size=15),
                    legendgroup=tool,
                    showlegend=False,
                ), row=1, col=col_i)

                fig.add_trace(go.Scatter(
                    x=[x] * len(tool_df),
                    y=tool_df['elapsed_minutes'],
                    mode='markers',
                    marker=dict(
                        color=color,
                        size=11,
                        opacity=0.65,
                        line=dict(width=1, color='white'),
                    ),
                    legendgroup=tool,
                    showlegend=False,
                    hovertemplate=f"Tool: {tool}<br>Tech: {tech}<br>Replicate runtime: %{{y:.2f}} min<extra></extra>",
                ), row=1, col=col_i)

            tick_positions = [(j - (n - 1) / 2) * bw for j in range(n)]
            xaxis = f"xaxis{col_i}" if col_i > 1 else "xaxis"
            yaxis = f"yaxis{col_i}" if col_i > 1 else "yaxis"
            fig.update_layout(**{
                xaxis: dict(
                    tickmode='array',
                    tickvals=tick_positions,
                    ticktext=present_tools,
                    tickfont=dict(size=15),
                    showgrid=False,
                    showline=True,
                    linewidth=2,
                    linecolor='black',
                ),
                yaxis: dict(
                    title='Run time (minutes)',
                    title_font_size=19,
                    tickfont_size=16,
                    showgrid=True,
                    gridcolor='lightgray',
                    showline=True,
                    linewidth=2,
                    linecolor='black',
                    zeroline=False,
                ),
            })

        fig.update_layout(
            title='Run time: inquiSTR vs Straglr',
            title_x=0.5,
            barmode='overlay',
            plot_bgcolor='white',
            font=dict(size=18),
            title_font_size=24,
            width=600,
            height=620,
            showlegend=False,
        )

        for annotation in fig['layout']['annotations']:
            annotation['font'] = dict(size=21)

        fig.write_html(output[0])


rule plot_straglr_chr21_memory:
    input:
        "tool_comparison/straglr_chr21_results.tsv"
    output:
        "tool_comparison/straglr_chr21_memory_plot.html"
    run:
        import pandas as pd
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        df = pd.read_csv(input[0], sep='\t')
        df['tool'] = df['tool'].replace({'STRaglr': 'Straglr'})

        subplots = [
            ('PacBio', 'pacbio'),
            ('ONT', 'ont'),
        ]

        tool_order = ['inquiSTR', 'Straglr']
        tool_colors = {
            'inquiSTR': '#1f77b4',
            'Straglr': '#d62728',
        }

        fig = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=[s[0] for s in subplots],
            horizontal_spacing=0.12,
        )

        for col_i, (label, tech) in enumerate(subplots, start=1):
            tech_df = df[df['technology'] == tech]
            present_tools = [t for t in tool_order if t in tech_df['tool'].values]
            n = len(present_tools)
            if n == 0:
                continue

            group_width = 0.8
            bw = group_width / n

            for j, tool in enumerate(present_tools):
                offset = (j - (n - 1) / 2) * bw
                x = offset
                color = tool_colors[tool]
                tool_df = tech_df[tech_df['tool'] == tool]
                mean_val = tool_df['max_memory_gb'].mean()

                fig.add_trace(go.Bar(
                    x=[x], y=[mean_val],
                    name=tool,
                    marker_color=color,
                    opacity=0.85,
                    width=bw * 0.9,
                    text=[f"{mean_val:.2f} GB"],
                    textposition='outside',
                    textfont=dict(size=15),
                    legendgroup=tool,
                    showlegend=False,
                ), row=1, col=col_i)

                fig.add_trace(go.Scatter(
                    x=[x] * len(tool_df),
                    y=tool_df['max_memory_gb'],
                    mode='markers',
                    marker=dict(
                        color=color,
                        size=11,
                        opacity=0.65,
                        line=dict(width=1, color='white'),
                    ),
                    legendgroup=tool,
                    showlegend=False,
                    hovertemplate=f"Tool: {tool}<br>Tech: {tech}<br>Replicate memory: %{{y:.2f}} GB<extra></extra>",
                ), row=1, col=col_i)

            tick_positions = [(j - (n - 1) / 2) * bw for j in range(n)]
            xaxis = f"xaxis{col_i}" if col_i > 1 else "xaxis"
            yaxis = f"yaxis{col_i}" if col_i > 1 else "yaxis"
            fig.update_layout(**{
                xaxis: dict(
                    tickmode='array',
                    tickvals=tick_positions,
                    ticktext=present_tools,
                    tickfont=dict(size=15),
                    showgrid=False,
                    showline=True,
                    linewidth=2,
                    linecolor='black',
                ),
                yaxis: dict(
                    title='Maximum Memory Usage (GB)',
                    title_font_size=19,
                    tickfont_size=16,
                    showgrid=True,
                    gridcolor='lightgray',
                    showline=True,
                    linewidth=2,
                    linecolor='black',
                    zeroline=False,
                ),
            })

        fig.update_layout(
            title='Memory usage: inquiSTR vs Straglr',
            title_x=0.5,
            barmode='overlay',
            plot_bgcolor='white',
            font=dict(size=18),
            title_font_size=24,
            width=1050,
            height=620,
            showlegend=False,
        )

        for annotation in fig['layout']['annotations']:
            annotation['font'] = dict(size=21)

        fig.write_html(output[0])


rule plot_polymorphic_pca:
    input:
        scores = "polymorphic/repeats_pca_scores.tsv",
        cohort = file_path
    output:
        "polymorphic/repeats_pca_colored.html"
    run:
        import pandas as pd
        import plotly.express as px

        scores = pd.read_csv(input.scores, sep='\t')
        cohort = pd.read_csv(input.cohort, sep='\t', usecols=['sample', 'Superpopulation code'])
        df = scores.merge(cohort, on='sample', how='left')

        superpop_colors = {
            'AFR': '#E41A1C',
            'AMR': '#FF7F00',
            'EAS': '#4DAF4A',
            'EUR': '#377EB8',
            'SAS': '#984EA3',
        }

        fig = px.scatter(
            df,
            x='PC1',
            y='PC2',
            color='Superpopulation code',
            color_discrete_map=superpop_colors,
            hover_data=['sample'],
            title='Polymorphic STR PCA \u2014 PC1 vs PC2',
        )
        fig.update_traces(marker_size=12)
        fig.update_layout(
            plot_bgcolor='white',
            font=dict(size=22),
            title_font_size=26,
            legend_title_text=None,
            legend=dict(
                yanchor='top', y=0.99,
                xanchor='right', x=0.99,
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='grey', borderwidth=0.5,
                font=dict(size=18),
                title_font_size=22,
            ),
            width=600,
            height=600,
        )
        fig.update_xaxes(showgrid=False, zeroline=False, showline=True, linewidth=2, linecolor='black', title_font_size=24, tickfont_size=20)
        fig.update_yaxes(showgrid=False, zeroline=False, showline=True, linewidth=2, linecolor='black', title_font_size=24, tickfont_size=20)
        fig.write_html(output[0])


rule plot_puretarget_heatmap:
    input:
        combined = "puretarget-calls/combined.tsv"
    output:
        html = "puretarget-calls/heatmap.html",
        png = "puretarget-calls/heatmap.png"
    run:
        import pandas as pd
        import plotly.graph_objects as go

        df = pd.read_csv(input.combined, sep='\t', comment='#')

        # Build a locus label from chrom/begin/end only (explicit no info field usage)
        def locus_label(row):
            return f"{row['chromosome']}:{row['begin']}-{row['end']}"

        df['locus'] = df.apply(locus_label, axis=1)
        # subset to loci and gene labels
        loci_of_interest = {
            "chr4:3074877-3074933": "HTT",
            "chr6:16327634-16327724": "ATXN1",
            "chr19:45770205-45770266": "DMPK",
            "chr4:39348425-39348483": "RFC1",
            "chr9:69037287-69037304": "FXN",
        }

        df = df[df['locus'].isin(loci_of_interest.keys())].copy()

        # use gene names as x-axis labels in the heatmap
        locus_to_gene = loci_of_interest

        # Detect sample columns: everything after the fixed metadata columns
        meta_cols = ['chromosome', 'begin', 'end', 'locus']
        sample_cols = [c for c in df.columns if c not in meta_cols]

        # For each sample two haplotype columns are expected (sampleX_H1, sampleX_H2)
        # Group them to get the max (longest allele) per sample
        samples = sorted(set(c.rsplit('_H', 1)[0] for c in sample_cols if '_H' in c))

        # Clean sample names by removing the prefix if present
        prefix = 'mapped_m84039_250829_182713_s3.hifi_reads.'
        cleaned_samples = [s.replace(prefix, '') for s in samples]

        ref_len = (df['end'] - df['begin']).astype(float)
        # removing duplicates from the visualization
        duplicates_to_skip = [
            "bc2010",
            "bc2011",
            "bc2012",
            "bc2013",
            "bc2014",
            "bc2015",
            "bc2018",
            "bc2019",
            "bc2022",
            "bc2023",
            "bc2033",
            "bc2034",
            "bc2035",
            "bc2037",
            "bc2038",
            "bc2039",
            "bc2042",
            "bc2043",
            "bc2046",
            "bc2047"
        ]
        # removing a set of unexpanded random donors that are not informative for the visualization
        unexpanded_to_skip = [
            "bc2033",
            "bc2034",
            "bc2035",
            "bc2036",
            "bc2037",
            "bc2038",
            "bc2039",
            "bc2040",
            "bc2041",
            "bc2042",
            "bc2043",
            "bc2044",
            "bc2045",
            "bc2046",
            "bc2047",
            "bc2048",
        ]

        barcode_to_sample = {
            "bc2009": "donor4888",
            "bc2016": "NA24385",
            "bc2017": "donor8375",
            "bc2020": "donor8359",
            "bc2021": "donor8360",
            "bc2024": "donor8361",
            "bc2025": "HM13509",
            "bc2026": "HM06926",
            "bc2027": "NA13537",
            "bc2028": "NA13536",
            "bc2029": "NA03697",
            "bc2030": "NA13509",
            "bc2031": "NA13515",
            "bc2032": "HG01175",
            "bc2049": "HG04228",
            "bc2050": "NA20752",
            "bc2051": "NA15848",
            "bc2052": "NA16212",
            "bc2053": "NA16202",
            "bc2054": "NA16213",
            "bc2055": "NA16237",
            "bc2056": "NA23629",
        }

        divisor_map = {gene: (5 if gene == 'RFC1' else 3) for gene in locus_to_gene.values()}

        records = []
        for sample, cleaned_sample in zip(samples, cleaned_samples):
            if cleaned_sample in duplicates_to_skip + unexpanded_to_skip:
                continue
            h1 = f"{sample}_H1"
            h2 = f"{sample}_H2"
            available = [c for c in [h1, h2] if c in df.columns]
            max_len = df[available].apply(pd.to_numeric, errors='coerce').max(axis=1)
            # add reference repeat length to allele length
            max_len_adjusted = max_len + ref_len

            # apply scaling by locus definition
            scaled_per_locus = []
            for i, locus in enumerate(df['locus']):
                gene_name = locus_to_gene.get(locus, locus)
                divisor = divisor_map.get(gene_name, 3)
                val = max_len_adjusted.iloc[i]
                scaled_per_locus.append(round(val / divisor) if not pd.isna(val) else np.nan)
            max_len_adjusted = pd.Series(scaled_per_locus, index=df.index)

            records.append(max_len_adjusted.rename(cleaned_sample))

        # Matrix: rows = loci, columns = samples
        matrix = pd.concat(records, axis=1)
        matrix.index = df['locus'].values

        # Label rows with gene names where available
        matrix.index = [locus_to_gene.get(l, l) for l in matrix.index]

        # Rename columns from barcodes to sample names
        matrix.columns = [barcode_to_sample.get(c, c) for c in matrix.columns]

        # Scale each row (locus) to [0, 1] for per-repeat colouring
        scaled = matrix.copy().astype(float)
        for row in scaled.index:
            row_min = scaled.loc[row].min()
            row_max = scaled.loc[row].max()
            if row_max > row_min:
                scaled.loc[row] = (scaled.loc[row] - row_min) / (row_max - row_min)
            else:
                scaled.loc[row] = 0.0

        # Build hover text with actual lengths
        hover = matrix.map(lambda v: f"{v:.0f} bp" if pd.notna(v) else "NA")

        # Render actual allele length labels in each cell
        cell_labels = matrix.applymap(lambda v: f"{int(v):d}" if pd.notna(v) else "")

        fig = go.Figure(go.Heatmap(
            z=scaled.values,
            x=scaled.columns.tolist(),
            y=scaled.index.tolist(),
            text=cell_labels.values,
            texttemplate="%{text}",
            textfont=dict(size=13, color='black'),
            customdata=hover.values,
            hovertemplate="Sample: %{y}<br>Locus: %{x}<br>Length: %{customdata}<extra></extra>",
            colorscale=[[0, 'white'], [1, '#d62728']],
            showscale=True,
            colorbar=dict(
                title=None,
                orientation='h',
                x=0.5,
                y=-0.42,
                xanchor='center',
                yanchor='top',
                tickvals=[0, 0.5, 1],
                ticktext=['short', 'mid', 'long'],
                tickfont=dict(size=17),
                len=0.8,
                thickness=18,
            ),
            zmin=0, zmax=1,
        ))

        fig.update_layout(
            title='Pathogenic TR lengths',
            title_x=0.5,
            title_font_size=26,
            plot_bgcolor='white',
            font=dict(size=17),
            xaxis=dict(
                title='',
                title_font_size=20,
                tickfont=dict(size=18),
                tickangle=45,
                showgrid=False,
                automargin=True,
            ),
            yaxis=dict(
                title='',
                title_font_size=20,
                tickfont=dict(size=18),
                showgrid=False,
                autorange='reversed',
                automargin=True,
            ),
            margin=dict(l=0, b=210, t=60, r=0),
            width=800,
            height=600,
        )

        # Draw boxes around the three sample groups
        samples_list = scaled.columns.tolist()
        n_loci = len(scaled.index)
        y0, y1 = -0.5, n_loci - 0.5

        # Group 1: non-carriers (start to donor8361)
        # Group 2: carriers (HM13509 to NA20752)
        # Group 3: FXN carriers (NA15848 to end)
        if 'donor8361' in samples_list:
            idx = samples_list.index('donor8361')
            fig.add_shape(type='rect', x0=-0.5, x1=idx + 0.5, y0=y0, y1=y1,
                          line=dict(color='black', width=2))
        if 'HM13509' in samples_list and 'NA20752' in samples_list:
            idx_start = samples_list.index('HM13509')
            idx_end = samples_list.index('NA20752')
            fig.add_shape(type='rect', x0=idx_start - 0.5, x1=idx_end + 0.5, y0=y0, y1=y1,
                          line=dict(color='black', width=2))
        if 'NA15848' in samples_list:
            idx_start = samples_list.index('NA15848')
            fig.add_shape(type='rect', x0=idx_start - 0.5, x1=len(samples_list) - 0.5, y0=y0, y1=y1,
                          line=dict(color='black', width=2))

        fig.write_html(output.html)
        # High-resolution raster for publication: scale=4 renders at ~4x the on-screen
        # pixel size (800x400 -> 3200x1600). Requires the kaleido package in the env.
        fig.write_image(output.png, scale=4)
