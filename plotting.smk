"""
Plotting rules for inquiSTR paper analysis
"""

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

        # Define ordering
        tool_order = ['inquiSTR', 'TRGT', 'inquiSTR+TRGT', 'LongTR', 'inquiSTR+LongTR']
        tech_order = ['pacbio', 'ont']
        tech_labels = {'pacbio': 'PacBio', 'ont': 'ONT'}

        # One color per tool
        tool_colors = {
            'inquiSTR':       '#1f77b4',
            'TRGT':           '#ff7f0e',
            'inquiSTR+TRGT':  '#d62728',
            'LongTR':         '#2ca02c',
            'inquiSTR+LongTR':'#9467bd',
        }

        group_width = 0.8
        tool_seen = set()
        fig = go.Figure()

        for tech_i, tech in enumerate(tech_order):
            tech_df = df[df['technology'] == tech]
            present_tools = [t for t in tool_order if t in tech_df['tool'].values]
            n = len(present_tools)
            bw = group_width / n

            for j, tool in enumerate(present_tools):
                offset = (j - (n - 1) / 2) * bw
                x = tech_i + offset
                color = tool_colors[tool]
                tool_df = tech_df[tech_df['tool'] == tool]
                mean_val = tool_df['elapsed_minutes'].mean()

                fig.add_trace(go.Bar(
                    x=[x],
                    y=[mean_val],
                    name=tool,
                    marker_color=color,
                    opacity=0.8,
                    width=bw * 0.9,
                    text=[f"{mean_val:.2f} min"],
                    textposition='outside',
                    legendgroup=tool,
                    showlegend=(tool not in tool_seen),
                ))
                tool_seen.add(tool)

                fig.add_trace(go.Scatter(
                    x=[x] * len(tool_df),
                    y=tool_df['elapsed_minutes'],
                    mode='markers',
                    marker=dict(color=color, size=7, opacity=0.6,
                                line=dict(width=1, color='white')),
                    legendgroup=tool,
                    showlegend=False,
                ))

        fig.update_layout(
            title='Tool Comparison: Runtime (mean with individual replicates)',
            xaxis_title='Technology',
            yaxis_title='Elapsed Time (minutes)',
            barmode='overlay',
            plot_bgcolor='white',
            xaxis=dict(
                tickmode='array',
                tickvals=list(range(len(tech_order))),
                ticktext=[tech_labels[t] for t in tech_order],
            ),
            legend=dict(yanchor='top', y=0.99, xanchor='right', x=0.99),
        )
        fig.update_xaxes(showgrid=False, showline=True, linewidth=2, linecolor='black')
        fig.update_yaxes(showgrid=True, gridcolor='lightgray', showline=True,
                         linewidth=2, linecolor='black')

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

        # Define ordering
        tool_order = ['inquiSTR', 'TRGT', 'inquiSTR+TRGT', 'LongTR', 'inquiSTR+LongTR']
        tech_order = ['pacbio', 'ont']
        tech_labels = {'pacbio': 'PacBio', 'ont': 'ONT'}

        # One color per tool
        tool_colors = {
            'inquiSTR':       '#1f77b4',
            'TRGT':           '#ff7f0e',
            'inquiSTR+TRGT':  '#d62728',
            'LongTR':         '#2ca02c',
            'inquiSTR+LongTR':'#9467bd',
        }

        group_width = 0.8
        tool_seen = set()
        fig = go.Figure()

        for tech_i, tech in enumerate(tech_order):
            tech_df = df[df['technology'] == tech]
            present_tools = [t for t in tool_order if t in tech_df['tool'].values]
            n = len(present_tools)
            bw = group_width / n

            for j, tool in enumerate(present_tools):
                offset = (j - (n - 1) / 2) * bw
                x = tech_i + offset
                color = tool_colors[tool]
                tool_df = tech_df[tech_df['tool'] == tool]
                mean_val = tool_df['max_memory_gb'].mean()

                fig.add_trace(go.Bar(
                    x=[x],
                    y=[mean_val],
                    name=tool,
                    marker_color=color,
                    opacity=0.8,
                    width=bw * 0.9,
                    text=[f"{mean_val:.2f} GB"],
                    textposition='outside',
                    legendgroup=tool,
                    showlegend=(tool not in tool_seen),
                ))
                tool_seen.add(tool)

                fig.add_trace(go.Scatter(
                    x=[x] * len(tool_df),
                    y=tool_df['max_memory_gb'],
                    mode='markers',
                    marker=dict(color=color, size=7, opacity=0.6,
                                line=dict(width=1, color='white')),
                    legendgroup=tool,
                    showlegend=False,
                ))

        fig.update_layout(
            title='Tool Comparison: Memory Usage (mean with individual replicates)',
            xaxis_title='Technology',
            yaxis_title='Maximum Memory Usage (GB)',
            barmode='overlay',
            plot_bgcolor='white',
            xaxis=dict(
                tickmode='array',
                tickvals=list(range(len(tech_order))),
                ticktext=[tech_labels[t] for t in tech_order],
            ),
            legend=dict(yanchor='top', y=0.99, xanchor='right', x=0.99),
        )
        fig.update_xaxes(showgrid=False, showline=True, linewidth=2, linecolor='black')
        fig.update_yaxes(showgrid=True, gridcolor='lightgray', showline=True,
                         linewidth=2, linecolor='black')

        fig.write_html(output[0])
