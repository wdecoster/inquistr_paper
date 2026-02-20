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

        # Define a fixed x-axis order for tools
        tool_order = ['inquiSTR', 'TRGT', 'LongTR', 'inquiSTR+TRGT', 'inquiSTR+LongTR']

        # Color and display name per technology
        tech_config = {
            'pacbio': {'color': 'purple', 'label': 'PacBio'},
            'ont':    {'color': 'steelblue', 'label': 'ONT'},
        }

        fig = go.Figure()

        for tech, cfg in tech_config.items():
            tech_df = df[df['technology'] == tech]
            mean_df = tech_df.groupby('tool')['elapsed_minutes'].mean().reset_index()

            # Bars for mean values
            fig.add_trace(go.Bar(
                x=mean_df['tool'],
                y=mean_df['elapsed_minutes'],
                name=cfg['label'],
                marker_color=cfg['color'],
                opacity=0.7,
                text=[f"{v:.2f} min" for v in mean_df['elapsed_minutes']],
                textposition='outside',
            ))

            # Scatter points for individual replicates
            for tool in tech_df['tool'].unique():
                tool_data = tech_df[tech_df['tool'] == tool]
                fig.add_trace(go.Scatter(
                    x=[tool] * len(tool_data),
                    y=tool_data['elapsed_minutes'],
                    mode='markers',
                    name=f'{cfg["label"]} replicates',
                    marker=dict(color=cfg['color'], size=8, opacity=0.6,
                                line=dict(width=1, color='white')),
                    showlegend=False,
                ))

        fig.update_layout(
            title='Tool Comparison: Runtime (mean with individual replicates)',
            xaxis_title='Tool',
            yaxis_title='Elapsed Time (minutes)',
            barmode='group',
            plot_bgcolor='white',
            xaxis=dict(
                categoryorder='array',
                categoryarray=tool_order,
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

        # Define a fixed x-axis order for tools
        tool_order = ['inquiSTR', 'TRGT', 'LongTR', 'inquiSTR+TRGT', 'inquiSTR+LongTR']

        # Color and display name per technology
        tech_config = {
            'pacbio': {'color': 'purple', 'label': 'PacBio'},
            'ont':    {'color': 'steelblue', 'label': 'ONT'},
        }

        fig = go.Figure()

        for tech, cfg in tech_config.items():
            tech_df = df[df['technology'] == tech]
            mean_df = tech_df.groupby('tool')['max_memory_gb'].mean().reset_index()

            # Bars for mean values
            fig.add_trace(go.Bar(
                x=mean_df['tool'],
                y=mean_df['max_memory_gb'],
                name=cfg['label'],
                marker_color=cfg['color'],
                opacity=0.7,
                text=[f"{v:.2f} GB" for v in mean_df['max_memory_gb']],
                textposition='outside',
            ))

            # Scatter points for individual replicates
            for tool in tech_df['tool'].unique():
                tool_data = tech_df[tech_df['tool'] == tool]
                fig.add_trace(go.Scatter(
                    x=[tool] * len(tool_data),
                    y=tool_data['max_memory_gb'],
                    mode='markers',
                    name=f'{cfg["label"]} replicates',
                    marker=dict(color=cfg['color'], size=8, opacity=0.6,
                                line=dict(width=1, color='white')),
                    showlegend=False,
                ))

        fig.update_layout(
            title='Tool Comparison: Memory Usage (mean with individual replicates)',
            xaxis_title='Tool',
            yaxis_title='Maximum Memory Usage (GB)',
            barmode='group',
            plot_bgcolor='white',
            xaxis=dict(
                categoryorder='array',
                categoryarray=tool_order,
            ),
            legend=dict(yanchor='top', y=0.99, xanchor='right', x=0.99),
        )
        fig.update_xaxes(showgrid=False, showline=True, linewidth=2, linecolor='black')
        fig.update_yaxes(showgrid=True, gridcolor='lightgray', showline=True,
                         linewidth=2, linecolor='black')

        fig.write_html(output[0])
