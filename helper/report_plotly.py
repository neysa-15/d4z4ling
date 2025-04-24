import pandas as pd
import numpy as np

import plotly
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
# import scipy

import argparse
# import sys

# Define a fixed marker size
MARKER_SIZE = 8

# Distribution plot
def get_distribution_plot(df, threshold=33000):

    # Convert column values to kilobases (divide by 1000)
    df["ReadLength"] = df["ReadLength"] # / 1000
    threshold_kb = threshold #/ 1000

    x_min, x_max = df["ReadLength"].min(), df["ReadLength"].max()
    
    # Create tick marks for better x-axis readability
    tick_vals = np.linspace(x_min, x_max, num=10)  # Adjust number of ticks as needed

    # Create the distribution plot
    fig = ff.create_distplot(
        [df["ReadLength"]],  # Data for distribution
        ["ReadLength"],  # Label
        show_hist=False,  
        show_rug=False,
        colors=["blue"]  # Apply color
    )

    # Fix x-axis to show correct ReadLength values
    fig.update_xaxes(
        title_text=f"ReadLength (kb)",  
        range=[x_min, x_max],  # Set x-axis range to match ReadLength values
        tickvals=tick_vals,  # Set tick values to align with ReadLength
        tickformat=".1f"  # Format to 1 decimal place
    )

    fig.update_yaxes(
        title_text="density"
    )

    # Get the maximum y-value of the density curve
    y_max = max(fig.data[0].y)  # Extract the highest density value
    y_limit = y_max + 3e-6  # Add 5 micro to extend vertical lines slightly above

    # Ensure the reference lines show up by adding as traces
    fig.add_trace(go.Scatter(
        x=[33000, 33000], y=[0, y_limit],
        mode="lines", line=dict(color="black", dash="dash"),
        name="Readlength 33kb threshold",
        showlegend=False
    ))

    # fig.add_vline(x=33, line_width=1.5, line_dash="dash", line_color="green", opacity=0.5)

    return fig

# ---------------------------------------------
# 1. D4Z4 Repeat Length Dot Plot with Read Classification
# ---------------------------------------------

def get_d4z4length_classification_dotplot_plotly(df, unique_labels, color_map):
    df = df.loc[df['GenomeCoords'].str.contains("chr4", na=False)].copy()
    df['Haplotype'] = df['Haplotype'].fillna('Unclassified')

    haplotype_categories = {
        "4qA": ["Complete 4qA", "Partial distal 4qA"],
        "4qB": ["Complete 4qB", "Partial distal 4qB"],
        "Unclassified": ["Partial proximal Unclassified", "Partial internal Unclassified", "Partial distal Unclassified"]
    }

    fig = make_subplots(rows=1, cols=3, shared_yaxes=True, subplot_titles=["4qA", "4qB", "Unclassified"])
    #color_map = {label: px.colors.qualitative.Set1[i] for i, label in enumerate(df["ReadLabel"].unique())}

    x_positions = {label: i for i, label in enumerate(sum(haplotype_categories.values(), []))}

    for i, haplotype in enumerate(["4qA", "4qB", "Unclassified"]):
        subset = df.loc[df["ReadLabel"].isin(haplotype_categories[haplotype])].copy()
        subset["X_Jittered"] = subset["ReadLabel"].map(x_positions) + np.random.uniform(-0.2, 0.2, size=len(subset))

        for read_label in haplotype_categories[haplotype]:
            filtered_df = subset[subset["ReadLabel"] == read_label]
            fig.add_trace(go.Scatter(
                x=filtered_df["X_Jittered"],
                y=filtered_df["MappedEstimatedCopies"],
                mode="markers",
                marker=dict(color=color_map.get(read_label, "gray"), size=MARKER_SIZE, opacity=0.5),
                name=read_label,
                legendgroup=read_label,  # Group legend entries
                showlegend=False,  # Hide duplicate entries
                customdata=filtered_df["ReadID"],
                hovertemplate="ReadID: %{customdata}<br>Copies: %{y}<extra></extra>"
            ), row=1, col=i+1)

    fig.update_yaxes(title_text="D4Z4 Copies")

    return fig

# ---------------------------------------------
# 2. Methylation Scatter Plot
# ---------------------------------------------

def get_methylation_plot_plotly(df, haplotype, unique_labels, color_map):
    df = df.loc[df['Haplotype'].str.contains(haplotype, na=False)].copy()

    # filter out copies under 3
    df = df[df["MappedEstimatedCopies"] >= 3]
 
    # Determine max number of copies and add 5
    max_copies = df["MappedEstimatedCopies"].max() + 5

    fig = px.scatter(
        df,
        x="pLAM_Methylation_Percentage",
        y="MappedEstimatedCopies",
        color="ReadLabel",
        title="D4Z4 Repeat Length vs Distal-Copy CpG Methylation",
        labels={"pLAM_Methylation_Percentage": "Distal D4Z4 Methylation (%)", "MappedEstimatedCopies": "D4Z4 Copies"},
        color_discrete_map=color_map,
        category_orders={"ReadLabel": unique_labels},
        custom_data=["ReadID"]
    )

    # Ensure the reference lines show up by adding as traces
    fig.add_trace(go.Scatter(
        x=[50, 50], y=[0, max_copies],
        mode="lines", line=dict(color="black", dash="dash"),
        name="50% Methylation"
    ))

    fig.add_trace(go.Scatter(
        x=[0, 100], y=[10, 10],
        mode="lines", line=dict(color="black", dash="dash"),
        name="10 Copies Threshold"
    ))

    fig.update_traces(marker=dict(size=MARKER_SIZE, opacity=0.5))

    fig.update_traces(
        hovertemplate="ReadID: %{customdata[0]}<br>Methylation: %{x}%<br>Copies: %{y}<extra></extra>"
    )

    return fig

# ---------------------------------------------
# 3. Restriction Site Profiles
# ---------------------------------------------

def get_4qA_restriction_site_plot_plotly(df, unique_labels, color_map):
    df = df.loc[df['Haplotype'].str.contains("4qA", na=False)].copy()

    # Find the max value for axis alignment
    max_range = max(df["XapI_Sensitive_Repeats"].max(), df["BlnI_Sensitive_Repeats"].max())


    fig = px.scatter(
        df, x="XapI_Sensitive_Repeats", y="BlnI_Sensitive_Repeats", color="ReadLabel",
        title="Restriction Site Profiles",
        labels={"XapI_Sensitive_Repeats": "XapI D4Z4 Copies", "BlnI_Sensitive_Repeats": "BlnI D4Z4 Copies"},
        color_discrete_map=color_map,
        category_orders={"ReadLabel": unique_labels}
    )

    # Add diagonal reference line
    fig.add_trace(go.Scatter(
        x=[0, max_range], y=[0, max_range],  # Ensure diagonal covers full range
        mode="lines", line=dict(dash="dash", color="black"),
        name="Diagonal Reference"
    ))

    # Set both axes to the same range
    fig.update_xaxes(range=[0, max_range])
    fig.update_yaxes(range=[0, max_range])
    
    fig.update_traces(marker=dict(size=MARKER_SIZE,  opacity=0.5))

    return fig

# ---------------------------------------------
# 4. Combine Plots in Grid Layout
# ---------------------------------------------

def combine_plots(df, dir_path, sample, unique_labels, color_map):
    fig = make_subplots(
        rows=5, cols=2,
        specs=[
            [{}, {}],
            [{"colspan": 2}, None],  # Row 1
            [{}, {}],
            [{}, {}],
            [{}, {}] 
        ],
        subplot_titles=["Read length distribution plot", "", "D4Z4 Repeat Length Dot Plot", "4qA Methylation vs D4Z4 Copies", "4qB Methylation vs D4Z4 Copies", "10qA Methylation vs D4Z4 Copies", "10qB Methylation vs D4Z4 Copies", "4qA Restriction Sites"],
        row_heights=[0.4, 0.4, 0.6, 0.6, 0.6],
        vertical_spacing=0.05
    )
   
    traces = []

    # Row 0 (Distribution Plot)
    for trace in get_distribution_plot(df).data:
        trace.showlegend = True if trace.name not in [t.name for t in traces] else False
        traces.append(trace)
        fig.add_trace(trace, row=1, col=1)

    # Row 1 (Classification Dotplot)
    for trace in get_d4z4length_classification_dotplot_plotly(df, unique_labels, color_map).data:
        trace.showlegend = True if trace.name not in [t.name for t in traces] else False  # Show legend only once
        traces.append(trace)
        fig.add_trace(trace, row=2, col=1)

    # Update x-axis labels in `combine_plots`
    new_labels = {
        "Complete 4qA": "Complete 4qA",
        "Partial distal 4qA": "Partial distal 4qA",
        "Complete 4qB": "Complete 4qB",
        "Partial distal 4qB": "Partial distal 4qB",
        "Partial proximal Unclassified": "Partial proximal Unclassified",
        "Partial internal Unclassified": "Partial internal Unclassified",
        "Partial distal Unclassified": "Partial distal Unclassified"
    }
    ordered_labels = list(new_labels.keys())
    x_positions = {label: i for i, label in enumerate(ordered_labels)}

    fig.update_xaxes(
        row=2, col=1,  # Apply to the correct subplot
        tickvals=list(x_positions.values()),
        ticktext=[new_labels[label] for label in ordered_labels]
    )

    # Row 2 (Methylation Scatter)
    for trace in get_methylation_plot_plotly(df, "4qA", unique_labels, color_map).data:
        trace.showlegend = True if trace.name not in [t.name for t in traces] else False
        traces.append(trace)
        fig.add_trace(trace, row=3, col=1)

    # Row 2 (Methylation Scatter)
    for trace in get_methylation_plot_plotly(df, "4qB", unique_labels, color_map).data:
        trace.showlegend = True if trace.name not in [t.name for t in traces] else False
        traces.append(trace)
        fig.add_trace(trace, row=3, col=2)

    # Row 2 (Methylation Scatter)
        for trace in get_methylation_plot_plotly(df, "10qA", unique_labels, color_map).data:
            trace.showlegend = True if trace.name not in [t.name for t in traces] else False
            traces.append(trace)
            fig.add_trace(trace, row=4, col=1)

    # Row 2 (Methylation Scatter)
    for trace in get_methylation_plot_plotly(df, "10qB", unique_labels, color_map).data:
        trace.showlegend = True if trace.name not in [t.name for t in traces] else False
        traces.append(trace)
        fig.add_trace(trace, row=4, col=2)

    # Row 2 (Restriction Sites)
    for trace in get_4qA_restriction_site_plot_plotly(df, unique_labels, color_map).data:
        trace.showlegend = True if trace.name not in [t.name for t in traces] else False
        traces.append(trace)
        fig.add_trace(trace, row=5, col=1)

    # **Add X-axis and Y-axis labels**
    fig.update_xaxes(title_text="Read length", row=1, col=1)
    fig.update_yaxes(title_text="Density", row=1, col=1)

    fig.update_yaxes(title_text="D4Z4 Repeat length", row=2, col=1)

    fig.update_yaxes(title_text="D4Z4 Copies", row=3, col=1)  # Methylation Scatter X-axis
    fig.update_xaxes(title_text="Methylation (%)", row=3, col=1)  # Methylation Scatter Y-axis

    fig.update_yaxes(title_text="D4Z4 Copies", row=3, col=2)  # Methylation Scatter X-axis
    fig.update_xaxes(title_text="Methylation (%)", row=3, col=2)  # Methylation Scatter Y-axis

    fig.update_yaxes(title_text="D4Z4 Copies", row=4, col=1)  # Methylation Scatter X-axis
    fig.update_xaxes(title_text="Methylation (%)", row=4, col=1)  # Methylation Scatter Y-axis

    fig.update_yaxes(title_text="D4Z4 Copies", row=4, col=2)  # Methylation Scatter X-axis
    fig.update_xaxes(title_text="Methylation (%)", row=4, col=2)  # Methylation Scatter Y-axis

    fig.update_xaxes(title_text="XapI Restriction Site", row=5, col=1)  # Restriction Sites X-axis
    fig.update_yaxes(title_text="BlnI Restriction Site", row=5, col=1)  # Restriction Sites Y-axis

    # Adjust layout
    fig.update_layout(
        title=dict(text=f"{sample} report", font=dict(size=40)),
        height=3500, 
        width=2200, 
        showlegend=True,
        plot_bgcolor="white",
        paper_bgcolor="white",
    )

    for row in range(1, 6):  
        for col in range(1, 6):  
            fig.update_xaxes(showline=True, linecolor="black", mirror=True, row=row, col=col)
            fig.update_yaxes(showline=True, linecolor="black", mirror=True, row=row, col=col)

    plotly.offline.plot(fig, filename=f"{dir_path}/{sample}_report.html")

# Run
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract alignment details and sequences from a BAM file.")
    parser.add_argument("--main_tsv", required=True, help="Main output TSV file to update with sensitive repeat counts.")
    args = parser.parse_args()

    file_path = args.main_tsv

    sample_name = file_path.split("/")[-1].split("_")[0]
    dir_path = "/".join(file_path.split("/")[:-1])

    # Load Data
    df = pd.read_csv(file_path, sep='\t')

    # Create a consistent color mapping
    unique_labels = df["ReadLabel"].dropna().unique()
    color_map = {label: px.colors.qualitative.Set1[i % len(px.colors.qualitative.Set1)] for i, label in enumerate(unique_labels)}

    combine_plots(df, dir_path, sample_name, unique_labels, color_map)

