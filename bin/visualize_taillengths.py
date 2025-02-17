from pathlib import Path
import numpy as np
import argparse
import polars as pl
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import plotly.graph_objects as go
import plotly
import pysam

# Takes output of taillength estimation
# and plots polyA taillength deviation on different 
# template-based algorithm associated rRNA intermediates.
# Horizontal barplots and violinplots are generated.

#####################################################################################################################################
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                    Argument parser                                                                #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#####################################################################################################################################

opt_parser = argparse.ArgumentParser()

opt_parser.add_argument(
    "-t",
    "--taillength_file",
    dest="taillength_name",
    help="Insert filepath with taillengths per read",
    metavar="FILE",
)
opt_parser.add_argument(
    "-a",
    "--template_file",
    dest="template_name",
    help="Insert filepath to filepath of template based analysis on per read basis",
    metavar="FILE",
)
opt_parser.add_argument(
    "-c",
    "--color",
    dest="color_sample",
    help="Mention which color the output plot should have",
    metavar="FILE",
)
opt_parser.add_argument(
    "-r",
    "--reference",
    dest="reference",
    help="Insert a path to the reference file",
    metavar="FILE",
)
opt_parser.add_argument(
    "-o",
    "--output_path",
    dest="output_path",
    help="Insert an output filepath to write to",
    metavar="FILE",
)

opt_parser.add_argument(
    "-f",
    "--fragment_df",
    dest="fragment_df_path",
    help="Insert a path to the fragment reference from literature",
    metavar="FILE",
)

opt_parser.add_argument(
    "-m",
    "--model_organism",
    dest="model_organism",
    help="Determine your model organism (Human,Yeast)",
    
)


options = opt_parser.parse_args()
taillength_file = options.taillength_name
template_file = options.template_name
reference_path = options.reference
output_path = options.output_path
color_sample = options.color_sample
fragment_df_path = options.fragment_df_path
model_organism = options.model_organism

#####################################################################################################################################
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                       Functions                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#####################################################################################################################################

def plot_reads_without_taillength(df:pd.DataFrame,fragment_df:pd.DataFrame, boundary_type: str):   
        if model_organism == "Human":
                fragment_df_plot = fragment_df[
                    fragment_df["Fragment"].isin(
                        ["18S", "28S", "5-8S", "5ETS", "ITS1", "ITS2", "ITS3", "3ETS"]
                    )
                ]
        elif model_organism == "Yeast":
                fragment_df_plot = fragment_df[
                    fragment_df["Fragment"].isin(
                        ["18S", "25S", "5-8S", "5ETS", "ITS1", "ITS2", "ITS3", "3ETS"]
                    )
                ]
        if boundary_type == "mean": 
            df_summary = df.groupby(
            by="Associated_Fragments_Overlap", as_index=False
            ).agg(
                {
                    "start": "mean",
                    "end": "mean",
                    "n_reads": "sum"
                }
            )
        elif boundary_type == "min_max":
            df_summary = df.groupby(
            by="Associated_Fragments_Overlap", as_index=False
            ).agg(
                {
                    "start": "min",
                    "end": "max",
                    "n_reads": "sum"
                }
            )
        elif boundary_type == "template": 
            df_summary = df.groupby(
            by="Associated_Fragments_Overlap", as_index=False
            ).agg(
                {
                    "start": "min",
                    "end": "max",
                    "n_reads": "sum"
                }
            )
            template_starts = []
            template_ends = []
            for assoc_fragment in df_summary["Associated_Fragments_Overlap"]:
                start = int(fragment_df.loc[fragment_df["Fragment"] == assoc_fragment,"Start"].values[0])
                end = int(fragment_df.loc[fragment_df["Fragment"] == assoc_fragment,"End"].values[0])
                template_starts.append(start)
                template_ends.append(end)
            df_summary["start"] = template_starts
            df_summary["end"] = template_ends
            
        df_summary = df_summary.sort_values(
            by=["end", "start"], ascending=[False, True]
        )
        df_summary = pd.DataFrame(df_summary).fillna(value=0)
        max_n_reads = df_summary["n_reads"].max()
        if max_n_reads == 0:
            max_n_reads = 1
        
        layout = go.Layout(height = 800)
        fig = go.Figure(layout=layout)

        fig.add_trace(
            go.Scatter(
                x=[i for i in range(0, len(reference_sequence))],
                y=[0 for i in range(0, len(reference_sequence))],
                name="O-line",
                line_color="white",
                line_width=0.1,
            )
        )

        for (
            associated_template,
            index,
            start,
            end,
            n_reads,
        ) in zip(
            df_summary["Associated_Fragments_Overlap"],
            [i for i in range(len(df_summary["Associated_Fragments_Overlap"]))],
            df_summary["start"],
            df_summary["end"],
            df_summary["n_reads"],
        ):
            fig.add_shape(
                x0=start,
                x1=end,
                y0=index + 1,
                y1=index + 1,
                opacity=min((n_reads / max_n_reads) + 0.1, 1),
                line=dict(color=color_dict[color_sample], width=5),
            )
            fig.update_layout(
                xaxis=dict(title="Position on reference", gridcolor="white",tickformat="d"),
                yaxis=dict(
                    title="Cluster",
                    gridcolor="white",
                    zeroline=True,
                    zerolinecolor="white",
                    zerolinewidth=0.1,
                    tickformat="d"
                ),
                plot_bgcolor="rgba(0,0,0,0)",
            )
        
        for index, row in fragment_df_plot.iterrows():
            fig.add_shape(
                x0=row["Start"],
                x1=row["Start"],
                y0=0,
                y1=df_summary.shape[0] + 1,
                opacity=min((n_reads / max_n_reads) + 0.9, 1),
                line=dict(color="red", width=0.1, dash="dash"),
            )
            fig.add_shape(
                x0=row["End"],
                x1=row["End"],
                y0=0,
                y1=df_summary.shape[0] + 1,
                opacity=min((n_reads / max_n_reads) + 0.9, 1),
                line=dict(color="red", width=0.1, dash="dash"),
            )
            fig.add_annotation(
                x=(row["Start"] + row["End"]) / 2,
                y=df_summary.shape[0] + 3,
                text=row["Fragment"],
                arrowhead=1,
                showarrow=False,
                ax=0,
                ay=-40,
                font=dict(size=12),
            )
            fig.add_annotation(
                x=(row["Start"] + row["End"]) / 2,
                y=-3,
                text=row["Fragment"],
                arrowhead=1,
                showarrow=False,
                ax=0,
                ay=-40,
                font=dict(size=12),
            )
        
        if boundary_type == "mean":
            plotly.io.write_html(fig, f"{output_path}/polyA_tails_intermediates_mean.html")
        elif boundary_type == "min_max":
            plotly.io.write_html(fig, f"{output_path}/polyA_tails_intermediates_min_max.html")
        elif boundary_type == "template":
            plotly.io.write_html(fig, f"{output_path}/polyA_tails_intermediates_template.html")
            
def plot_reads_with_taillength(df:pd.DataFrame,fragment_df:pd.DataFrame, boundary_type: str):   
        if model_organism == "Human":
                fragment_df_plot = fragment_df[
                    fragment_df["Fragment"].isin(
                        ["18S", "28S", "5-8S", "5ETS", "ITS1", "ITS2", "ITS3", "3ETS"]
                    )
                ]
        elif model_organism == "Yeast":
                fragment_df_plot = fragment_df[
                    fragment_df["Fragment"].isin(
                        ["18S", "25S", "5-8S", "5ETS", "ITS1", "ITS2", "ITS3", "3ETS"]
                    )
                ]
        df_short = df.loc[
            :,
            [
                "Associated_Fragments_Overlap",
                "taillength",
                "taillength",
                "taillength",
                "Refstart",
                "Refend",
                "n_Reads",
                ],
            ]
        df_short.columns = [
            "Associated_Fragments_Overlap",
            "mean_taillength",
            "median_taillength",
            "std_taillength",
            "start",
            "end",
            "n_reads",
        ]
        df_summary_test = df_short.groupby(
            by="Associated_Fragments_Overlap", 
            as_index=False
            ).agg(
                {
                    "mean_taillength": "mean",
                    "median_taillength": "median",
                    "std_taillength": "std",
                    "start": "mean",
                    "end": "mean",
                    "n_reads": "sum",
                }
            )
        fragments_to_clean = df_summary_test["Associated_Fragments_Overlap"][df_summary_test["mean_taillength"] != 0]
        for fragment_to_clean in np.unique(fragments_to_clean):
            df_short = df_short[~((df_short["Associated_Fragments_Overlap"] == fragment_to_clean) & (df_short["mean_taillength"] == 0))]
        if boundary_type == "mean": 
            df_summary = df_short.groupby(
                by="Associated_Fragments_Overlap", as_index=False
            ).agg(
                {
                    "mean_taillength": "mean",
                    "median_taillength": "median",
                    "std_taillength": "std",
                    "start": "mean",
                    "end": "mean",
                    "n_reads": "sum",
                }
            )
            
        elif boundary_type == "min_max":
            df_summary = df_short.groupby(
                by="Associated_Fragments_Overlap", as_index=False
            ).agg(
                {
                    "mean_taillength": "mean",
                    "median_taillength": "median",
                    "std_taillength": "std",
                    "start": "min",
                    "end": "max",
                    "n_reads": "sum",
                }
            )

            df_summary.to_csv(
                f"{output_path}/taillength_per_intermediate_min_max.csv",
                sep="\t",
                header=True,
                index=None,
            )
        
        elif boundary_type == "template": 
            df_summary = df_short.groupby(
                by="Associated_Fragments_Overlap", as_index=False
            ).agg(
                {
                    "mean_taillength": "mean",
                    "median_taillength": "median",
                    "std_taillength": "std",
                    "start": "mean",
                    "end": "mean",
                    "n_reads": "sum",
                }
            )
            
            template_starts = []
            template_ends = []
            for assoc_fragment in df_summary["Associated_Fragments_Overlap"]:
                start = int(fragment_df.loc[fragment_df["Fragment"] == assoc_fragment,"Start"].values[0])
                end = int(fragment_df.loc[fragment_df["Fragment"] == assoc_fragment,"End"].values[0])
                template_starts.append(start)
                template_ends.append(end)

            df_summary["start"] = template_starts
            df_summary["end"] = template_ends
        
        df_summary = df_summary.sort_values(
        by=["end", "start"], ascending=[False, True]
        )
        df_summary = pd.DataFrame(df_summary).fillna(value=0)
            
        max_n_reads = df_summary["n_reads"].max()
        if max_n_reads == 0:
            max_n_reads = 1
        
        
        layout = go.Layout(height = 800)
        fig = go.Figure(layout=layout)
        
        for (
        associated_template,
        index,
        mean_taillength,
        median_taillength,
        std_taillength,
        start,
        end,
        n_reads,
        ) in zip(
            df_summary["Associated_Fragments_Overlap"],
            [i for i in range(len(df_summary["mean_taillength"]))],
            df_summary["mean_taillength"],
            df_summary["median_taillength"],
            df_summary["std_taillength"],
            df_summary["start"],
            df_summary["end"],
            df_summary["n_reads"],
        ):
            fig.add_shape(
                x0=start,
                x1=end,
                y0=index + 1,
                y1=index + 1,
                opacity=min((n_reads / max_n_reads) + 0.1, 1),
                line=dict(color=color_dict[color_sample], width=5),
            )
            fig.add_shape(
                x0=end,
                x1=end + mean_taillength,
                y0=index + 1,
                y1=index + 1,
                opacity=min((n_reads / max_n_reads) + 0.9, 1),
                line=dict(color="burlywood", width=4),
            )
            fig.add_shape(
                x0=end + mean_taillength - (std_taillength),
                x1=end + mean_taillength + (std_taillength),
                y0=index + 1,
                y1=index + 1,
                opacity=min((n_reads / max_n_reads) + 0.9, 1),
                line=dict(color="brown", width=0.7),
            )
            fig.add_trace(
                go.Scatter(
                    x=[end + mean_taillength + (std_taillength)],
                    y=[index + 1],  # Use the same x position and y at the top of the line
                    mode="markers",
                    marker=dict(size=10, color="rgba(0,0,0,0)"),  # Invisible marker
                    showlegend=False,
                    hovertext=f"Fragment:{associated_template},\nStart:{int(start)}\nEnd:{int(end)}\n,Taillength:{int(mean_taillength)}\n,Stdd taillength:{int(std_taillength)}\n#Reads:{n_reads}\n",  # Hover text
                )
            )
            fig.update_layout(
                xaxis=dict(title="Position on reference", gridcolor="white",tickformat="d"),
                yaxis=dict(
                    title="Cluster",
                    gridcolor="white",
                    zeroline=True,
                    zerolinecolor="white",
                    zerolinewidth=0.1,
                    tickformat="d"
                ),
                plot_bgcolor="rgba(0,0,0,0)",
            )

        for index, row in fragment_df_plot.iterrows():
            fig.add_shape(
                x0=row["Start"],
                x1=row["Start"],
                y0=0,
                y1=df_summary.shape[0] + 1,
                opacity=min((n_reads / max_n_reads) + 0.9, 1),
                line=dict(color="red", width=0.1, dash="dash"),
            )
            fig.add_shape(
                x0=row["End"],
                x1=row["End"],
                y0=0,
                y1=df_summary.shape[0] + 1,
                opacity=min((n_reads / max_n_reads) + 0.9, 1),
                line=dict(color="red", width=0.1, dash="dash"),
            )
            fig.add_annotation(
                x=(row["Start"] + row["End"]) / 2,
                y=df_summary.shape[0] + 3,
                text=row["Fragment"],
                arrowhead=1,
                showarrow=False,
                ax=0,
                ay=-40,
                font=dict(size=12),
            )
            fig.add_annotation(
                x=(row["Start"] + row["End"]) / 2,
                y=-3,
                text=row["Fragment"],
                arrowhead=1,
                showarrow=False,
                ax=0,
                ay=-40,
                font=dict(size=12),
            )
        
        if boundary_type == "mean":
            plotly.io.write_html(fig, f"{output_path}/polyA_tails_intermediates_mean.html")
        elif boundary_type == "min_max":
            plotly.io.write_html(fig, f"{output_path}/polyA_tails_intermediates_min_max.html")
        elif boundary_type == "template":
            plotly.io.write_html(fig, f"{output_path}/polyA_tails_intermediates_template.html")     
            
            
def filter_outliers(
    df, column, lower_percentile=0.1, upper_percentile=0.9, multiplier=2
    ):
        """
        Filters out outliers from a DataFrame based on the Interquartile Range (IQR) method.
        
        Parameters:
        df (pd.DataFrame): The input DataFrame.
        column (str): The name of the column for which outliers are to be filtered.
        lower_percentile (float): The lower percentile to calculate Q1 (default is 0.25).
        upper_percentile (float): The upper percentile to calculate Q3 (default is 0.75).
        multiplier (float): The multiplier for the IQR to define the cutoff boundaries (default is 2).
        
        Returns:
        pd.DataFrame: A DataFrame with outliers filtered out based on the specified criteria.
        
        Example:
        >>> data = {'Category': ['A', 'B', 'C', 'A', 'B', 'C'], 'Value': [1, 5, 3, 2, 6, 4]}
        >>> df = pd.DataFrame(data)
        >>> filtered_df = filter_outliers(df, 'Value')
        >>> print(filtered_df)
        """
        # Calculate the IQR
        Q1 = df[column].quantile(lower_percentile)
        Q3 = df[column].quantile(upper_percentile)
        IQR = Q3 - Q1

        # Define the cutoff boundaries
        lower_bound = Q1 - (IQR * multiplier)
        upper_bound = Q3 + (IQR * multiplier)

        # Filter the data
        filtered_df = df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]
        return filtered_df, Q1, Q3   

#####################################################################################################################################
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                           Script                                                                  #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#####################################################################################################################################

fasta_file = pysam.FastaFile(reference_path)
reference = fasta_file.references[0]
reference_sequence = str(fasta_file.fetch(reference))

taillength_df = pl.read_csv(
    taillength_file, separator=";", columns=["read_id", "taillength"], has_header=True
)

template_df = pl.read_csv(
    template_file,
    separator=";",
    columns=[
        "ID",
        "Overlap",
        "Associated_Fragments_Overlap",
        "Refstart",
        "Refend",
        "n_Reads",
    ],
    has_header=True,
)

joined_df = template_df.join(
    taillength_df, left_on="ID", right_on="read_id", how="inner"
).to_pandas()
joined_df = joined_df.dropna()

fragment_df = pd.read_csv(fragment_df_path, sep=";", header=0, index_col=None)

if model_organism == "Human":
    fragment_df_plot = fragment_df[
        fragment_df["Fragment"].isin(
            ["18S", "28S", "5-8S", "5ETS", "ITS1", "ITS2", "ITS3", "3ETS"]
        )
    ]
elif model_organism == "Yeast":
        fragment_df_plot = fragment_df[
            fragment_df["Fragment"].isin(
                ["18S", "25S", "5-8S", "5ETS", "ITS1", "ITS2", "ITS3", "3ETS"]
            )
        ]




color_dict = {}
for sample, color in zip(
    ["Nucleus", "Cytoplasm", "SN1", "SN2", "SN3"],
    ["navy", "darkred", "seagreen", "darkorange", "darkviolet"],
):
    color_dict[sample] = color
for sample, color in zip(
    ["blue", "red", "green", "orange", "purple"],
    ["navy", "darkred", "seagreen", "darkorange", "darkviolet"],
):
    color_dict[sample] = color

if model_organism == "Human":
    formatter_list = [
        "18S-E",
        "21S-C",
        "21S",
        "26S",
        "30S",
        "36S",
        "32S",
        "28.5S",
        "12S",
        "7S",
        "5-8S",
        "18S",
        "28S",
    ]  # ,"5ETS","5ETS-01","ITS1","ITS2","3ETS"]
elif model_organism == "Yeast":
    formatter_list = [
        "25S",
        "23S",
        "22S",
        "21S",
        "20S",
        "32S",
        "27SA",
        "27SB",
        "25.5S",
        "7S",
        "6S",
        "5-8S",
        "18S",
        "25S"
    ]  # ,"5ETS","5ETS-01","ITS1","ITS2","3ETS"]

if joined_df.shape[0] == 0:
    joined_df = pd.DataFrame({"Associated_Fragments_Overlap":[""],"ID":["0:0:0"], "position":[0], "start":[0], "end":[0], "n_reads":[0], "rel_n_reads":[0],"taillength":[0]})    
    template_df = template_df.to_pandas()
    template_df.columns = [
        "ID",
        "start",
        "end",
        "n_reads",
        "Overlap",
        "Associated_Fragments_Overlap",
    ]
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=500)
    sns.violinplot(
        data=template_df,
        x="Associated_Fragments_Overlap",
        y=[0 for i in template_df["Associated_Fragments_Overlap"]],
        order=formatter_list,
        ax=ax,
        inner_kws=dict(box_width=2, whis_width=1, color="0.5"),
        color=color_dict[color_sample],
    )
    
    group_sizes = (
        template_df["Associated_Fragments_Overlap"].value_counts().sort_index().reset_index()
    )
    n_reads_per_group = []
    for element in formatter_list:
        n_reads_for_group = list(
            group_sizes.loc[group_sizes["Associated_Fragments_Overlap"] == element, "count"]
        )
        if n_reads_for_group == []:
            n_reads_for_group = 0
        else:
            n_reads_for_group = int(n_reads_for_group[0])
        n_reads_per_group.append(n_reads_for_group)

    for i, v in enumerate(n_reads_per_group):
        ax.text(i, int(max(joined_df["taillength"]) * 1.02), f"n={v}", ha="center", va="bottom", c="black", fontdict=dict(size=5, rotation=45))
    ax.set_xlabel("Ribosomal intermediate")
    ax.set_ylabel("PolyA-tail length")
    fig.savefig(f"{output_path}/violinplot_taillength_per_intermediate.png", format="png")
    
    for boundary_type in ["mean","min_max","template"]:
        plot_reads_without_taillength(template_df,fragment_df,boundary_type)
    
else:
    joined_df = template_df.join(
    taillength_df, left_on="ID", right_on="read_id", how="full"
    ).to_pandas()
    joined_df = joined_df.fillna(0)
    joined_df = joined_df[joined_df["Associated_Fragments_Overlap"] != 0]
    joined_df, lower_quantile, upper_quantile = filter_outliers(joined_df, "taillength")
    joined_df_summary_test = joined_df.groupby(
        by="Associated_Fragments_Overlap", 
        as_index=False
        ).agg(
            {
                "taillength": "mean"
            }
        )
    fragments_to_clean = joined_df_summary_test["Associated_Fragments_Overlap"][joined_df_summary_test["taillength"] != 0]
    for fragment_to_clean in np.unique(fragments_to_clean):
        joined_df = joined_df[~((joined_df["Associated_Fragments_Overlap"] == fragment_to_clean) & (joined_df["taillength"] == 0))]
    fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=500)
    sns.violinplot(
        data=joined_df,
        x="Associated_Fragments_Overlap",
        y="taillength",
        order=formatter_list,
        ax=ax,
        inner_kws=dict(box_width=2, whis_width=1, color="0.5"),
        color=color_dict[color_sample],
    )
    group_sizes = (
        joined_df["Associated_Fragments_Overlap"].value_counts().sort_index().reset_index()
    )
    n_reads_per_group = []
    for element in formatter_list:
        n_reads_for_group = list(
            group_sizes.loc[group_sizes["Associated_Fragments_Overlap"] == element, "count"]
        )
        if n_reads_for_group == []:
            n_reads_for_group = 0
        else:
            n_reads_for_group = int(n_reads_for_group[0])
        n_reads_per_group.append(n_reads_for_group)

    for i, v in enumerate(n_reads_per_group):
        ax.text(i, int(max(joined_df["taillength"]) * 1.02), f"n={v}", ha="center", va="bottom", c="black", fontdict=dict(size=5, rotation=45))
    ax.set_ylim(0, int(max(joined_df["taillength"]) * 1.1))
    ax.set_yticks([i for i in range(0, int(max(joined_df["taillength"]) * 1.1) , 5)])
    ax.set_yticklabels([str(i) if (i % 10) == 0 else "" for i in range(0, int(max(joined_df["taillength"]) * 1.1), 5)], fontsize=5)
    ax.set_xlabel("Ribosomal intermediate")
    ax.set_ylabel("PolyA-tail length")
    fig.savefig(f"{output_path}/violinplot_taillength_per_intermediate.png", format="png")
    
    joined_df = template_df.join(
    taillength_df, left_on="ID", right_on="read_id", how="full"
    ).to_pandas()
    joined_df = joined_df.fillna(0)
    joined_df = joined_df[joined_df["Associated_Fragments_Overlap"] != 0]
    joined_df, lower_quantile, upper_quantile = filter_outliers(joined_df, "taillength")
    
    for boundary_type in ["mean","min_max","template"]:
        plot_reads_with_taillength(joined_df,fragment_df,boundary_type)