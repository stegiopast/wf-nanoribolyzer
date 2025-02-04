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
# template-free clustered read groups.
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
    "-o",
    "--output_path",
    dest="output_path",
    help="Insert an output filepath to write to",
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
output_path = options.output_path
color_sample = options.color_sample
reference_path = options.reference
fragment_df_path = options.fragment_df_path
model_organism = options.model_organism

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

fragment_df = pd.read_csv(fragment_df_path, sep=";", header=0, index_col=None)

if model_organism == "Human":
    fragment_df = fragment_df[
        fragment_df["Fragment"].isin(
            ["18S", "28S", "5-8S", "5ETS", "ITS1", "ITS2", "ITS3", "3ETS"]
        )
    ]
elif model_organism == "Yeast":
        fragment_df = fragment_df[
            fragment_df["Fragment"].isin(
                ["18S", "25S", "5-8S", "5ETS", "ITS1", "ITS2", "ITS3", "3ETS"]
            )
        ]

pre_template_df = pl.read_csv(
    template_file,
    separator=";",
    columns=["ID", "IDS", "Refstart", "Refend", "n_Reads", "rel_n_Reads"],
    schema_overrides={
        "ID": pl.Utf8,
        "IDS": pl.Utf8,
        "Refend": pl.Int32,
        "n_Reads": pl.Int32,
        "rel_n_Reads": pl.Float32,
    },
    has_header=True,
)


pre_template_df = pre_template_df.sort(by=["rel_n_Reads"], descending=[True])[
    0 : min(300, pre_template_df.shape[0])
]
pre_template_df = pre_template_df.sort(by=["Refend"], descending=[False])

template_dict = []
for position, ids, refstart, refend, n_reads, rel_n_reads in zip(
    pre_template_df["ID"],
    pre_template_df["IDS"],
    pre_template_df["Refstart"],
    pre_template_df["Refend"],
    pre_template_df["n_Reads"],
    pre_template_df["rel_n_Reads"],
):
    id_list = eval(ids)
    for id in id_list:
        temp_dict = {
            "ID": id,
            "position": position,
            "start": refstart,
            "end": refend,
            "n_reads": n_reads,
            "rel_n_reads": rel_n_reads,
        }
        template_dict.append(temp_dict)

template_df = pl.DataFrame(template_dict)


joined_df = template_df.join(
    taillength_df, left_on="ID", right_on="read_id", how="inner"
).to_pandas()
joined_df = joined_df.dropna()


# def filter_outliers(
#     df, column, lower_percentile=0.25, upper_percentile=0.75, multiplier=2
# ):
#     """
#     Filters out outliers from a DataFrame based on the Interquartile Range (IQR) method.
    
#     Parameters:
#     df (pd.DataFrame): The input DataFrame.
#     column (str): The name of the column for which outliers are to be filtered.
#     lower_percentile (float): The lower percentile to calculate Q1 (default is 0.25).
#     upper_percentile (float): The upper percentile to calculate Q3 (default is 0.75).
#     multiplier (float): The multiplier for the IQR to define the cutoff boundaries (default is 2).
    
#     Returns:
#     pd.DataFrame: A DataFrame with outliers filtered out based on the specified criteria.
    
#     Example:
#     >>> data = {'Category': ['A', 'B', 'C', 'A', 'B', 'C'], 'Value': [1, 5, 3, 2, 6, 4]}
#     >>> df = pd.DataFrame(data)
#     >>> filtered_df = filter_outliers(df, 'Value')
#     >>> print(filtered_df)
#     """
#     # Calculate the IQR
#     Q1 = df[column].quantile(lower_percentile)
#     Q3 = df[column].quantile(upper_percentile)
#     IQR = Q3 - Q1

#     # Define the cutoff boundaries
#     lower_bound = Q1 - (IQR * multiplier)
#     upper_bound = Q3 + (IQR * multiplier)

#     # Filter the data
#     filtered_df = df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]
#     return filtered_df, Q1, Q3


#joined_df, lower_quantile, upper_quantile = filter_outliers(joined_df, "taillength")
# print(joined_df)

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


ids_for_violin_plot = [id for id in pre_template_df["ID"][0:20]]
violin_joined_df = joined_df[joined_df["position"].isin(ids_for_violin_plot)]
violin_joined_df = violin_joined_df.sort_values(by=["end","start"],ascending=[True,True])

if violin_joined_df.shape[0] == 0:
    violin_joined_df = pd.DataFrame({"ID":["0:0:0"], "position":[0], "start":[0], "end":[0], "n_reads":[0], "rel_n_reads":[0],"taillength":[0]})


fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=500)
sns.violinplot(
    data=violin_joined_df,
    x="position",
    y="taillength",
    ax=ax,
    inner_kws=dict(box_width=2, whis_width=1, color="0.5"),
    color=color_dict[color_sample],
)
group_sizes = violin_joined_df["position"].value_counts().sort_index().reset_index()


n_reads_per_group = []
for element in violin_joined_df["position"].unique():
    n_reads_for_group = list(
        group_sizes.loc[group_sizes["position"] == element, "count"]
    )
    if n_reads_for_group == []:
        n_reads_for_group = 0
    else:
        n_reads_for_group = int(n_reads_for_group[0])
    n_reads_per_group.append(n_reads_for_group)

for i, v in enumerate(n_reads_per_group):
    ax.text(
        i,
        int(max(violin_joined_df["taillength"]) * 1.02),
        f"n={v}",
        ha="center",
        va="bottom",
        c="black",
        fontdict=dict(size=4, rotation=45),
    )


ax.set_ylim(0,  int(max(violin_joined_df["taillength"]) * 1.1))
ax.set_yticks([i for i in range(0,  int(max(violin_joined_df["taillength"]) * 1.1), 5)])
ax.set_yticklabels([str(i) if (i % 10) == 0 else "" for i in range(0,  int(max(violin_joined_df["taillength"]) * 1.1), 5)], fontsize=5)
ax.set_xticks(ax.get_xticks())
ax.set_xticklabels(ax.get_xticklabels(), fontsize=5, rotation=45)
ax.set_xlabel("Ribosomal intermediate")
ax.set_ylabel("PolyA-tail length")
fig.savefig(f"{output_path}/violinplot_taillength_per_cluster.png", format="png")
plt.close()


joined_df_short = joined_df.loc[
    :,
    [
        "position",
        "taillength",
        "taillength",
        "taillength",
        "start",
        "end",
        "rel_n_reads",
        "n_reads",
    ],
]
joined_df_short.columns = [
    "Associated_Fragments_Overlap",
    "mean_taillength",
    "median_taillength",
    "std_taillength",
    "start",
    "end",
    "rel_n_reads",
    "n_reads",
]

joined_df_summary = joined_df_short.groupby(
    by="Associated_Fragments_Overlap", as_index=False
).agg(
    {
        "mean_taillength": "mean",
        "median_taillength": "median",
        "std_taillength": "std",
        "start": "first",
        "end": "first",
        "rel_n_reads": "first",
        "n_reads": "first",
    }
)


joined_df_summary.to_csv(
    f"{output_path}/taillength_per_cluster.csv", sep="\t", header=True, index=None
)

layout = go.Layout(height = 800)
fig = go.Figure(layout=layout)

fig.add_trace(
    go.Scatter(
        x=[i for i in range(0, len(reference_sequence))],
        y=[0 for i in range(0, len(reference_sequence))],
        name="",
        line_color="white",
        line_width=0.1,
    )
)

max_n_reads = joined_df_summary["n_reads"].max()

joined_df_summary = joined_df_summary.sort_values(
    by=["end", "start"], ascending=[False, False]
)

joined_df_summary = pd.DataFrame(joined_df_summary).fillna(value=0)

for (
    index,
    mean_taillength,
    median_taillength,
    std_taillength,
    start,
    end,
    rel_n_reads,
    n_reads,
) in zip(
    [i for i in range(len(joined_df_summary["mean_taillength"]))],
    joined_df_summary["mean_taillength"],
    joined_df_summary["median_taillength"],
    joined_df_summary["std_taillength"],
    joined_df_summary["start"],
    joined_df_summary["end"],
    joined_df_summary["rel_n_reads"],
    joined_df_summary["n_reads"],
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
            hovertext=f"Start:{start}\nEnd:{end}\n,Taillength:{int(mean_taillength)}\n,Stdd taillength:{int(std_taillength)}\n#Reads:{n_reads}\n",  # Hover text
        )
    )
    fig.update_layout(
        xaxis=dict(title="Position on reference", gridcolor="white", tickformat="d"),
        yaxis=dict(
            title="Cluster",
            gridcolor="white",
            zeroline=True,
            zerolinecolor="white",
            zerolinewidth=0.1,
            tickformat="d",
        ),
        plot_bgcolor="rgba(0,0,0,0)",
    )

for index, row in fragment_df.iterrows():
    fig.add_shape(
        x0=row["Start"],
        x1=row["Start"],
        y0=0,
        y1=joined_df_summary.shape[0] + 1,
        opacity=min((n_reads / max_n_reads) + 0.9, 1),
        line=dict(color="red", width=0.1, dash="dash"),
    )
    fig.add_shape(
        x0=row["End"],
        x1=row["End"],
        y0=0,
        y1=joined_df_summary.shape[0] + 1,
        opacity=min((n_reads / max_n_reads) + 0.9, 1),
        line=dict(color="red", width=0.1, dash="dash"),
    )
    fig.add_annotation(
        x=(row["Start"] + row["End"]) / 2,
        y=joined_df_summary.shape[0] + 3,
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

plotly.io.write_html(fig, f"{output_path}/polyA_tails_clustering.html")
