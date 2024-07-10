import argparse
import pandas as pd
import numpy as np
import logging
import sys
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
import matplotlib
import seaborn as sb
import plotly.graph_objects as go
import plotly
import pysam

# Argument Parser
opt_parser = argparse.ArgumentParser()

opt_parser.add_argument(
    "-a",
    "--alignment_df",
    dest="alignment_df_name",
    help="Insert an alignment dataframe",
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
    "-c",
    "--color",
    dest="color_sample",
    help="Mention which color the output plot should have",
    metavar="FILE",
)

opt_parser.add_argument(
    "-o",
    "--output_path",
    dest="output",
    help="Insert an output directory to write to",
    metavar="FILE",
)

opt_parser.add_argument(
    "-r",
    "--reference",
    dest="reference",
    help="Insert a path to the reference file",
    metavar="FILE",
)

options = opt_parser.parse_args()

alignment_df_name = options.alignment_df_name
output = options.output
fragment_df_path = options.fragment_df_path
color_sample = options.color_sample
reference_path = options.reference

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

# Logger construction
logger = logging.getLogger(__name__)
stream = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(
    "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

stream.setFormatter(formatter)
file_stream = logging.FileHandler(f"{output}fragment_analysis.log")
file_stream.setFormatter(formatter)
logger.addHandler(stream)
logger.addHandler(file_stream)
logger.setLevel("INFO")


def plot_total_coverage_single_sample(
    a_df: pl.DataFrame,
    scale: str,
    length_of_reference: int,
    sample: str,
    fragment_df: pd.DataFrame,
):
    fragment_df = fragment_df[fragment_df["Fragment"].isin(["5-8S", "18S", "28S"])]
    coverage_counter_list = np.array([0 for nucleotide in range(length_of_reference)])
    coverage_index_list = [i + 1 for i in range(length_of_reference)]
    for start_position, end_position in zip(a_df["Refstart"], a_df["Refend"]):
        coverage_counter_list[start_position:end_position] += 1
    if scale == "absolute":
        mapper = map(lambda x: x, coverage_counter_list)
    elif scale == "relative":
        mapper = map(lambda x: x / int(a_df.shape[0]), coverage_counter_list)
    rel_coverage_counter_list = tuple(mapper)

    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 10))
    ax.plot(
        coverage_index_list,
        rel_coverage_counter_list,
        color=color_dict[sample],
        linewidth=2,
    )
    # ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), title="Sample")
    ax.set_xlabel("Position on reference")
    ax.set_ylabel(f"{scale.capitalize()} coverage")
    if scale == "absolute":
        ax.set_ylim(0, a_df.shape[0])
        ax.vlines(
            x=[
                [start, end]
                for start, end in zip(fragment_df["Start"], fragment_df["End"])
            ],
            ymin=0,
            colors=["red", "red"],
            ymax=a_df.shape[0],
            ls="--",
            lw=1,
            alpha=1,
        )
    elif scale == "relative":
        ax.set_ylim(0, 1)
        ax.vlines(
            x=[
                [start, end]
                for start, end in zip(fragment_df["Start"], fragment_df["End"])
            ],
            ymin=0,
            colors=["red", "red"],
            ymax=1,
            ls="--",
            lw=1,
            alpha=1,
        )
    plt.savefig(f"{output}/coverage_total_sample_{scale}.png", format="png")
    return fig


def plot_fragment_dynamics_single_sample(
    alignment_df: pl.DataFrame,
    length_of_reference: int,
    scale: str,
    fragment_df: pd.DataFrame,
):
    list_of_fragments = [
        "47S",
        "45S",
        "43S",
        "30S+1",
        "30S",
        "26S",
        "21S",
        "21S-C",
        "18S-E",
        "36S",
        "36S-C",
        "32S",
        "28.5S",
        "12S",
        "5-8S",
        "18S",
        "28S",
    ]
    half = len(list_of_fragments) // 2
    cmap = cm.get_cmap("tab20b").colors
    # list_of_fragments = np.unique(alignment_df["Associated_Fragments_Overlap"])
    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(len(list_of_fragments), hspace=1)
    axs = gs.subplots()
    for index, fragment in enumerate(list_of_fragments):
        temp_alignment_df = alignment_df.filter(
            pl.col("Associated_Fragments_Overlap") == fragment
        )

        if temp_alignment_df.shape[0] != 0:
            coverage_counter_list = np.array(
                [0 for nucleotide in range(length_of_reference)]
            )
            coverage_index_list = np.array([i + 1 for i in range(length_of_reference)])
            for start_position, end_position in zip(
                temp_alignment_df["Refstart"], temp_alignment_df["Refend"]
            ):
                coverage_counter_list[start_position:end_position] += 1
            if scale == "absolute":
                mapper = map(lambda x: x, coverage_counter_list)
            elif scale == "relative":
                mapper = map(
                    lambda x: x / int(temp_alignment_df.shape[0]), coverage_counter_list
                )
            elif scale == "absolute_all":
                mapper = map(
                    lambda x: x, coverage_counter_list
                )
            rel_coverage_counter_list = tuple(mapper)

            axs[index].plot(
                coverage_index_list,
                rel_coverage_counter_list,
                color=cmap[index],
                label=fragment,
            )
            axs[index].set_title(fragment)
            if scale == "absolute":
                axs[index].set_ylim(
                    0, temp_alignment_df.shape[0] + temp_alignment_df.shape[0] * 0.1
                )
                axs[index].vlines(
                    x=[
                        int(
                            fragment_df[fragment_df["Fragment"] == fragment][
                                "Start"
                            ].iloc[0]
                        ),
                        int(
                            fragment_df[fragment_df["Fragment"] == fragment][
                                "End"
                            ].iloc[0]
                        ),
                    ],
                    ymin=0,
                    colors=["red", "red"],
                    ymax=temp_alignment_df.shape[0] + temp_alignment_df.shape[0] * 0.1,
                    ls="--",
                    lw=1,
                    alpha=1,
                )
            elif scale == "relative":
                axs[index].set_ylim(0, 1.1)
                axs[index].vlines(
                    x=[
                        int(
                            fragment_df[fragment_df["Fragment"] == fragment][
                                "Start"
                            ].iloc[0]
                        ),
                        int(
                            fragment_df[fragment_df["Fragment"] == fragment][
                                "End"
                            ].iloc[0]
                        ),
                    ],
                    ymin=0,
                    colors=["red", "red"],
                    ymax=1.1,
                    ls="--",
                    lw=1,
                    alpha=1,
                )
            elif scale == "absolute_all":
                axs[index].set_ylim(0, alignment_df.shape[0])
                axs[index].vlines(
                    x=[
                        int(
                            fragment_df[fragment_df["Fragment"] == fragment][
                                "Start"
                            ].iloc[0]
                        ),
                        int(
                            fragment_df[fragment_df["Fragment"] == fragment][
                                "End"
                            ].iloc[0]
                        ),
                    ],
                    ymin=0,
                    colors=["red", "red"],
                    ymax=alignment_df.shape[0],
                    ls="--",
                    lw=1,
                    alpha=1,
                )

    for index, ax in enumerate(axs.flat):
        if index == half:
            mid_ax = ax

    for ax in fig.get_axes():
        ax.label_outer()
        last_ax = ax
    y_label = scale.capitalize().replace("_"," ")
    mid_ax.set(ylabel=f"{y_label} fragment coverage")
    last_ax.set(xlabel="Position on reference")
    plt.savefig(f"{output}/coverage_fragments_{scale}.png", format="png")


alignment_df = pl.read_csv(
    alignment_df_name,
    separator=";",
    columns=[
        "ID",
        "Length",
        "Refstart",
        "Refend",
        "n_Reads",
        "Overlap",
        "Matches",
        "Associated_Fragments_Overlap",
    ],
    has_header=True,
)
fasta_file = pysam.FastaFile(reference_path)
reference = fasta_file.references[0]
reference_sequence = str(fasta_file.fetch(reference))

fragment_df = pd.read_csv(fragment_df_path, sep=";", header=0, index_col=None)

plot_total_coverage_single_sample(
    alignment_df, "absolute", len(reference_sequence), color_sample, fragment_df
)

plot_total_coverage_single_sample(
    alignment_df, "relative", len(reference_sequence), color_sample, fragment_df
)

plot_fragment_dynamics_single_sample(
    alignment_df=alignment_df,
    length_of_reference=len(reference_sequence),
    scale="absolute",
    fragment_df=fragment_df,
)

plot_fragment_dynamics_single_sample(
    alignment_df=alignment_df,
    length_of_reference=len(reference_sequence),
    scale="absolute_all",
    fragment_df=fragment_df,
)

plot_fragment_dynamics_single_sample(
    alignment_df=alignment_df,
    length_of_reference=len(reference_sequence),
    scale="relative",
    fragment_df=fragment_df,
)
