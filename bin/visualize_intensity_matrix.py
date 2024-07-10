import argparse
import pandas as pd
import numpy as np
import logging
import sys
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sb
import plotly.graph_objects as go
import plotly


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
    "-t",
    "--template_df_name",
    dest="template_df",
    help="Insert dataframe with ribosomal RNA templates",
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


options = opt_parser.parse_args()

alignment_df_name = options.alignment_df_name
output = options.output
template_df_name = options.template_df
color_sample = options.color_sample


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


#####################################################################################################################################
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                          Functions                                                                #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#####################################################################################################################################


def create_intensity_matrix(fusion_alignment_df: pl.DataFrame):
    """
    Create an intensity matrix and an ID dictionary from a given fusion alignment DataFrame.

    Parameters:
    -----------
    fusion_alignment_df : pl.DataFrame
        A Polars DataFrame containing fusion alignment data with columns "Refstart", "Refend", and "ID".

    Returns:
    --------
    dbscan_matrix : np.array
        A NumPy array representing the DBSCAN matrix with columns: start points, end points, and intensity values.
    id_dict : dict
        A dictionary where keys are string representations of start and end points in the format "start:end", 
        and values are lists of IDs corresponding to these start and end points.

    This function performs the following steps:
    1. Initializes an intensity matrix with dimensions based on the maximum value of "Refend" from the DataFrame.
    2. Creates an empty dictionary to store IDs corresponding to each unique start and end point pair.
    3. Iterates over the rows of the DataFrame to update the intensity matrix and the ID dictionary.
    4. Constructs a DBSCAN matrix from the non-zero values of the intensity matrix.
    5. Returns the DBSCAN matrix and the ID dictionary.

    Example:
    --------
    >>> fusion_alignment_df = pl.DataFrame({
            "Refstart": [1, 2, 3],
            "Refend": [4, 5, 6],
            "ID": ["a", "b", "c"]
        })
    >>> dbscan_matrix, id_dict = create_intensity_matrix(fusion_alignment_df)
    """
    intensity_matrix = np.zeros(
        (
            max(fusion_alignment_df["Refend"]) + 1,
            max(fusion_alignment_df["Refend"]) + 1,
        ),
        dtype=int,
    )
    id_dict = {}
    for x, y, id in zip(
        fusion_alignment_df["Refstart"],
        fusion_alignment_df["Refend"],
        fusion_alignment_df["ID"],
    ):
        if f"{x}:{y}" not in id_dict.keys():
            id_dict[f"{x}:{y}"] = []
        intensity_matrix[x][y] += 1
        id_dict[f"{x}:{y}"].append(id)
    intensity_matrix.sum()
    id_dict.keys()
    dbscan_matrix = np.array(
        [[x, y, intensity_matrix[x][y]] for x, y in zip(*np.nonzero(intensity_matrix))]
    )
    return dbscan_matrix, id_dict


def plot_matrix(
    dbscan_matrix: np.array, color_sample: str, template_df_name: str, id_dict: dict
):
    """
    Plot the DBSCAN matrix with overlayed template fragments and save the plot and data.

    Parameters:
    -----------
    dbscan_matrix : np.array
        A NumPy array representing the DBSCAN matrix with columns: start points, end points, and intensity values.
    color_sample : str
        The sample type to determine the color of the scatter plot points.
    template_df_name : str
        The file name of the template DataFrame in CSV format.
    id_dict : dict
        A dictionary where keys are string representations of start and end points in the format "start:end", 
        and values are lists of IDs corresponding to these start and end points.

    This function performs the following steps:
    1. Reads the template DataFrame from the provided CSV file.
    2. Filters the template DataFrame to include specific fragments.
    3. Constructs a DataFrame from the DBSCAN matrix and calculates alpha values for scatter plot transparency.
    4. Sets up color mapping for the scatter plot based on the sample type.
    5. Plots the start and end points with varying transparency based on intensity values.
    6. Overlays horizontal and vertical dashed lines to indicate fragment boundaries from the template DataFrame.
    7. Annotates the plot with fragment names.
    8. Saves the plot as a PNG file and the DBSCAN data with IDs as a CSV file.

    Example:
    --------
    >>> dbscan_matrix = np.array([[1, 2, 10], [3, 4, 20], [5, 6, 15]])
    >>> color_sample = "Nucleus"
    >>> template_df_name = "template.csv"
    >>> id_dict = {"1:2": ["id1", "id2"], "3:4": ["id3"], "5:6": ["id4", "id5"]}
    >>> plot_matrix(dbscan_matrix, color_sample, template_df_name, id_dict)
    """
    template_df = pd.read_csv(template_df_name, sep=";", header=0, index_col=None)
    template_df = template_df[
        template_df["Fragment"].isin(
            ["18S", "28S", "5-8S", "5ETS", "ITS1", "ITS2", "ITS3", "3ETS"]
        )
    ]
    dbscan_df = pd.DataFrame(dbscan_matrix)
    start_points = dbscan_df.iloc[:, 0]
    end_points = dbscan_df.iloc[:, 1]
    maximum = dbscan_df.iloc[:, 2].max()
    alphas = dbscan_df.iloc[:, 2] / maximum
    alphas = alphas + 0.01
    alphas[alphas > 1] = 1

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
    fig, ax = plt.subplots(1, 1, figsize=(10, 10), dpi=500)
    ax.scatter(
        x=start_points, y=end_points, alpha=alphas, s=1, c=color_dict[color_sample]
    )
    for template, start, end in zip(
        template_df["Fragment"], template_df["Start"], template_df["End"]
    ):
        ax.hlines(
            y=end, xmin=start, xmax=end, color="red", linewidth=1, linestyles="dashed"
        )
        ax.vlines(
            x=start, ymin=start, ymax=end, color="red", linewidth=1, linestyles="dashed"
        )
        ax.hlines(
            y=start, xmin=start, xmax=end, color="red", linewidth=1, linestyles="dashed"
        )
        ax.vlines(
            x=end, ymin=start, ymax=end, color="red", linewidth=1, linestyles="dashed"
        )
        text_x_position = int(start + ((end - start) * 0.5)) - (len(template) * 50)
        text_y_position = start - 250
        ax.text(x=text_x_position, y=text_y_position, s=template, fontsize="x-small")
    ax.set_xlabel("Start sites")
    ax.set_ylabel("End sites")
    plt.savefig(f"{output}/intensity_matrix.png", format="png")
    dbscan_df.columns = ["start", "end", "n_reads"]
    ids = []
    for x, y in zip(dbscan_df["start"], dbscan_df["end"]):
        ids.append(list(id_dict[f"{x}:{y}"]))
    dbscan_df["ids"] = ids
    dbscan_df.to_csv(f"{output}/intensity_matrix.csv", sep=";", index=None)

    filtered_dbscan_df = dbscan_df.sort_values(by="n_reads", ascending=False)[
        0 : min(25000, dbscan_df.shape[0])
    ]
    start_points = filtered_dbscan_df.iloc[:, 0]
    end_points = filtered_dbscan_df.iloc[:, 1]
    maximum = filtered_dbscan_df.iloc[:, 2].max()
    alphas = filtered_dbscan_df.iloc[:, 2] / maximum
    alphas = alphas + 0.01
    alphas[alphas > 1] = 1

    layout = go.Layout(height = 800)
    fig2 = go.Figure(layout=layout)
    fig2.add_trace(
        go.Scatter(
            x=start_points,
            y=end_points,
            mode="markers",
            text=[f"n_reads: {n}" for n in filtered_dbscan_df["n_reads"]],
            marker=dict(
                color=f"rgb{mcolors.to_rgb(color_dict[color_sample])}", opacity=alphas
            ),
        )
    )
    for template, start, end in zip(
        template_df["Fragment"], template_df["Start"], template_df["End"]
    ):
        fig2.add_shape(
            x0=start,
            x1=end,
            y0=start,
            y1=end,
            line=dict(color="red", width=1, dash="dash"),
        )

    fig2.update_layout(
        xaxis=dict(
            title="Start site", gridcolor="white", zerolinecolor="white", tickformat="d"
        ),
        yaxis=dict(
            title="End site",
            gridcolor="white",
            zeroline=True,
            zerolinecolor="white",
            zerolinewidth=0.1,
            tickformat="d",
        ),
        plot_bgcolor="rgba(0,0,0,0)",
    )

    #plotly.offline.plot(fig2, filename=f"{output}/intensity_matrix.html")
    plotly.io.write_html(fig2, f"{output}/intensity_matrix.html")


#####################################################################################################################################
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                             Script                                                                #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#####################################################################################################################################


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
        "Associated_Fragments_Overlap"
    ],
    has_header=True,
)
dbscan_matrix, id_dict = create_intensity_matrix(alignment_df)
plot_matrix(dbscan_matrix, color_sample, template_df_name, id_dict)
