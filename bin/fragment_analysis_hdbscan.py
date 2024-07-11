import argparse
import os
from tqdm import tqdm
import math
import pysam
import pandas as pd
import numpy as np
import operator
from multiprocessing import get_context
from multiprocessing import Pool
from itertools import repeat
import json
import logging
import sys
import polars as pl
from sklearn.cluster import HDBSCAN
import pyarrow as pa


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
    "-i",
    "--input_file",
    dest="bamfile_name",
    help="Insert a sample bam file",
    metavar="FILE",
)
opt_parser.add_argument(
    "-r",
    "--reference_file",
    dest="fasta_name",
    help="Insert a reference fasta file",
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
    "-m",
    "--min_neighbors",
    dest="min_neighbors",
    help="Define a minimal number of neighbors for HDBSCAN",
    metavar="FILE",
)
opt_parser.add_argument(
    "-t",
    "--identity_threshold",
    dest="identity",
    help="Define a minimal matching percentage between reads and reference",
    metavar="FILE",
)
opt_parser.add_argument(
    "-c", "--cores", dest="cores", help="Insert number of cores", metavar="FILE"
)
opt_parser.add_argument(
    "-s",
    "--sample_type",
    dest="sample_type",
    help="Mention which type of sample you have",
    metavar="FILE",
)

opt_parser.add_argument(
    "-d",
    "--demand",
    dest="demand",
    help="Use high or low demand algorithms ? (Results can slightly vary)",
    metavar="FILE",
)

options = opt_parser.parse_args()

bamfile_name = options.bamfile_name
fasta_name = options.fasta_name
output = options.output
min_neighbors = int(options.min_neighbors)
identity = float(options.identity)
sample_type = str(options.sample_type)
cores = int(options.cores)
demand = str(options.demand)

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


def reconstruct_alignment(samfile: pysam.AlignmentFile):
    counter_forward = 0
    counter_reverse = 0
    temp_list = []
    alignment_df = pl.DataFrame(
        schema={
            "ID": pl.String,
            "Sequence": pl.String,
            "Proportion_Sequence": pl.List(pl.List(pl.Int32)),
            "Length": pl.Int32,
            "Refstart": pl.Int32,
            "Refend": pl.Int32,
            "n_Reads": pl.Int64,
            "IDS": pl.List(pl.String),
            "alignment_probability": pl.Float64,
        }
    )
    for read in tqdm(samfile.fetch(), total=number_of_reads):
        if read.is_forward:
            counter_forward += 1
            read_sequence = read.get_forward_sequence()
        else:
            counter_reverse += 1
            read_sequence = read.query_sequence
        c_tuple_read = read.cigartuples
        string = ""
        position_read = 0
        for c_tuple in c_tuple_read:
            operation = c_tuple[0]
            length = c_tuple[1]
            if operation == 0:
                string += read_sequence[position_read : (position_read + length)]
                position_read += length
            elif operation == 1:
                position_read += length
            elif operation == 2:
                string += "-" * length
            elif operation == 3:
                string += read_sequence[position_read] * length
            elif operation == 4:
                position_read += length
        read_data = {
            "ID": read.query_name,
            "Sequence": string,
            "Proportion_Sequence": [],
            "Length": int(read.reference_end - read.reference_start),
            "Refstart": read.reference_start,
            "Refend": read.reference_end,
            "n_Reads": 1,
            "IDS": [read.query_name],
            "alignment_probability": 1-(10**-(read.mapping_quality/10))
        }
        temp_list.append(read_data)
        if len(temp_list) >= 1000000:
            chunk_alignment_df = pl.DataFrame(
                temp_list,
                schema={
                    "ID": pl.String,
                    "Sequence": pl.String,
                    "Proportion_Sequence": pl.List(pl.List(pl.Int32)),
                    "Length": pl.Int32,
                    "Refstart": pl.Int32,
                    "Refend": pl.Int32,
                    "n_Reads": pl.Int64,
                    "IDS": pl.List(pl.String),
                    "alignment_probability": pl.Float64,
                },
            )
            alignment_df = pl.concat([alignment_df, chunk_alignment_df])
            temp_list = []
    chunk_alignment_df = pl.DataFrame(
        temp_list,
        schema={
            "ID": pl.String,
            "Sequence": pl.String,
            "Proportion_Sequence": pl.List(pl.List(pl.Int32)),
            "Length": pl.Int32,
            "Refstart": pl.Int32,
            "Refend": pl.Int32,
            "n_Reads": pl.Int64,
            "IDS": pl.List(pl.String),
            "alignment_probability": pl.Float64,
        },
    )
    alignment_df = pl.concat([alignment_df, chunk_alignment_df])
    temp_list = []
    return alignment_df, counter_forward, counter_reverse


def create_intensity_matrix(fusion_alignment_df: pl.DataFrame):
    """
    Creates an intensity matrix from a DataFrame of fusion alignments and returns a matrix for DBSCAN clustering.

    Parameters:
    -----------
    fusion_alignment_df : pl.DataFrame
        A polars DataFrame containing fusion alignment information. The DataFrame must include the following columns:
        - "Refstart" : The reference start position of the alignment.
        - "Refend" : The reference end position of the alignment.
        - "ID" : The read/query name.

    Returns:
    --------
    tuple
        A tuple containing:
        - dbscan_matrix: np.ndarray, a 2D array of coordinates with non-zero intensity values. This array is suitable
          for use with DBSCAN clustering.
        - id_dict: dict, a dictionary mapping "Refstart:Refend" string keys to lists of read/query IDs that correspond
          to those start and end positions.

    The function constructs a 2D intensity matrix where the value at position (x, y) represents the number of alignments
    that start at position x and end at position y. It also creates a dictionary that maps each unique "Refstart:Refend"
    combination to a list of IDs that have those start and end positions.
    
    Example:
    --------
    >>> fusion_alignment_df = pl.DataFrame({
            "Refstart": [1, 2, 3],
            "Refend": [5, 6, 7],
            "ID": ["read1", "read2", "read3"]
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
    dbscan_matrix = np.array([[x, y] for x, y in zip(*np.nonzero(intensity_matrix))])
    return dbscan_matrix, id_dict


def hdbscan_clustering(
    fusion_alignment_df: pl.DataFrame,
    dbscan_matrix: np.array,
    id_dict: dict,
    min_number_of_neighbours: int,
):
    """
    Performs HDBSCAN clustering on fusion alignment data and returns grouped read IDs and a DataFrame of single reads.

    Parameters:
    -----------
    fusion_alignment_df : pl.DataFrame
        A polars DataFrame containing fusion alignment information.
    dbscan_matrix : np.array
        A 2D numpy array of coordinates with non-zero intensity values, used for clustering.
    id_dict : dict
        A dictionary mapping "Refstart:Refend" string keys to lists of read/query IDs that correspond to those start and end positions.
    min_number_of_neighbours : int
        The minimum number of neighboring points required for HDBSCAN to consider a point as part of a cluster.

    Returns:
    --------
    tuple
        A tuple containing:
        - list_read_groups: list, a list of lists where each inner list contains read IDs grouped into a cluster.
        - single_reads_df: pl.DataFrame, a DataFrame containing reads that were not clustered.

    The function performs the following steps:
    1. Applies HDBSCAN clustering on the provided `dbscan_matrix`.
    2. Constructs a DataFrame from the clustering results, associating each coordinate with a cluster label.
    3. Uses the labels to group read IDs into clusters or mark them as single reads.
    4. Compiles a list of grouped read IDs and a DataFrame of single reads.
    
    Example:
    --------
    >>> fusion_alignment_df = pl.DataFrame({
            "Refstart": [1, 2, 3],
            "Refend": [5, 6, 7],
            "ID": ["read1", "read2", "read3"]
        })
    >>> dbscan_matrix, id_dict = create_intensity_matrix(fusion_alignment_df)
    >>> list_read_groups, single_reads_df = hdbscan_clustering(fusion_alignment_df, dbscan_matrix, id_dict, 5)
    """
    result = HDBSCAN(min_cluster_size=min_number_of_neighbours).fit(dbscan_matrix)
    dbscan_dataframe = pd.DataFrame(dbscan_matrix, columns=["Start", "End"])
    dbscan_dataframe["Labels"] = np.array(result.labels_).astype(int)
    cluster_dict = {}
    additional_label = max(result.labels_) + 1
    single_read_ids = []
    list_read_groups = []
    for Start, End, Label in tqdm(
        zip(dbscan_dataframe.Start, dbscan_dataframe.End, dbscan_dataframe.Labels),
        total=len(dbscan_dataframe.Start),
    ):
        # This is necessary since unclustered reads can still be composed of many reads in the intensity matrix
        if len(id_dict[f"{Start}:{End}"]) <= 1:
            if Label == -1 or Label == "-1":
                id = id_dict[f"{Start}:{End}"][0]
                single_read_ids.append(id)
                continue
        if Label == -1 or Label == "-1":
            Label = additional_label
            additional_label += 1
        if Label not in cluster_dict.keys():
            cluster_dict[Label] = []
        for id in id_dict[f"{Start}:{End}"]:
            cluster_dict[Label].append(id)
    for cluster_number in tqdm(cluster_dict.keys(), total=len(cluster_dict.keys())):
        cluster_list = list(cluster_dict[cluster_number])
        list_read_groups.append(cluster_list)
    single_reads_df = (
        fusion_alignment_df.lazy().filter(pl.col("ID").is_in(single_read_ids)).collect()
    )
    return list_read_groups, single_reads_df


def fusion_read_groups_low_demand(temp_df: pl.DataFrame, reference_dict: dict):
    if not temp_df.is_empty():
        mean_refstart = int(min(temp_df["Refstart"]))
        mean_refend = int(max(temp_df["Refend"]))
        final_consensus_ids = [id for id in temp_df["ID"]]
        number_of_reads = sum(temp_df["n_Reads"])
        string = reference_dict["Sequence"][mean_refstart:mean_refend]
        min_max_length = mean_refend - mean_refstart
        absolute_base_count_array = []
        output_dict = {
            "ID": f"{mean_refstart}:{mean_refend}:{number_of_reads}",
            "Sequence": string,
            "Proportion_Sequence": absolute_base_count_array,
            "Length": min_max_length,
            "Refstart": mean_refstart,
            "Refend": mean_refend,
            "n_Reads": number_of_reads,
            "IDS": final_consensus_ids,
        }
        return output_dict


def fusion_read_groups(temp_df: pl.DataFrame, reference_dict: dict):
    """
    Generates a consensus sequence and related information from a group of fusion reads.

    Parameters:
    -----------
    temp_df : pl.DataFrame
        A polars DataFrame containing the fusion read information. The DataFrame must include the following columns:
        - "ID": Read/query name.
        - "Sequence": Nucleotide sequence of the read.
        - "Refstart": Start position of the read on the reference.
        - "Refend": End position of the read on the reference.
        - "n_Reads": Number of reads contributing to this entry.
    reference_dict : dict
        A dictionary containing reference information, with a key "Length" indicating the length of the reference sequence.

    Returns:
    --------
    dict
        A dictionary with the following keys:
        - "ID": A string representing the range of the consensus sequence and the number of reads, formatted as "min_refstart:max_refend:number_of_reads".
        - "Sequence": The consensus nucleotide sequence.
        - "Proportion_Sequence": A list of lists representing the absolute base count at each position in the consensus sequence.
        - "Length": The length of the consensus sequence.
        - "Refstart": The start position of the consensus sequence.
        - "Refend": The end position of the consensus sequence.
        - "n_Reads": The total number of reads contributing to the consensus sequence.
        - "IDS": A list of read/query IDs contributing to the consensus sequence.

    The function performs the following steps:
    1. Determines the minimum start and maximum end positions from the provided DataFrame.
    2. Initializes a position array to count nucleotide occurrences at each position.
    3. Iterates through each read in the DataFrame, updating the position array based on the nucleotide counts.
    4. Constructs the consensus sequence by selecting the most frequent nucleotide at each position.
    5. Adjusts the consensus sequence to fit within the reference length, padding with gaps as needed.
    6. Calculates the absolute base counts for each position in the consensus sequence.
    7. Determines the true start and end positions of the consensus sequence by trimming leading and trailing gaps.
    8. Returns a dictionary containing the consensus sequence and related information.
    
    Example:
    --------
    >>> temp_df = pl.DataFrame({
            "ID": ["read1", "read2"],
            "Sequence": ["ATCG", "ATGG"],
            "Refstart": [0, 1],
            "Refend": [4, 5],
            "n_Reads": [1, 1]
        })
    >>> reference_dict = {"Length": 6}
    >>> output_dict = fusion_read_groups(temp_df, reference_dict)
    """
    reference_length = int(reference_dict["Length"])
    if not temp_df.is_empty():
        min_refstart = min(temp_df["Refstart"])
        max_refend = max(temp_df["Refend"])
        final_consensus_ids = [id for id in temp_df["ID"]]
        number_of_reads = sum(temp_df["n_Reads"])
        position_array = [
            {"A": 0, "T": 0, "C": 0, "G": 0, "U": 0, "-": 0}
            for location in range(min_refstart, max_refend + 1)
        ]
        for row_seq, rowstart, rowend, row_n_reads in zip(
            temp_df["Sequence"],
            temp_df["Refstart"],
            temp_df["Refend"],
            temp_df["n_Reads"],
        ):
            for j in range(min_refstart, max_refend + 1):
                index_seq = j - rowstart
                index_position_array = j - min_refstart
                if j >= rowstart and j < rowend:
                    nucleotide = row_seq[index_seq]
                    position_array[index_position_array][nucleotide] = (
                        position_array[index_position_array][nucleotide] + row_n_reads
                    )
                else:
                    position_array[index_position_array]["-"] = (
                        position_array[index_position_array]["-"] + row_n_reads
                    )
        string = ""
        for k in position_array:
            string += max(k.items(), key=operator.itemgetter(1))[0]
        temp_string = (
            "-" * min_refstart + string + "-" * (reference_length - max_refend)
        )
        absolute_base_count_array = []
        for item in position_array:
            item_conv_list = [0, 0, 0, 0, 0, 0]
            item_conv_list[0] = item["A"]
            item_conv_list[1] = item["T"]
            item_conv_list[2] = item["C"]
            item_conv_list[3] = item["G"]
            item_conv_list[4] = item["U"]
            item_conv_list[5] = item["-"]
            absolute_base_count_array.append(item_conv_list)
        min_refstart = len(temp_string)
        max_refend = 0
        for j, val in enumerate(temp_string):
            if val != "-":
                if j < min_refstart:
                    min_refstart = j
                if j > max_refend:
                    max_refend = j
        min_max_length = max_refend - min_refstart
        string = temp_string[min_refstart:max_refend]
        output_dict = {
            "ID": f"{min_refstart}:{max_refend}:{number_of_reads}",
            "Sequence": string,
            "Proportion_Sequence": absolute_base_count_array,
            "Length": min_max_length,
            "Refstart": min_refstart,
            "Refend": max_refend,
            "n_Reads": number_of_reads,
            "IDS": final_consensus_ids,
        }
        return output_dict



def percentage_of_fits_parrallel(
    index: int,
    alignment_df_row: list,
    reference_length: int,
    reference_sequence: str,
    identity: float,
):
    """
    Calculates the percentage of bases in an alignment row that match the reference sequence and checks if it meets a given identity threshold.

    Parameters:
    -----------
    index : int
        The index of the current alignment row being processed.
    alignment_df_row : list
        A list containing the sequence information for the alignment row. The expected format is [sequence, ref_start, ref_end].
    reference_length : int
        The length of the reference sequence.
    reference_sequence : str
        The reference nucleotide sequence against which the alignment is compared.
    identity : float
        The identity threshold for the proportion of matching bases required to consider the alignment as fitting.

    Returns:
    --------
    tuple
        A tuple containing:
        - int: The index of the alignment row if the percentage of fitting bases meets or exceeds the identity threshold, otherwise -1.
        - float: The percentage of fitting bases if it meets the identity threshold, otherwise -1.

    The function performs the following steps:
    1. Constructs the full sequence of the alignment row, including padding with gaps to match the reference length.
    2. Iterates through the reference sequence and the alignment row sequence to count the total bases and the number of matching bases.
    3. Calculates the percentage of matching bases.
    4. Checks if the percentage of matching bases meets the identity threshold.
    5. Returns the index and percentage of matching bases if the threshold is met, otherwise returns -1 for both values.
    
    Example:
    --------
    >>> index = 0
    >>> alignment_df_row = ["ATCG", 2, 6]
    >>> reference_length = 10
    >>> reference_sequence = "NNATCGNNNN"
    >>> identity = 0.75
    >>> percentage_of_fits_parrallel(index, alignment_df_row, reference_length, reference_sequence, identity)
    (0, 1.0)
    """
    row_sequence = (
        "-" * alignment_df_row[1]
        + alignment_df_row[0]
        + "-" * (reference_length - alignment_df_row[2])
    )
    n_of_bases = 0
    n_fitting_bases_to_reference = 0

    for i, val in enumerate(reference_sequence):
        if row_sequence[i] in ["A", "C", "G", "T", "U"]:
            n_of_bases += 1
            if row_sequence[i] == val:
                n_fitting_bases_to_reference += 1
    if n_of_bases == 0:
        return -1, -1
    percentage_of_fitting_bases = n_fitting_bases_to_reference / n_of_bases
    if percentage_of_fitting_bases >= identity:
        return index, percentage_of_fitting_bases
    else:
        return -1, -1


def create_colored_bed(
    table_name: str = "",
    output_folder: str = output,
    output_file: str = "",
    sample_type: str = "",
):
    """
    Generates BED files with colored annotations based on read counts from a CSV file.

    Parameters:
    -----------
    table_name : str
        The path to the input CSV file containing genomic data.
    output_folder : str
        The folder where the output BED files will be saved.
    output_file : str
        The name of the output BED file.
    sample_type : str
        The sample type which determines the color scheme to be applied. 
        Acceptable values are "Nucleus", "blue", "Cytoplasm", "red", "SN1", "green", 
        "SN2", "orange", "SN3", and "purple".

    The function reads the input CSV file, assigns colors to each row based on the 
    `rel_n_Reads` value and the specified sample type, and writes two BED files:
    1. A full BED file with all data.
    2. A filtered BED file containing only rows with `rel_n_Reads` greater than 0.001.

    The color assignment follows a predefined scheme for each sample type:
    - "Nucleus" / "blue" : various shades of blue.
    - "Cytoplasm" / "red" : various shades of red.
    - "SN1" / "green" : various shades of green.
    - "SN2" / "orange" : various shades of orange.
    - "SN3" / "purple" : various shades of purple.
    """
    df = pd.read_csv(table_name, sep=";", header=0)
    colors_list = []
    if sample_type == "Nucleus" or sample_type == "blue":
        for i, row in df.iterrows():
            if row["rel_n_Reads"] > 0.05:
                colors_list.append("4,65,127")
            elif row["rel_n_Reads"] > 0.01:
                colors_list.append("38,107,176")
            elif row["rel_n_Reads"] > 0.005:
                colors_list.append("79,144,209")
            elif row["rel_n_Reads"] > 0.001:
                colors_list.append("153,204,255")
            elif row["rel_n_Reads"] > 0.0005:
                colors_list.append("144,196,249")
            else:
                colors_list.append("190,222,255")
    elif sample_type == "Cytoplasm" or sample_type == "red":
        for i, row in df.iterrows():
            if row["rel_n_Reads"] > 0.05:
                colors_list.append("152,19,19")
            elif row["rel_n_Reads"] > 0.01:
                colors_list.append("198,51,51")
            elif row["rel_n_Reads"] > 0.005:
                colors_list.append("228,86,86")
            elif row["rel_n_Reads"] > 0.001:
                colors_list.append("230,100,100")
            elif row["rel_n_Reads"] > 0.0005:
                colors_list.append("255,129,129")
            else:
                colors_list.append("255,204,204")
    elif sample_type == "SN1" or sample_type == "green":
        for i, row in df.iterrows():
            if row["rel_n_Reads"] > 0.05:
                colors_list.append("0,51,0")
            elif row["rel_n_Reads"] > 0.01:
                colors_list.append("0,102,0")
            elif row["rel_n_Reads"] > 0.005:
                colors_list.append("0,204,0")
            elif row["rel_n_Reads"] > 0.001:
                colors_list.append("0,255,0")
            elif row["rel_n_Reads"] > 0.0005:
                colors_list.append("153,255,153")
            else:
                colors_list.append("204,255,204")
    elif sample_type == "SN2" or sample_type == "orange":
        for i, row in df.iterrows():
            if row["rel_n_Reads"] > 0.05:
                colors_list.append("51,25,0")
            elif row["rel_n_Reads"] > 0.01:
                colors_list.append("153,76,0")
            elif row["rel_n_Reads"] > 0.005:
                colors_list.append("204,102,0")
            elif row["rel_n_Reads"] > 0.001:
                colors_list.append("255,153,51")
            elif row["rel_n_Reads"] > 0.0005:
                colors_list.append("255,204,153")
            else:
                colors_list.append("255,229,204")
    elif sample_type == "SN3" or sample_type == "purple":
        for i, row in df.iterrows():
            if row["rel_n_Reads"] > 0.05:
                colors_list.append("25,0,51")
            elif row["rel_n_Reads"] > 0.01:
                colors_list.append("76,0,153")
            elif row["rel_n_Reads"] > 0.005:
                colors_list.append("127,0,255")
            elif row["rel_n_Reads"] > 0.001:
                colors_list.append("178,102,255")
            elif row["rel_n_Reads"] > 0.0005:
                colors_list.append("204,153,255")
            else:
                colors_list.append("229,204,255")
    df["Colors"] = colors_list
    df["Zeros"] = 0
    df["Strand"] = "."
    df["Thikstart"] = df["Refstart"]
    df["Thikend"] = df["Refend"]
    df["Score"] = df["rel_n_Reads"] * 1000000
    bed_df = df[
        [
            "ID",
            "Refstart",
            "Refend",
            "ID",
            "Score",
            "Strand",
            "Thikstart",
            "Thikend",
            "Colors",
        ]
    ]
    bed_df.iloc[:, 0] = "RNA45SN1"
    with open(f"{output}{output_file}", "w") as fp:
        fp.write(bed_df.to_csv(sep="\t", header=False, index=False))

    most_abundant_df = df[df["rel_n_Reads"] > 0.001]
    bed_df = most_abundant_df[
        [
            "ID",
            "Refstart",
            "Refend",
            "ID",
            "Score",
            "Strand",
            "Thikstart",
            "Thikend",
            "Colors",
        ]
    ]
    bed_df.iloc[:, 0] = "RNA45SN1"
    with open(f"{output}most_abundant_{output_file}", "w") as fp:
        fp.write(bed_df.to_csv(sep="\t", header=False, index=False))


def hdbscan_fusion(
    fusion_alignment_df: pl.DataFrame = pl.DataFrame(),
    min_number_of_neighbours: int = 5,
    demand: str = "low",
):
    """
    Performs hierarchical density-based spatial clustering of applications with noise (HDBSCAN) on fusion alignment data to identify clusters of reads,
    and generates consensus sequences for each cluster.

    Parameters:
    -----------
    fusion_alignment_df : pl.DataFrame, optional
        A DataFrame containing fusion alignment data, by default an empty DataFrame.
    min_number_of_neighbours : int, optional
        The minimum number of neighbors for a point to be considered a core point in HDBSCAN, by default 5.

    Returns:
    --------
    consensus_df : pd.DataFrame
        A DataFrame containing the consensus sequences for each cluster, along with their metadata and the percentage of bases fitting the reference.
    single_reads_df : pl.DataFrame
        A DataFrame containing the reads that were not clustered.

    The function performs the following steps:
    1. Creates a reference dictionary with metadata about the reference sequence.
    2. Creates an intensity matrix and ID dictionary from the fusion alignment DataFrame.
    3. Performs HDBSCAN clustering on the intensity matrix.
    4. Generates consensus sequences for each cluster.
    5. Calculates the percentage of bases in the consensus sequences that match the reference sequence.
    6. Returns the consensus sequences DataFrame and the single reads DataFrame.
    
    Example:
    --------
    >>> fusion_alignment_df = pl.DataFrame({"ID": ["read1", "read2"], "Refstart": [0, 5], "Refend": [10, 15], "Sequence": ["ATCG", "CGTA"], "n_Reads": [1, 1]})
    >>> min_number_of_neighbours = 3
    >>> consensus_df, single_reads_df = hdbscan_fusion(fusion_alignment_df, min_number_of_neighbours)
    """
    reference_dict = {
        "ID": str(reference),
        "Sequence": str(fasta_file.fetch(reference)),
        "Proportion_Sequence": [],
        "Length": len(fasta_file.fetch(reference)),
        "Refstart": 0,
        "Refend": len(fasta_file.fetch(reference)),
        "n_Reads": 1,
        "IDS": [],
    }
    logger.info("Create intensity matrix")
    dbscan_matrix, id_dict = create_intensity_matrix(fusion_alignment_df)
    logger.info("Run hdbscan clustering")
    list_read_groups, single_reads_df = hdbscan_clustering(
        fusion_alignment_df, dbscan_matrix, id_dict, min_number_of_neighbours
    )
    consensus_rows = []
    logger.info("Find consensus of defined clusters")
    # temp_fusion_alignment_df = fusion_alignment_df
    if demand == "high":
        for id_list in tqdm(list_read_groups, total=len(list_read_groups)):
            out_dict = fusion_read_groups(
                fusion_alignment_df.lazy()
                .filter(pl.col("ID").is_in(id_list))
                .collect(),
                reference_dict,
            )
            # temp_fusion_alignment_df = temp_fusion_alignment_df.lazy().filter(~pl.col("ID").is_in(id_list)).collect()
            consensus_rows.append(out_dict)
    else:
        temp_fusion_alignment_df = fusion_alignment_df.select([
            "ID",
            "Length",
            "Refstart",
            "Refend",
            "n_Reads",
            "alignment_probability"
        ])
        for id_list in tqdm(list_read_groups, total=len(list_read_groups)):
            out_dict = fusion_read_groups_low_demand(
                temp_fusion_alignment_df.lazy()
                .filter(pl.col("ID").is_in(id_list))
                .collect(),
                reference_dict,
            )
            # temp_fusion_alignment_df = temp_fusion_alignment_df.lazy().filter(~pl.col("ID").is_in(id_list)).collect()
            consensus_rows.append(out_dict)
    consensus_df = pd.DataFrame.from_dict(consensus_rows)
    consensus_df = consensus_df.sort_values(
        by=["Refstart", "Length"], ascending=[True, False]
    ).reset_index(drop=True)

    percentage_of_fitting_bases = []
    for i, row_sequence, rowstart, rowend in tqdm(
        zip(
            range(consensus_df.shape[0]),
            consensus_df["Sequence"],
            consensus_df["Refstart"],
            consensus_df["Refend"],
        ),
        total=consensus_df.shape[0],
    ):
        row_sequence = (
            "-" * rowstart + row_sequence + "-" * (reference_dict["Length"] - rowend)
        )
        if i == 0:
            percentage_of_fitting_bases.append(1)
            continue
        n_of_bases = 0
        n_fitting_bases_to_reference = 0
        for index, val in enumerate(reference_dict["Sequence"]):
            if row_sequence[index] in ["A", "C", "G", "T", "U"]:
                n_of_bases += 1
                if row_sequence[index] == val:
                    n_fitting_bases_to_reference += 1
        if n_of_bases == 0:
            percentage_of_fitting_bases.append(0)
        else:
            percentage_of_fitting_bases.append(
                n_fitting_bases_to_reference / n_of_bases
            )
    consensus_df["Percentage_of_fits"] = percentage_of_fitting_bases
    return consensus_df, single_reads_df


#####################################################################################################################################
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                          Script                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#####################################################################################################################################

logger.info("Similarity driven analysis")
logger.info("Load Samfile")
samfile = pysam.AlignmentFile(bamfile_name, "rb")
number_of_reads = samfile.count()

logger.info("Load Reference in fasta format")
fasta_file = pysam.FastaFile(fasta_name)
reference = fasta_file.references[0]

# Create a list of dictionairies with reference data as first entry, the samples will follow in further steps
temp_list = list()
ref_data = {
    "ID": str(reference),
    "Sequence": str(fasta_file.fetch(reference)),
    "Proportion_Sequence": [],
    "Length": len(fasta_file.fetch(reference)),
    "Refstart": 0,
    "Refend": len(fasta_file.fetch(reference)),
    "n_Reads": 1,
    "IDS": [],
}


# Reconstruct alignment of each read to reference by using cigarstring of bam files
logger.info("Create alignment with sequence and cigarstring of bam file")
alignment_df, counter_forward, counter_reverse = reconstruct_alignment(samfile)


logger.info("Quality check for reference and query identity")
indices_list = [index for index in range(alignment_df.shape[0] - 1)]
intermediate_list = list(
    zip(alignment_df["Sequence"], alignment_df["Refstart"], alignment_df["Refend"])
)
read_list_df = [intermediate_list[i] for i in tqdm(indices_list)]
reference_len = len(fasta_file.fetch(reference))
reference_seq = str(fasta_file.fetch(reference))

# Alignment Quality control
if demand == "high":    
    with Pool(cores) as p:
        pool_output = p.starmap(
            percentage_of_fits_parrallel,
            zip(
                indices_list,
                read_list_df,
                repeat(reference_len),
                repeat(reference_seq),
                repeat(identity),
            ),
        )
    selected_indices = []
    for objects in pool_output:
        if objects[0] >= 0:
            selected_indices.append(objects[0])
    read_list_df = []
    p.close()

    alignment_df = alignment_df[selected_indices]
if demand == "low": 
    alignment_df = alignment_df.filter(pl.col("alignment_probability") >= identity)
alignment_df = alignment_df.sort(["Refstart", "Length"], descending=[False, True])

logger.info("Perform hdbscan clustering")
final_df, single_reads_df = hdbscan_fusion(alignment_df, min_neighbors, demand)
final_df["rel_n_Reads"] = final_df["n_Reads"] / number_of_reads
final_df.to_csv(f"{output}fragment_df.csv", sep=";")
total_number_of_clustered_reads = sum(final_df["n_Reads"])
total_number_of_unclustered_reads = single_reads_df.shape[0]
total_number_of_processed_reads = alignment_df.shape[0]
single_reads_df = single_reads_df.to_pandas(use_pyarrow_extension_array=True)
single_reads_df.to_csv(f"{output}unclustered_reads_df.csv", sep=";")

logger.info(
    "---                                  Write results                                     ---"
)
simple_consensus_df = pd.DataFrame(
    final_df[
        [
            "ID",
            "Refstart",
            "Refend",
            "Length",
            "n_Reads",
            "rel_n_Reads",
            "Percentage_of_fits",
        ]
    ]
)
simple_consensus_df.to_csv(f"{output}fragment_df_simple.csv", sep=";")
bed_consensus_df = final_df[["ID", "Refstart", "Refend", "ID"]]
bed_consensus_df.iloc[:, 0] = "RNA45SN1"
create_colored_bed(
    table_name=f"{output}fragment_df.csv",
    output_folder=output,
    output_file="no_template.bed",
    sample_type=sample_type,
)
alignment_df = alignment_df.to_pandas(use_pyarrow_extension_array=True)
alignment_df.to_csv(f"{output}alignment_df.csv", sep=";")

logger.info("Done")
logger.info("------------------")
logger.info(f"Metainformation")
logger.info(f"Number of reads: {number_of_reads}")
logger.info(f"Number of forward reads: {counter_forward}")
logger.info(f"Number of reverse reads: {counter_reverse}\n")
logger.info(f"Number of processed reads: {total_number_of_processed_reads}\n")
logger.info(f"Number of clustered reads: {total_number_of_clustered_reads}\n")
logger.info(f"Number of unclustered reads: {total_number_of_unclustered_reads}\n")
logger.info("Results after hdbscan clustering")
logger.info("Clusters with at least 0.01% of the reads")
logger.info(simple_consensus_df[simple_consensus_df.rel_n_Reads > 0.0001])
logger.info("Clusters with at least 0.1% of the reads")
logger.info(simple_consensus_df[simple_consensus_df.rel_n_Reads > 0.001])
logger.info("Clusters with at least 0.5% of the reads")
logger.info(simple_consensus_df[simple_consensus_df.rel_n_Reads > 0.005])
logger.info("Clusters with at least 1% of the reads")
logger.info(simple_consensus_df[simple_consensus_df.rel_n_Reads > 0.01])
logger.info("############################################\n")
