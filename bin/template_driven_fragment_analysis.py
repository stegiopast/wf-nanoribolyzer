import argparse
from tqdm import tqdm
import argparse
import pysam
import pandas as pd
from multiprocessing import Pool
import numpy as np
from itertools import repeat
import logging
import sys
import polars as pl

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

# Argument parser
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
    "-f",
    "--fragment_file",
    dest="fragment_name",
    help="Insert a fragment reference file",
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
    "-c", "--cores", dest="cores", help="Insert number of cores", metavar="FILE"
)
opt_parser.add_argument(
    "-t",
    "--identity_threshold",
    dest="identity",
    help="Define a minimal matching percentage between reads and reference",
    metavar="FILE",
)
opt_parser.add_argument(
    "-s",
    "--sample_type",
    dest="sample_type",
    help="Mention which type of sample (Cytoplasm,Nucleus,SN1,SN2,SN3) you have or just name a colorset (blue,red,green,orange,purple)",
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
fragment_name = options.fragment_name
output = options.output
identity = float(options.identity)
sample_type = options.sample_type
cores = int(options.cores)
demand = str(options.demand)

# Logger construction
logger = logging.getLogger(__name__)
stream = logging.StreamHandler(sys.stdout)
file_stream = logging.FileHandler(f"{output}template_driven_analysis.log")
formatter = logging.Formatter(
    "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
stream.setFormatter(formatter)
logger.addHandler(stream)
logger.addHandler(file_stream)
file_stream.setFormatter(formatter)
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


"""
For each read (R) and each known fragment(F) an overlap in start-end and end-start direction is calculated.
Similarly the base matches between all R and all F are calculated. 
To calculate the overlap the mean of start-end and end-start direction overlap percentage is calculated.
The minimal overlap of both direction is taken in consdideration in the next step.
The argmax overlapping fragment is considered as best fitting fragment.
"""


def argmax_of_min_overlap_fragment_association_parallel(
    query_row_seq: str,
    query_row_refstart: int,
    query_row_refend: int,
    query_row_length: int,
    query_row_id: str,
    fragment_df: pd.DataFrame,
):
    """
    Identifies the fragment with the maximum overlap with the query sequence and calculates the fitting statistics.

    Parameters:
    -----------
    query_row_seq : str
        The sequence of the query row.
    query_row_refstart : int
        The reference start position of the query row.
    query_row_refend : int
        The reference end position of the query row.
    query_row_length : int
        The length of the query row sequence.
    query_row_id : str
        The ID of the query row.
    fragment_df : pd.DataFrame
        A DataFrame containing fragment information. Must include columns "Sequence", "Start", "End", "Length", and "Fragment".

    Returns:
    --------
    fitting_stats : list of lists
        A list where each element is a list containing query row IDs that best fit the corresponding fragment.
    fitting_stat_overlap : float
        The overlap ratio between the query sequence and the fragment with the maximum overlap.
    fitting_stat_matches : float
        The match ratio between the query sequence and the fragment with the maximum overlap.
    fitting_stat_fragment : Any
        The fragment that has the maximum overlap with the query sequence.
    query_row_id : str
        The ID of the query row.

    This function performs the following steps:
    1. Initializes lists to store fitting statistics.
    2. Constructs the full query sequence by padding with '*' characters based on reference start and end positions.
    3. Iterates over each fragment in the fragment DataFrame.
    4. Calculates the overlap and match ratios between the query sequence and each fragment.
    5. Determines the fragment with the maximum overlap.
    6. Returns the fitting statistics and the query row ID.

    Example:
    --------
    >>> fragment_df = pd.DataFrame({
    >>>     "Sequence": ["ATCG", "CGTA", "TACG"],
    >>>     "Start": [0, 5, 10],
    >>>     "End": [10, 15, 20],
    >>>     "Length": [10, 10, 10],
    >>>     "Fragment": ["frag1", "frag2", "frag3"]
    >>> })
    >>> result = argmax_of_min_overlap_fragment_association_parallel(
    >>>     "ATCG", 0, 10, 10, "query1", fragment_df
    >>> )
    """
    fitting_stats = [[] for i in range(0, len(fragment_df))]
    fitting_stat_overlap = []
    fitting_stat_matches = []
    fitting_stat_fragment = []
    fits = {"overlap": [], "matches": []}
    query_sequence = (
        "*" * query_row_refstart
        + query_row_seq
        + "*" * (ref_data["Length"] - query_row_refend)
    )
    for frag_row_seq, frag_row_start, frag_row_end, frag_row_length in zip(
        fragment_df["Sequence"],
        fragment_df["Start"],
        fragment_df["End"],
        fragment_df["Length"],
    ):
        frag_sequence = frag_row_seq
        temp_overlap_forward = 0
        temp_overlap_backward = 0
        temp_matches = 0
        max_start = max(query_row_refstart, frag_row_start)
        min_end = min(query_row_refend, frag_row_end)
        if min_end - max_start < 0:
            temp_overlap_choice = 0
            temp_matches = 0
        else:
            temp_overlap_forward = (min_end - max_start) / query_row_length
            temp_overlap_backward = (min_end - max_start) / frag_row_length
            temp_overlap_choice = min(temp_overlap_forward, temp_overlap_backward)
            for base_in_query, base_in_frag in zip(
                query_sequence[max_start:min_end], frag_sequence[max_start:min_end]
            ):
                if base_in_query == base_in_frag:
                    temp_matches += 1
            temp_matches = temp_matches / query_row_length
        fits["overlap"].append(temp_overlap_choice)
        fits["matches"].append(temp_matches)
    overlap = fits["overlap"][np.argmax(fits["overlap"])]
    matches = fits["matches"][np.argmax(fits["overlap"])]
    fitting_stats[np.argmax(fits["overlap"])].append(query_row_id)
    matching_fragment = fragment_df["Fragment"][np.argmax(fits["overlap"])]
    fitting_stat_overlap = overlap
    fitting_stat_matches = matches
    fitting_stat_fragment = matching_fragment
    return (
        fitting_stats,
        fitting_stat_overlap,
        fitting_stat_matches,
        fitting_stat_fragment,
        query_row_id,
    )


"""
The function creates bed files with fragments, which can be visualized in igv. Different colorscheme were applied for different conditions
"""


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
    df["Thikstart"] = df["Start"]
    df["Thikend"] = df["End"]
    df["Score"] = df["rel_n_Reads"] * 1000000

    bed_df = df[
        [
            "Fragment",
            "Start",
            "End",
            "Fragment",
            "Score",
            "Strand",
            "Thikstart",
            "Thikend",
            "Colors",
        ]
    ]
    bed_df.iloc[:, 0] = "RNA45SN1"
    template = 'track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"'
    with open(f"{output}{output_file}", "w") as fp:
        fp.write(bed_df.to_csv(sep="\t", header=False, index=False))

    most_abundant_df = df[df["rel_n_Reads"] > 0.001]
    bed_df = most_abundant_df[
        [
            "Fragment",
            "Start",
            "End",
            "Fragment",
            "Score",
            "Strand",
            "Thikstart",
            "Thikend",
            "Colors",
        ]
    ]
    bed_df.iloc[:, 0] = "RNA45SN1"
    template = 'track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"'
    with open(f"{output}most_abundant_{output_file}", "w") as fp:
        fp.write(bed_df.to_csv(sep="\t", header=False, index=False))


"""
General cut sites are determined by using the bamfile information of alignment start and end sites. Along the reference the start and end sites are quantified by addition on the reference position. 
The mean and standarddeviation between all positions with start- or end-sites is determined and end-sites that occur significantly often are considered to be true cut sites  
"""


def determine_general_cut_sites(alignment_df: pl.DataFrame, output: str):
    """
    Determines general cut sites based on alignment data and saves results in CSV and BED formats.

    Parameters:
    -----------
    alignment_df : pl.DataFrame
        DataFrame containing alignment information with columns "Refstart" and "Refend".
    output : str
        Path to the output directory where results will be saved.

    This function performs the following steps:
    1. Initializes lists to store start and end sites.
    2. Constructs a list of sequence bases from the reference.
    3. Iterates through each alignment row to count start and end sites.
    4. Computes mean and standard deviation of start and end sites.
    5. Determines thresholds for identifying significant start and end points.
    6. Flags sequences as starting or ending points based on thresholds.
    7. Saves the results in a CSV file named "cutting_sites_general.csv".
    8. Iterates through flagged start and end sites to create BED format data.
    9. Saves start sites in a BED file named "start_sites_general.bed".
    10. Saves end sites in a BED file named "end_sites_general.bed".

    Example:
    --------
    >>> alignment_data = pl.DataFrame({
    >>>     "Refstart": [10, 20, 30],
    >>>     "Refend": [15, 25, 35]
    >>> })
    >>> determine_general_cut_sites(alignment_data, "/path/to/output/")
    """
    start_sites = [0 for i in range(len(fasta_file.fetch(reference)))]
    end_sites = [0 for i in range(len(fasta_file.fetch(reference)))]
    sequence_list = [base for base in str(fasta_file.fetch(reference))]

    for row_refstart, row_refend in tqdm(
        zip(alignment_df["Refstart"], alignment_df["Refend"]),
        total=alignment_df.shape[0],
    ):
        start_sites[row_refstart - 1] += 1
        end_sites[row_refend - 1] += 1
    start_end_df = pd.DataFrame(
        {"Starts": start_sites, "Ends": end_sites, "Base": sequence_list},
        index=range(0, len(fasta_file.fetch(reference))),
    )

    starts_mean = start_end_df["Starts"].mean()
    starts_std = start_end_df["Starts"].std()
    starts_threshold = starts_mean + (2 * starts_std)
    ends_mean = start_end_df["Ends"].mean()
    ends_std = start_end_df["Ends"].std()
    ends_threshold = ends_mean + (2 * ends_std)
    start_end_df["Starting_points"] = [
        1 if i > starts_threshold else 0 for i in start_end_df["Starts"]
    ]
    start_end_df["End_points"] = [
        1 if i > ends_threshold else 0 for i in start_end_df["Ends"]
    ]
    start_end_df.to_csv(output + "cutting_sites_general.csv", sep=";")

    start_sites = []
    end_sites = []
    for index, row in tqdm(start_end_df.iterrows(), total=start_end_df.shape[0]):
        if row["Starting_points"] == 1:
            start_sites.append([index, row["Starts"]])
        if row["End_points"] == 1:
            end_sites.append([index, row["Ends"]])
    start_sites_list = []
    end_sites_list = []
    for i in tqdm(start_sites, total=len(start_sites)):
        temp_dict = {
            "Reference": "RNA45SN1",
            "Start": i[0],
            "End": i[0],
            "Fragment": "Startcut",
            "Score": i[1] / number_of_reads * 1000000,
            "Strand": ".",
            "Thikstart": i[0],
            "Thikend": i[0],
            "Colors": "0,153,0",
        }
        start_sites_list.append(temp_dict)
    start_sites_bed = pd.DataFrame.from_dict(start_sites_list)
    start_sites_bed.to_csv(
        output + "start_sites_general.bed", sep="\t", header=False, index=False
    )
    for j in tqdm(end_sites, total=len(end_sites)):
        temp_dict = {
            "Reference": "RNA45SN1",
            "Start": j[0],
            "End": j[0],
            "Fragment": "Endcut",
            "Score": j[1] / number_of_reads * 1000000,
            "Strand": ".",
            "Thikstart": j[0],
            "Thikend": j[0],
            "Colors": "0,153,153",
        }
        end_sites_list.append(temp_dict)
    end_sites_bed = pd.DataFrame.from_dict(end_sites_list)
    end_sites_bed.to_csv(
        output + "end_sites_general.bed", sep="\t", header=False, index=False
    )


"""
Fragment based cut sites are determined in a similar way like general cut sites, although the mean and stdd are calculated for each literature base fragment separately.
Intersecting fragment site values are determined by using the mean and the mean of the stdd between two the intersecting fragments.
"""


def determine_fragment_based_cut_sites(
    alignment_df: pl.DataFrame, fragment_df: pd.DataFrame, output: str
):
    """
    Determines fragment-based cut sites based on alignment and fragment data and saves results in CSV and BED formats.

    Parameters:
    -----------
    alignment_df : pl.DataFrame
        DataFrame containing alignment information with columns "Refstart" and "Refend".
    fragment_df : pd.DataFrame
        DataFrame containing fragment information with columns "Start", "End", "Fragment", "Length".
    output : str
        Path to the output directory where results will be saved.

    This function performs the following steps:
    1. Filters fragment DataFrame to include only specified fragments.
    2. Initializes lists to store start and end sites.
    3. Constructs a list of sequence bases from the reference.
    4. Iterates through each alignment row to count start and end sites.
    5. Calculates mean and standard deviation of start and end sites for each fragment.
    6. Computes thresholds for identifying significant start and end points based on fragment statistics.
    7. Flags sequences as starting or ending points based on computed thresholds.
    8. Saves the results in a CSV file named "cutting_sites_fragment_based.csv".
    9. Iterates through flagged start and end sites to create BED format data.
    10. Saves start sites in a BED file named "start_sites_fragment_based.bed".
    11. Saves end sites in a BED file named "end_sites_fragment_based.bed".

    Example:
    --------
    >>> alignment_data = pl.DataFrame({
    >>>     "Refstart": [10, 20, 30],
    >>>     "Refend": [15, 25, 35]
    >>> })
    >>> fragment_data = pd.DataFrame({
    >>>     "Start": [5, 10, 15],
    >>>     "End": [8, 13, 18],
    >>>     "Fragment": ["5ETS", "ITS1", "ITS2"],
    >>>     "Length": [4, 3, 4]
    >>> })
    >>> determine_fragment_based_cut_sites(alignment_data, fragment_data, "/path/to/output/")
    """
    cut_sides_df = fragment_df.loc[
        fragment_df["Fragment"].isin(
            ["5ETS", "ITS1", "ITS2", "3ETS", "18S", "5-8S", "28S"]
        )
    ]
    start_sites = [0 for i in range(len(fasta_file.fetch(reference)))]
    end_sites = [0 for i in range(len(fasta_file.fetch(reference)))]
    sequence_list = [base for base in str(fasta_file.fetch(reference))]

    for row_refstart, row_refend in tqdm(
        zip(alignment_df["Refstart"], alignment_df["Refend"]),
        total=alignment_df.shape[0],
    ):
        start_sites[row_refstart - 1] += 1
        end_sites[row_refend - 1] += 1
    start_end_df = pd.DataFrame(
        {"Starts": start_sites, "Ends": end_sites, "Base": sequence_list},
        index=range(0, len(fasta_file.fetch(reference))),
    )

    starts_means_over_reference = [[] for i in range(len(fasta_file.fetch(reference)))]
    starts_std_over_reference = [[] for i in range(len(fasta_file.fetch(reference)))]
    for index, row in cut_sides_df.iterrows():
        fragment_mean = start_end_df["Starts"][row["Start"] : row["End"]].mean()
        fragment_std = start_end_df["Starts"][row["Start"] : row["End"]].std()
        for index2 in range(row["Start"], row["End"]):
            starts_means_over_reference[index2].append(fragment_mean)
            starts_std_over_reference[index2].append(fragment_std)
    starts_means_over_reference_summ = [
        0 for i in range(len(fasta_file.fetch(reference)))
    ]
    starts_std_over_reference_summ = [
        0 for i in range(len(fasta_file.fetch(reference)))
    ]

    for index, element in enumerate(starts_means_over_reference):
        if len(element) > 1:
            starts_means_over_reference_summ[index] = np.array(element).mean()
        elif len(element) == 1:
            starts_means_over_reference_summ[index] = element[0]
        else:
            starts_means_over_reference_summ[index] = 0
    for index, element in enumerate(starts_std_over_reference):
        if len(element) > 1:
            starts_std_over_reference_summ[index] = np.array(element).mean()
        elif len(element) == 1:
            starts_std_over_reference_summ[index] = element[0]
        else:
            starts_std_over_reference_summ[index] = 0
    starts_threshold = [
        mean + (2 * std_dev)
        for mean, std_dev in zip(
            starts_means_over_reference_summ, starts_std_over_reference_summ
        )
    ]

    ends_means_over_reference = [[] for i in range(len(fasta_file.fetch(reference)))]
    ends_std_over_reference = [[] for i in range(len(fasta_file.fetch(reference)))]
    for index, row in cut_sides_df.iterrows():
        fragment_mean = start_end_df["Ends"][row["Start"] : row["End"]].mean()
        fragment_std = start_end_df["Ends"][row["Start"] : row["End"]].std()
        for index2 in range(row["Start"], row["End"]):
            ends_means_over_reference[index2].append(fragment_mean)
            ends_std_over_reference[index2].append(fragment_std)
    ends_means_over_reference_summ = [
        0 for i in range(len(fasta_file.fetch(reference)))
    ]
    ends_std_over_reference_summ = [0 for i in range(len(fasta_file.fetch(reference)))]

    for index, element in enumerate(ends_means_over_reference):
        if len(element) > 1:
            ends_means_over_reference_summ[index] = np.array(element).mean()
        elif len(element) == 1:
            ends_means_over_reference_summ[index] = element[0]
        else:
            ends_means_over_reference_summ[index] = 0
    for index, element in enumerate(ends_std_over_reference):
        if len(element) > 1:
            ends_std_over_reference_summ[index] = np.array(element).mean()
        elif len(element) == 1:
            ends_std_over_reference_summ[index] = element[0]
        else:
            ends_std_over_reference_summ[index] = 0
    ends_threshold = [
        mean + (2 * std_dev)
        for mean, std_dev in zip(
            ends_means_over_reference_summ, ends_std_over_reference_summ
        )
    ]

    start_end_df["Start_thresholds"] = starts_threshold
    start_end_df["End_thresholds"] = ends_threshold
    start_candidates = []
    for index, threshold in enumerate(starts_threshold):
        if start_end_df["Starts"][index] > threshold:
            start_candidates.append(1)
        else:
            start_candidates.append(0)
    start_end_df["Starting_points"] = start_candidates

    end_candidates = []
    for index, threshold in enumerate(ends_threshold):
        if start_end_df["Ends"][index] > threshold:
            end_candidates.append(1)
        else:
            end_candidates.append(0)
    start_end_df["End_points"] = end_candidates
    start_end_df.to_csv(output + "cutting_sites_fragment_based.csv", sep=";")

    start_sites = []
    end_sites = []
    for index, row in tqdm(start_end_df.iterrows(), total=start_end_df.shape[0]):
        if row["Starting_points"] == 1:
            start_sites.append([index, row["Starts"]])
        if row["End_points"] == 1:
            end_sites.append([index, row["Ends"]])

    start_sites_list = []
    end_sites_list = []
    for i in tqdm(start_sites, total=len(start_sites)):
        temp_dict = {
            "Reference": "RNA45SN1",
            "Start": i[0],
            "End": i[0],
            "Fragment": "Startcut",
            "Score": i[1] / number_of_reads * 1000000,
            "Strand": ".",
            "Thikstart": i[0],
            "Thikend": i[0],
            "Colors": "0,153,0",
        }
        start_sites_list.append(temp_dict)
    start_sites_bed = pd.DataFrame.from_dict(start_sites_list)
    start_sites_bed.to_csv(
        output + "start_sites_fragment_based.bed", sep="\t", header=False, index=False
    )

    for j in tqdm(end_sites, total=len(end_sites)):
        temp_dict = {
            "Reference": "RNA45SN1",
            "Start": j[0],
            "End": j[0],
            "Fragment": "Endcut",
            "Score": j[1] / number_of_reads * 1000000,
            "Strand": ".",
            "Thikstart": j[0],
            "Thikend": j[0],
            "Colors": "0,153,153",
        }
        end_sites_list.append(temp_dict)
    end_sites_bed = pd.DataFrame.from_dict(end_sites_list)
    end_sites_bed.to_csv(
        output + "end_sites_fragment_based.bed", sep="\t", header=False, index=False
    )


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


logger.info("Template driven analysis")
logger.info("Load Samfile")
samfile = pysam.AlignmentFile(bamfile_name, "rb")
number_of_reads = samfile.count()


logger.info("Load Reference in fasta format")
fasta_file = pysam.FastaFile(fasta_name)


logger.info("Template fragment dataframe")
fragment_df = pd.read_csv(fragment_name, sep=";", header=0, index_col=None)
fragment_df["Length"] = fragment_df["End"] - fragment_df["Start"]
fragment_df = fragment_df.sort_values(
    by=["Start", "Length"], ascending=[True, True]
).reset_index(drop=True)
reference = fasta_file.references[0]

# Create a list of dictionairies with reference data as first entry, the samples will follow in further steps
temp_list = list()
ref_data = {
    "ID": str(reference),
    "Sequence": str(fasta_file.fetch(reference)),
    "Length": len(fasta_file.fetch(reference)),
    "Refstart": 0,
    "Refend": len(fasta_file.fetch(reference)),
}

# Create a preliminary dictionairy for each read the Bam file contains
# Read Sequence is primed with dashes in legnth of reference string
# Length of the sequence is primed by reference_end - reference_start
# Create alignment dataframe will primed placeholders

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

    #Alignment quality checkup
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


logger.info("Reconstruct sequence for rRNA fragments known in Literature")
fragment_sequences = []
for index, row in fragment_df.iterrows():
    new_fragment = "-" * int(ref_data["Length"])
    new_string = ref_data["Sequence"][int(row["Start"]) : int(row["End"])]
    new_fragment = (
        new_fragment[: int(row["Start"])] + new_string + new_fragment[int(row["End"]) :]
    )
    fragment_sequences.append(new_fragment)
fragment_df["Sequence"] = fragment_sequences


logger.info("Associate querys to argmax of overlap with literature fragments")
with Pool(cores) as p:
    pool_output = p.starmap(
        argmax_of_min_overlap_fragment_association_parallel,
        zip(
            alignment_df["Sequence"],
            alignment_df["Refstart"],
            alignment_df["Refend"],
            alignment_df["Length"],
            alignment_df["ID"],
            repeat(fragment_df),
        ),
    )
fitting_stats = [[] for i in range(0, len(fragment_df))]
fitting_stats_overlap = []
fitting_stats_matches = []
fitting_stats_fragment = []
fitting_stats_ids = []
for objects in pool_output:
    for index, element in enumerate(objects[0]):
        if element != []:
            fitting_stats[index].append(element[0])
    fitting_stats_overlap.append(objects[1])
    fitting_stats_matches.append(objects[2])
    fitting_stats_fragment.append(objects[3])
    fitting_stats_ids.append(objects[4])
p.close()


alignment_df = alignment_df.hstack([pl.Series("Overlap", fitting_stats_overlap)])
alignment_df = alignment_df.hstack([pl.Series("Matches", fitting_stats_matches)])
alignment_df = alignment_df.hstack(
    [pl.Series("Associated_Fragments_Overlap", fitting_stats_fragment)]
)
alignment_df = alignment_df.filter(pl.col("ID") != "RNA45SN1")

logger.info("Calculate general start and end sites")
determine_general_cut_sites(alignment_df, output)
logger.info("Write results for general start and end sites")

logger.info("Calculate fragment based start and end sites")
determine_fragment_based_cut_sites(alignment_df, fragment_df, output)
logger.info("Write results for fragment based start and end sites")

logger.info(" Write Results of query association")
alignment_df = alignment_df.to_pandas(use_pyarrow_extension_array=True)
alignment_df.to_csv(f"{output}template_alignment_df.csv", sep=";")

fragment_df["IDS"] = fitting_stats
fragment_df["n_Reads"] = [len(i) for i in fitting_stats]
fragment_df["rel_n_Reads"] = fragment_df["n_Reads"] / number_of_reads
fragment_df.to_csv(f"{output}template_fragment_df.csv", sep=";")


over_zero_fragment_df = fragment_df.loc[fragment_df["n_Reads"] > 0]
bed_fragment_df = pd.DataFrame(
    over_zero_fragment_df[["Reference", "Start", "End", "Fragment"]]
)
create_colored_bed(
    table_name=f"{output}template_fragment_df.csv",
    output_folder=output,
    output_file="template_driven.bed",
    sample_type=sample_type,
)


logger.info("Metainformation")
logger.info(f"Number of reads: {number_of_reads}")
logger.info(f"Number of forward reads: {counter_forward}")
logger.info(f"Number of reverse reads: {counter_reverse}")
logger.info("Argmax overlap fragment associated queries")
logger.info(fragment_df.sort_values(by="Length", ascending=True).reset_index())
