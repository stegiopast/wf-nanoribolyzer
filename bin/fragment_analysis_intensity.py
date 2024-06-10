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

options = opt_parser.parse_args()

bamfile_name = options.bamfile_name
fasta_name = options.fasta_name
output = options.output
identity = float(options.identity)
sample_type = str(options.sample_type)
cores = int(options.cores)

# Logger construction
logger = logging.getLogger(__name__)
stream = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(
    "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
stream.setFormatter(formatter)
if os.path.exists(f"{output}fragment_analysis.log"):
    with open(f"{output}fragment_analysis.log", "w") as file:
        file.write("")
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

"""
An alignment between Reference and read sequences is reconstructed in a pandas dataframe using the syntax documented in pysam.
"""


def reconstruct_alignment(samfile: pysam.AlignmentFile):
    counter_forward = 0
    counter_reverse = 0
    temp_list = []
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
        }
        temp_list.append(read_data)
    alignment_df = pl.DataFrame(temp_list)
    return alignment_df, counter_forward, counter_reverse


"""
Intensity matrix spans over reference_length x reference_length, which is initialized with 0. Start and Edn point of reads are used to increment the x,y position on the matrix 
to obtain an start-end site intensity matrix. Similarly the read_position is stored in an start,end dictionairy using the coordiantes as keys.
"""


def create_intensity_matrix(fusion_alignment_df: pl.DataFrame):
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
        [[x, y, intensity_matrix[x, y]] for x, y in zip(*np.nonzero(intensity_matrix))]
    )
    return dbscan_matrix, id_dict


"""
The intensity matrix is used to determine which reads have exactly the same start an end sites. Perfectly matching reads will be considered as fragment clusters. 
"""


def intensity_clustering(
    fusion_alignment_df: pl.DataFrame, dbscan_matrix: np.array, id_dict: dict
):
    dbscan_dataframe = pd.DataFrame(dbscan_matrix, columns=["Start", "End", "Counts"])
    single_read_ids = []
    list_read_groups = []
    for Start, End, Count in tqdm(
        zip(dbscan_dataframe.Start, dbscan_dataframe.End, dbscan_dataframe.Counts),
        total=len(dbscan_dataframe.Start),
    ):
        # This is necessary since unclustered reads can still be composed of many reads in the intensity matrix
        if Count <= 1:
            single_read_ids.append(id_dict[f"{Start}:{End}"][0])
            continue
        else:
            # read_group_df = fusion_alignment_df.lazy().filter(pl.col("Refstart") == Start).filter(pl.col("Refend") == End).collect()["ID"]
            read_group_list = id_dict[f"{Start}:{End}"]
            list_read_groups.append(read_group_list)
    single_reads_df = fusion_alignment_df.filter(pl.col("ID").is_in(single_read_ids))
    return list_read_groups, single_reads_df


"""
The output of the intensity_clustering function determines which reads can be considered to origin from the same fragment. 
The function performs an additive approach to determine a consensus sequence of all reads belonging to a common fragment. 
The final sequence is determined by the most common base at a given position from the minstart to the maxend of the determined reads.
"""


def fusion_read_groups(temp_df: pl.DataFrame, reference_dict: dict):
    reference_length = int(reference_dict["Length"])
    # temp_df = pl.read_json(_read_group_path)
    if not temp_df.is_empty():
        min_refstart = min(temp_df["Refstart"])
        max_refend = max(temp_df["Refend"])
        final_consensus_ids = [id for id in temp_df["ID"]]
        number_of_reads = sum(temp_df["n_Reads"])
        position_array = [
            {"A": 0, "T": 0, "C": 0, "G": 0, "U": 0, "-": 0}
            for location in range(min_refstart, max_refend + 1)
        ]
        for row_seq, rowstart, rowend in zip(
            temp_df["Sequence"], temp_df["Refstart"], temp_df["Refend"]
        ):
            for j in range(min_refstart, max_refend + 1):
                index_seq = j - rowstart
                index_position_array = j - min_refstart
                if j >= rowstart and j < rowend:
                    nucleotide = row_seq[index_seq]
                    position_array[index_position_array][nucleotide] += 1
                else:
                    position_array[index_position_array]["-"] += 1
        string = ""
        absolute_base_count_array = []
        for item in position_array:
            string += max(item.items(), key=operator.itemgetter(1))[0]
            # for item in position_array:
            item_conv_list = [0, 0, 0, 0, 0, 0]
            item_conv_list[0] = item["A"]
            item_conv_list[1] = item["T"]
            item_conv_list[2] = item["C"]
            item_conv_list[3] = item["G"]
            item_conv_list[4] = item["U"]
            item_conv_list[5] = item["-"]
            absolute_base_count_array.append(item_conv_list)
        temp_string = (
            "-" * min_refstart + string + "-" * (reference_length - max_refend)
        )
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


"""
Read (r) and reference (ref) sequence are compared base by base to verify the alignment quality
"""


def percentage_of_fits_parrallel(
    index: int,
    alignment_df_row: list,
    reference_length: int,
    reference_sequence: str,
    identity: float,
):
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
The function creates bed files with fragments, which can be visualized in igv. Different colorscheme were applied for different conditions
"""


def create_colored_bed(
    table_name: str = "",
    output_folder: str = output,
    output_file: str = "",
    sample_type: str = "",
):
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

    mean = df["rel_n_Reads"].mean()
    over_mean_df = df[df["rel_n_Reads"] > mean]
    bed_df = over_mean_df[
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
    with open(f"{output}over_mean_{output_file}", "w") as fp:
        fp.write(bed_df.to_csv(sep="\t", header=False, index=False))


"""
This function mediates the intensity matrix construction and consensus sequence determination processes.  
"""


def intensity_fusion(fusion_alignment_df=pl.DataFrame()):
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
    logger.info("Run intensity clustering")
    list_read_groups, single_reads_df = intensity_clustering(
        fusion_alignment_df, dbscan_matrix, id_dict
    )
    consensus_rows = []
    logger.info("Find consensus of defined clusters")
    # temp_fusion_alignment_df = fusion_alignment_df
    for id_list in tqdm(list_read_groups, total=len(list_read_groups)):
        out_dict = fusion_read_groups(
            fusion_alignment_df.lazy().filter(pl.col("ID").is_in(id_list)).collect(),
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
alignment_df = alignment_df.sort(["Refstart", "Length"], descending=[False, True])

logger.info("Perform intensity clustering")
final_df, single_reads_df = intensity_fusion(alignment_df)
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
