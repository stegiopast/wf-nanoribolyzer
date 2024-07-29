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
from tqdm import tqdm
import os

# Argument Parser
opt_parser = argparse.ArgumentParser()

opt_parser.add_argument(
    "-b",
    "--bamfile",
    dest="bamfile_name",
    help="Insert a path to a bamfile",
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

bamfile_name = options.bamfile_name
output = options.output
reference_path = options.reference



def plot_reads_with_deletions(bamfile_name: str, minimal_deletion_size: int, length_of_reference: int):
    counter_forward = 0
    counter_reverse = 0
    
    coverage_counter_list = np.array([0 for nucleotide in range(length_of_reference)])
    coverage_index_list = [i + 1 for i in range(length_of_reference)]

    samfile = pysam.AlignmentFile(bamfile_name, "rb")
    number_of_reads = samfile.count()

    with pysam.AlignmentFile(f"{output}/deletions.bam", "wb", header=samfile.header) as outfile:
        for read in tqdm(samfile.fetch(), total=number_of_reads):
            if read.is_forward:
                counter_forward += 1
            else:
                counter_reverse += 1
            c_tuple_read = read.cigartuples
            for c_tuple in c_tuple_read:
                operation = c_tuple[0]
                length = c_tuple[1]
                if operation == 2 and length >= minimal_deletion_size:
                    outfile.write(read)
            start_position = read.reference_start
            end_position = read.reference_end
            coverage_counter_list[start_position:end_position] += 1

    deletion_bamfile_name = f"{output}/deletions.bam"
    os.system(f"samtools index {deletion_bamfile_name}")
    deletion_samfile = pysam.AlignmentFile(deletion_bamfile_name, "rb")
    coverage_counter_list_without_deletions = np.array([0 for nucleotide in range(length_of_reference)])
    coverage_counter_list_deletions = np.array([0 for nucleotide in range(length_of_reference)])
    counter_forward = 0
    counter_reverse = 0
    number_of_reads = deletion_samfile.count()
    for read in tqdm(deletion_samfile.fetch(), total=number_of_reads):    
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

        for nucleotide,position in zip(string,[index for index in range(read.reference_start,read.reference_end)]):
            if nucleotide != "-":
                coverage_counter_list_deletions[position] += 1
            coverage_counter_list_without_deletions[position] += 1
    
    print(coverage_counter_list)
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 10))
    # ax.plot(
    #     coverage_index_list,
    #     coverage_counter_list,
    #     color="black",
    #     linewidth=2,
    # )
    ax.plot(
        coverage_index_list,
        coverage_counter_list_without_deletions,
        color="blue",
        linewidth=2,
    )
    ax.plot(
        coverage_index_list,
        coverage_counter_list_deletions,
        color="red",
        linewidth=2,
        alpha = 0.5
    )
    #ax.set_ylim(0,1)
    # ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), title="Sample")
    ax.set_xlabel("Position on reference")
    plt.show()

                

                    




fasta_file = pysam.FastaFile(reference_path)
reference = fasta_file.references[0]
reference_sequence = str(fasta_file.fetch(reference))

plot_reads_with_deletions(bamfile_name,1000,len(reference_sequence))