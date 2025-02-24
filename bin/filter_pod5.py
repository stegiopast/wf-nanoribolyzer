import pandas as pd
import os
import pysam
import shutil
import argparse
from tqdm import tqdm
import pod5 as p5

# Argument Parser
opt_parser = argparse.ArgumentParser()
opt_parser.add_argument(
    "-i",
    "--input_bam_file",
    dest="bamfile_name",
    help="Insert a sample bam file",
    metavar="FILE",
)
opt_parser.add_argument(
    "-o",
    "--output_pod5_directory",
    dest="pod5_output_dir_name",
    help="Insert an input pod5 directory with sample bam file",
    metavar="FILE",
)

options = opt_parser.parse_args()
bamfile_name = options.bamfile_name
destination_directory_pod5 = options.pod5_output_dir_name

# Initialize bamfile
samfile = pysam.AlignmentFile(bamfile_name, "rb")
alignment_ids = []

# Write filenames of successfully aligned bamfile into textfile
for index, read in tqdm(enumerate(samfile.fetch())):
    if index == 0:
        alignment_ids.append(read.query_name)
    else:
        alignment_ids.append(f"\n{read.query_name}")

with open(f"{destination_directory_pod5}/filtered_reads.txt", "w") as file:
    file.write("")

with open(f"{destination_directory_pod5}/filtered_reads.txt", "a") as file:
    file.writelines(alignment_ids)

# filter pod5 file to obtain pod5 file with reads only aligning to rRNA (45SN1 of hg38)
os.system(
    f"cat {destination_directory_pod5}/filtered_reads.txt | sort > {destination_directory_pod5}/sorted_filtered_reads.txt"
)
