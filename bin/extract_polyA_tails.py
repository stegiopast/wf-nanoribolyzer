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
    "--output_directory",
    dest="output_dir_name",
    help="Insert an input pod5 directory with sample bam file",
    metavar="FILE",
)

#Set variables
options = opt_parser.parse_args()
bamfile_name = options.bamfile_name
output_dir = options.output_dir_name

# Read samfile
bamfile = pysam.AlignmentFile(bamfile_name, "rb", check_sq=False)

# Fetch taillength from bamfile
# Construct dataframe
list_read_data = []
list_read_ids = []
for read in bamfile.fetch(until_eof=True):
    read_id = read.query_name
    for i in read.get_tags():
        if i[0] == "pt":
            taillength = i[1]
    try:
        list_read_data.append([read_id, taillength])
    #    list_read_ids.append(read_id)
    except:
    #    list_read_ids.append(read_id)
        continue
    
#if len(list_read_data) == 0:
#    for read_id in list_read_ids:
#        list_read_data.append([read_id, 0])
    
df_read_data = pd.DataFrame(list_read_data, columns=["read_id", "taillength"])
df_read_data.index = df_read_data["read_id"]

# Write Dataframe
df_read_data.to_csv(f"{output_dir}/tail_estimation.csv", sep=";")
