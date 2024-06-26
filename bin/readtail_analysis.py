import pysam
from tqdm import tqdm
from pathlib import Path
import operator
import os
import pandas as pd
import json
import argparse
import logging
import sys

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

opt_parser.add_argument("-i", "--bamfile_name", dest="bamfile_name", help="Insert a sample bam file", metavar="FILE")
opt_parser.add_argument("-r", "--df_name",dest="df_name", help="Insert a template csv file", metavar="FILE")
opt_parser.add_argument("-o", "--output_path",dest="output_path", help="Insert an output directory to write to", metavar="FILE")

options = opt_parser.parse_args()


bamfile_name = options.bamfile_name
bamfile_name = Path(bamfile_name)
df_name = options.df_name
df_name = Path(df_name)
output_path = options.output_path
output_path = Path(output_path)

# Logger construction
logger = logging.getLogger(__name__)
stream = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(
    "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
stream.setFormatter(formatter)
if os.path.exists(f"{output_path}readtail_analysis.log"):
    with open(f"{output_path}readtail_analysis.log", "w") as file:
        file.write("")
file_stream = logging.FileHandler(f"{output_path}readtail_analysis.log")
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

def fragment_based_readtail_analysis(df_name,bamfile_name,output_path,logger):
    """
    Perform polyU analysis on fragments from a DataFrame based on reads in a BAM file.

    Parameters:
    -----------
    df_name : str
        Path to the DataFrame CSV file containing columns: "IDS", "Fragment", "n_Reads", "Start", "End".
    bamfile_name : str
        Path to the BAM file containing aligned reads.
    output_path : str
        Path to the directory where analysis results will be saved.

    This function performs the following steps:
    1. Creates a directory for polyU analysis results if it doesn't exist.
    2. Reads the DataFrame containing fragment information.
    3. Iterates through each fragment, extracting relevant reads from the BAM file.
    4. For each read in the fragment, performs polyU analysis:
       - Identifies the last 20 bases with their nucleotide frequencies.
       - Determines the polyU tail position and sequences.
       - Writes sequences and analysis results to respective files.
    5. Saves the last 20 bases nucleotide frequencies in a JSON file.
    6. Saves a consensus sequence of the most frequent base in the last 20 bases in a TXT file.

    Example:
    --------
    >>> fragment_based_readtail_analysis("fragment_data.csv", "aligned_reads.bam", "/path/to/output/")
    """
    logger.info("Read Samfile")
    bamfile = pysam.AlignmentFile(bamfile_name, "rb")
    df = pd.read_csv(df_name, sep=";")
    
    logger.info("Extract last 20 nucleotides of reads")
    for fragment_id_list,fragment,n_reads,start,end in tqdm(zip(df["IDS"],df["Fragment"],df["n_Reads"],df["Start"],df["End"]),total=len(df["Fragment"])):
        logger.info(f"Extraction from reads of {fragment}")
        last_20 = [{"A":0,"C":0,"G":0,"T":0} for i in range(20)]
        fragment_id_list = eval(fragment_id_list)
        fragment_id_dict = {}
        for fragment_id in fragment_id_list:
            fragment_id_dict[fragment_id] = fragment_id
        with open(f"{output_path}/{fragment}.fasta","w") as fasta_output:
            fasta_output.write("")
        with open(f"{output_path}/{fragment}_tail_comparison.fasta","w") as fasta_output:
            fasta_output.write("")
        for read in tqdm(bamfile.fetch("RNA45SN1",start=max(0,start-2000),end=min(13350,end+2000))):
            if read.query_name in fragment_id_dict:
                if read.is_forward:
                    read_sequence = read.get_forward_sequence()
                else:
                    read_sequence = read.query_sequence
                reversed_read_sequence = read_sequence[::-1]
                position = 0
                first_A_met = False
                nucleotide_counter = 0
                adenin_counter = 0
                for nuc_index,nucleotide in enumerate(reversed_read_sequence):
                    if first_A_met:
                        if nucleotide == "A":
                            nucleotide_counter += 1
                            adenin_counter += 1
                            position += 1
                        else:
                            nucleotide_counter += 1
                            position += 1
                            difference = position - nucleotide_counter
                            if adenin_counter/max(position-difference-1,1) <= 0.7:
                                #print(position)
                                difference = position - nucleotide_counter
                                while (reversed_read_sequence[position] != "A" or (adenin_counter/max(position - difference,1)) <= 0.95):
                                    if position == 0:
                                        break
                                    polyA_sequence = read_sequence[-(position):-(difference)]
                                    adenin_counter = 0
                                    for i in polyA_sequence:
                                        if i == "A":
                                            adenin_counter += 1
                                    position -= 1 
                                if position != 0: 
                                    position += 1
                                break
                    else:
                        if nuc_index + 2 == len(reversed_read_sequence):
                            position = 0
                            break
                        if nuc_index > 50:
                            position = 0
                            break
                        if nucleotide == "A" and reversed_read_sequence[nuc_index + 1] == "A" and reversed_read_sequence[nuc_index + 2] == "A":
                            nucleotide_counter += 1
                            adenin_counter += 1
                            first_A_met = True
                            position = nuc_index 
                cutoff_sequence = read_sequence[-(position + 20):-(position)]
                long_cutoff_sequence = read_sequence[-(position + 20):]
                for index,base in enumerate(cutoff_sequence):
                    last_20[index][base] += 1
                if len(cutoff_sequence) == 20:
                    with open(f"{output_path}/{fragment}.fasta","a") as fasta_output:
                        if cutoff_sequence != "":
                            fasta_output.write(f">{read.query_name}\n")
                            fasta_output.write(f"{cutoff_sequence}\n")
                with open(f"{output_path}/{fragment}_tail_comparison.fasta","a") as fasta_output:
                    if cutoff_sequence != "":
                        fasta_output.write(f">{read.query_name}\n")
                        fasta_output.write(f"{cutoff_sequence}\n")
                        fasta_output.write(f"{long_cutoff_sequence}\n")
                        fasta_output.write(f"{adenin_counter/max(position - difference,1)}\n")
        string = ""
        for dictionairy in last_20:
            string += max(dictionairy.items(), key=operator.itemgetter(1))[0]
        
        with open(f"{output_path}/{fragment}.json", "w") as outfile: 
            json.dump(last_20, outfile)
        
        with open(f"{output_path}/{fragment}.txt","w") as fasta_output:
                fasta_output.write(string)
    logger.info("Extraction finished")
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
logger.info("Readtail analysis")
logger.info("Start tail analysis")
fragment_based_readtail_analysis(df_name,bamfile_name,output_path,logger)
logger.info("Tail analysis finished")
