import argparse
import os
from tqdm import tqdm
import math
import pysam
import pandas as pd
import numpy as np
import matplotlib as plt
import json
import logging
import sys
import random
import ctypes


class MutableString(object):
    def __init__(self, data):
        self.data = list(data)
    def __repr__(self):
        return "".join(self.data)
    def __setitem__(self, index, value):
        self.data[index] = value
    def __getitem__(self, index):
        if type(index) == slice:
            return "".join(self.data[index])
        return self.data[index]
    def __delitem__(self, index):
        del self.data[index]
    def __add__(self, other):
        self.data.extend(list(other))
    def __len__(self):
        return len(self.data)

opt_parser = argparse.ArgumentParser()

opt_parser.add_argument("-r", "--reference_file",dest="fasta_name", help="Insert a reference fasta file", metavar="FILE")
opt_parser.add_argument("-c", "--number_of_reads",dest="n_reads", help="Insert a reference fasta file", metavar="FILE")
opt_parser.add_argument("-o", "--output_path",dest="output", help="Insert an output directory to write to", metavar="FILE")

options = opt_parser.parse_args()


fasta_name = options.fasta_name
output = options.output
n_reads = int(options.n_reads)

if not os.path.exists(output):
    os.mkdir(output)


fasta_file = pysam.FastaFile(fasta_name)
reference = fasta_file.references[0]
reference_string = fasta_file.fetch(reference)
start_index = 0
end_index = len(reference_string)


n_extracted_reads = 0



cluster_list = [[3654,5524],[3654,5524],[3654,5524],[3654,5524],[3654,5524],[3654,5524],[3654,5524],\
                [6601,6756],[6601,6756],[6601,6756],[6601,6756],\
                [7925,12990],[7925,12990],[7925,12990],[7925,12990],[7925,12990],[7925,12990],[7925,12990],\
                [0,13349],[0,13349]]

cluster_index = random.randrange(0,len(cluster_list),1)
start = cluster_list[cluster_index][0]
end = cluster_list[cluster_index][1]

current_start_coordinate = random.randrange(start, end,1)
current_end_coordinate = random.randrange(start, end,1)

iterator = 0
while current_end_coordinate <= current_start_coordinate + 50:
    current_end_coordinate = random.randrange(current_start_coordinate, end,1)
    iterator += 1
    if iterator >= 10:
        current_start_coordinate = random.randrange(start, current_end_coordinate,1)

dictionairy_of_reads = {}
iteration = 1

with open(f"{output}/validation.fastq", 'w') as file:
    file.write("")

while n_extracted_reads < n_reads:
    fractionizer = [50000] * 25 + [25000] * 50 + [12000] * 100 + [6000] * 50 + [3000] * 25 + [100] * 10 
    fractionizer_index = random.randrange(0,len(fractionizer))
    current_number_of_reads = n_reads // fractionizer[fractionizer_index] 
    #print(current_number_of_reads)
    if (n_extracted_reads + current_number_of_reads) > n_reads:
        current_number_of_reads = n_reads - n_extracted_reads
    n_extracted_reads += current_number_of_reads
    for index in range(current_number_of_reads):
        cut_start = 2 #random.randrange(1,2,1)
        cut_end = 2 #random.randrange(1,2,1)
        current_string = str(reference_string[current_start_coordinate:current_end_coordinate])
        random_changes_indices = [random.randrange(0,len(current_string),1) for i in range(int(len(current_string)*0.05))]
        change_type = [random.randrange(1,4,1) for i in range(len(random_changes_indices))]
        #print(change_type)
        changes_on_string = [0 for i in range(len(current_string))]
        for random_index,change_type in zip(random_changes_indices,change_type):
            changes_on_string[random_index] = change_type
        #print(changes_on_string)
        changed_string = ""
        skip = False
        counter = 0
        to_skip = 0
        value_list = []
        for index2, value in enumerate(changes_on_string):
            value_list.append(value)
            if counter != 0:#Part of deletion
                if counter <= to_skip:
                    counter += 1
                    continue
                else:
                    counter = 0
                    continue
            if skip == True:
                continue
            if value == 0: #Normal
                changed_string = changed_string + current_string[index2]
                continue
            elif value == 1: #Substitution
                changed_string = changed_string + str(random.choice("ACGT"))
                continue
            elif value == 2:#Deletion
                skip == True
                to_skip = random.randrange(1,5,1) 
                continue
            elif value == 3:#Insertion
                for i in range(random.randrange(1,5,1)):
                    changed_string = changed_string + random.choice("ACGT")
                changed_string = changed_string + current_string[index2]
                continue
        if cut_start > 1:
            cut_size_start = random.randrange(1,max(min(len(changed_string)//5,200),2),1)
            changed_string = changed_string[cut_size_start:]
        if cut_end > 1:
            cut_size_end = random.randrange(1,max(min(len(changed_string)//5,200),2),1)
            changed_string = changed_string[:-cut_size_end]
        read_name = f"@iteration_{iteration}_read_{index}"
        read_direction = "+"
        read_sequence = changed_string
        read_quality = "I" * len(read_sequence)
        with open(f"{output}/validation.fastq", 'a') as file:
            file.write(f"{read_name}\n")
            file.write(f"{read_sequence}\n")
            file.write(f"{read_direction}\n")
            file.write(f"{read_quality}\n")
        
    dictionairy_of_reads[iteration] = [iteration,current_number_of_reads,current_string,current_start_coordinate,current_end_coordinate]
    cluster_index = random.randrange(0,len(cluster_list)-1,1)
    start = cluster_list[cluster_index][0]
    end = cluster_list[cluster_index][1]
    iterator = 0

    current_start_coordinate = random.randrange(start, end ,1)
    current_end_coordinate = random.randrange(start, end,1)
    
    while current_end_coordinate <= current_start_coordinate + 50:
        current_end_coordinate = random.randrange(current_start_coordinate, end,1)
        iterator += 1
        if iterator >= 10:
            current_start_coordinate = random.randrange(start, current_end_coordinate,1)
    iteration += 1

dataframe_of_reads = pd.DataFrame.from_dict(dictionairy_of_reads,orient='index',columns=["iteration","n_reads","sequence","start","end"])
dataframe_of_reads["length"] = dataframe_of_reads["end"] - dataframe_of_reads["start"]
dataframe_of_reads = dataframe_of_reads.loc[:,dataframe_of_reads.columns != "sequence"].sort_values(by=['start','length'], ascending=[True,False]).reset_index(drop = True)
dataframe_of_reads.to_csv(f"{output}validation_dataset.csv",index = 0,sep="\t")
dataframe_of_reads["reference"] = "RNA45SN1"
dataframe_of_reads["strand"] = "+"
dataframe_of_reads["score"] = (dataframe_of_reads["n_reads"] / sum(dataframe_of_reads["n_reads"])) * 1000000
bed_of_reads = dataframe_of_reads[["reference","start","end","iteration","score","strand","start","end"]]
with open(f"{output}validation_dataset.bed", 'w') as fp:
    fp.write(bed_of_reads.to_csv(sep="\t",header=False,index=False))

os.mkdir(f"{output}/basecalling_output/")
os.system(f"minimap2 -ax map-ont ~/NanoRibolyzer/references/RNA45SN1.fasta {output}validation.fastq | samtools sort | samtools view -hbS -F 3884 > {output}/basecalling_output/filtered.bam")
os.system(f"samtools index {output}/basecalling_output/filtered.bam")
os.system(f"bamToBed -i {output}/basecalling_output/filtered.bam > {output}validation.bed")