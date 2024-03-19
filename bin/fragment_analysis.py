import argparse
import os
from tqdm import tqdm
import math
import pysam
import pandas as pd
import numpy as np
import operator
import matplotlib as plt
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from multiprocessing import Pool
from itertools import repeat
import json
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

opt_parser.add_argument("-i", "--input_file", dest="bamfile_name", help="Insert a sample bam file", metavar="FILE")
opt_parser.add_argument("-r", "--reference_file",dest="fasta_name", help="Insert a reference fasta file", metavar="FILE")
opt_parser.add_argument("-o", "--output_path",dest="output", help="Insert an output directory to write to", metavar="FILE")
opt_parser.add_argument("-t","--threshold", dest= "threshold", help="Insert a overlap Threshold between 0-1 (1 -> 100 percent overlap)",metavar="FILE")
opt_parser.add_argument("-c","--cores", dest= "cores", help="Insert number of cores",metavar="FILE")
opt_parser.add_argument("-s", "--sample_type",dest="sample_type", help="Mention which type of sample you have", metavar="FILE")

options = opt_parser.parse_args()

bamfile_name = options.bamfile_name
fasta_name = options.fasta_name
output = options.output
overlap_threshold = float(options.threshold)
sample_type = str(options.sample_type)
cores = int(options.cores)

logger = logging.getLogger(__name__)
stream = logging.StreamHandler(sys.stdout)
file_stream = logging.FileHandler(f"{output}fragment_analysis.log")
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


def size_bin_selection(template_length, reduce_gt_1000 = 100, reduce_gt_100 = 50, reduce_gt_10 = 25):
    while template_length % 1000 != 0:
        template_length += 1
    list_of_lengths = []
    temporary_length = template_length
    while temporary_length > 5:
        if temporary_length > 1000:
            list_of_lengths.append([temporary_length - 1,(temporary_length - reduce_gt_1000)]) 
            temporary_length -= reduce_gt_1000
        elif temporary_length > 100:
            list_of_lengths.append([temporary_length - 1,(temporary_length - reduce_gt_100)])
            temporary_length -= reduce_gt_100
        else:
            list_of_lengths.append([temporary_length - 1,(temporary_length - reduce_gt_10)]) 
            temporary_length -= reduce_gt_10
    return list_of_lengths


def pearson_correlation(_temp_alignment_df_path,_reference_dict):
    read_group_path = _temp_alignment_df_path.replace(".json", "")
    temp_list_read_groups = []
    temp_list_read_singles = []
    temp_single_reads_df = pd.DataFrame()
    with open(_temp_alignment_df_path) as json_read:
        _temp_alignment_df = json.load(json_read)  
    if _temp_alignment_df != "[]":
        _temp_alignment_df = pd.DataFrame.from_dict(_temp_alignment_df)
        if _temp_alignment_df.shape[0] > 2:
            pearson_corr_matrix = [[0 for j in range(_temp_alignment_df.shape[0])] for i in range(len(_reference_dict["Sequence"])+1)]
            temp_iterator = 0
            for row, rowstart, rowend in zip(_temp_alignment_df["Sequence"],_temp_alignment_df["Refstart"],_temp_alignment_df["Refend"]):
                for j,val in enumerate(row):
                    j = j + rowstart
                    if val != "-" and j < rowend:
                        pearson_corr_matrix[j][temp_iterator] = 1 
                temp_iterator += 1
            pearson_df = pd.DataFrame(pearson_corr_matrix,columns=list(_temp_alignment_df.index))
            pearson_corr_result_df = pearson_df.corr()
            pearson_corr_result_df = pearson_corr_result_df.iloc[0:pearson_corr_result_df.shape[0]+1,0:pearson_corr_result_df.shape[1]+1]
            dissimilarity = (pearson_corr_result_df - 1) / (-2)
            Tree = linkage(squareform(dissimilarity), 'complete')
            threshold = 1 - overlap_threshold
            labels = fcluster(Tree, threshold, criterion='distance')
            temp_fusion_candidates_list = [[] for i in range(max(labels)+1)]
            for i,value in enumerate(labels):
                temp_fusion_candidates_list[value].append(i)
            non_fusion_candidates_list = [list(np.unique(np.array(i))) for i in temp_fusion_candidates_list if len(i) < 2]
            fusion_candidates_list = [list(np.unique(np.array(i))) for i in temp_fusion_candidates_list if len(i) > 1]
            
            
            for index,fusion_group_list in enumerate(fusion_candidates_list):
                fusion_df = _temp_alignment_df.iloc[fusion_group_list,:]
                fusion_json = fusion_df.to_dict()
                json_path = f"{read_group_path}_read_group_{index}.json"
                with open(json_path,"w") as json_read:
                    json.dump(fusion_json, json_read)
                temp_list_read_groups.append(json_path)
            for non_fusion_group_list in non_fusion_candidates_list:
                non_fusion_df = _temp_alignment_df.iloc[non_fusion_group_list,:]
                temp_list_read_singles.append(non_fusion_df)
                temp_single_reads_df = pd.concat([temp_single_reads_df,non_fusion_df],axis = 0)
        else:
            temp_single_reads_df = pd.concat([temp_single_reads_df,_temp_alignment_df], axis = 0)
        os.remove(_temp_alignment_df_path)
        return temp_list_read_groups,temp_list_read_singles,temp_single_reads_df


def fusion_read_groups(_read_group_path,_reference_dict):
    with open(_read_group_path,"r") as json_read:
        temp_df = json.load(json_read)
    #temp_df = [dict(i) for i in temp_df]
    if temp_df != "[]":
        temp_df = pd.DataFrame.from_dict(temp_df)
        #print(temp_df.columns)
        #print(temp_df["Proportion_Sequence"])
        min_refstart = min(temp_df["Refstart"]) #// len(temp_df["Refstart"]) 
        max_refend = max(temp_df["Refend"]) #// len(temp_df["Refend"])
        final_consensus_ids = []
        for list_id in temp_df["IDS"]:
            for single_id in list_id:
                final_consensus_ids.append(single_id)
        number_of_reads = sum(temp_df["n_Reads"])
        position_array = [{"A": 0, "T": 0,"C": 0,"G": 0,"U": 0,"-": 0} for location in range(_reference_dict["Length"])]
        for row_seq,row_prop_seq,rowstart,rowend,row_n_reads in zip(temp_df["Sequence"],temp_df["Proportion_Sequence"],temp_df["Refstart"],temp_df["Refend"],temp_df["n_Reads"]):
            if row_prop_seq == []:
                for j in range(_reference_dict["Length"]):
                    if  j >= rowstart and j < rowend:
                        nucleotide = row_seq[j-rowstart]
                        position_array[j][nucleotide] = position_array[j][nucleotide] + row_n_reads  
                    else:
                        position_array[j]["-"] = position_array[j]["-"] + row_n_reads
            else:
                for j in range(_reference_dict["Length"]):
                    if  j >= rowstart and j < rowend:
                        position_array[j]["A"] = position_array[j]["A"] + row_prop_seq[j-rowstart][0]
                        position_array[j]["T"] = position_array[j]["T"] + row_prop_seq[j-rowstart][1]
                        position_array[j]["C"] = position_array[j]["C"] + row_prop_seq[j-rowstart][2]
                        position_array[j]["G"] = position_array[j]["G"] + row_prop_seq[j-rowstart][3]
                        position_array[j]["U"] = position_array[j]["U"] + row_prop_seq[j-rowstart][4]
                        position_array[j]["-"] = position_array[j]["-"] + row_prop_seq[j-rowstart][5]
                    else:
                        position_array[j]["-"] = position_array[j]["-"] + row_n_reads
        string = ""
        for k in position_array[min_refstart:max_refend]:
            string += max(k.items(), key=operator.itemgetter(1))[0]
        temp_string = "-" * min_refstart + string + "-" * (_reference_dict["Length"] - max_refend)
        absolute_base_count_array = []
        for item in position_array:
            item_conv_list = [0,0,0,0,0,0]
            item_conv_list[0] = item["A"]
            item_conv_list[1] = item["T"]
            item_conv_list[2] = item["C"]
            item_conv_list[3] = item["G"]
            item_conv_list[4] = item["U"]
            item_conv_list[5] = item["-"]
            absolute_base_count_array.append(item_conv_list)
        min_refstart = len(temp_string)
        max_refend = 0
        for j,val in enumerate(temp_string):
            if val != "-":
                if j < min_refstart:
                    min_refstart = j
                if j > max_refend:
                    max_refend = j
        min_max_length = max_refend - min_refstart
        string = temp_string[min_refstart:max_refend]
        output_dict = {"ID" : temp_df["ID"][0], "Sequence" : string, "Proportion_Sequence": absolute_base_count_array[min_refstart:max_refend], "Length": min_max_length, "Refstart": min_refstart, "Refend": max_refend,"n_Reads": number_of_reads,"IDS": final_consensus_ids}
        output_df = pd.DataFrame()
        return output_dict, output_df


def percentage_of_fits_parrallel(index,alignment_df_row, reference_length, reference_sequence):
    row_sequence = "-" * alignment_df_row[1] + alignment_df_row[0] + "-" * (reference_length - alignment_df_row[2])
    n_of_bases = 0
    n_fitting_bases_to_reference = 0
    
    for i,val in enumerate(reference_sequence):
        if row_sequence[i] in ["A","C","G","T","U"]:
            n_of_bases += 1
            if row_sequence[i] == val:
                n_fitting_bases_to_reference += 1
    if n_of_bases == 0:
        return -1, -1
    percentage_of_fitting_bases = n_fitting_bases_to_reference / n_of_bases
    if percentage_of_fitting_bases >= 0.9:
        return index, percentage_of_fitting_bases
    else:
        return -1, -1
    




def create_colored_bed(table_name = "", output_folder = output, output_file = "", sample_type = ""):
    df = pd.read_csv(table_name, sep=";",header=0)

    colors_list = []
    if sample_type == "Nucleus" or sample_type == "blue": 
        for i,row in df.iterrows():
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
        for i,row in df.iterrows():
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
        for i,row in df.iterrows():
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
        for i,row in df.iterrows():
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
        for i,row in df.iterrows():
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
    print(df)

    bed_df = df[["ID","Refstart","Refend","ID","Score","Strand","Thikstart","Thikend","Colors"]]
    bed_df.iloc[:,0] = "RNA45SN1"
    template = 'track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"'

    with open(output + output_file, 'w') as fp:
        fp.write(bed_df.to_csv(sep="\t",header=False,index=False))
    print(bed_df) 


    most_abundant_df = df[df["rel_n_Reads"] > 0.001]

    bed_df = most_abundant_df[["ID","Refstart","Refend","ID","Score","Strand","Thikstart","Thikend","Colors"]]
    bed_df.iloc[:,0] = "RNA45SN1"
    template = 'track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"'

    with open(output + "most_abundant_" + output_file, 'w') as fp:
        fp.write(bed_df.to_csv(sep="\t",header=False,index=False))
    print(bed_df) 
        



# Select reads in size bins to perform a pearson correlation on a binary alignment signal. 
# If a read has a base at a given position it converts to 1 otherwise it converts to 0
# If binsize exceeds 100 reads the bins will be split and the pearson correlation will be computed seperately.
# Similarity will be min max normalized ((pearson_corr_matrix -1)/(-2)). Dissimilarites below certain threshold (0.2 is default) will be clustered.
# For clustering a list of pearson correlated dataframes will be stored and clustered in a follow up processing step


def pearson_correlation_fusion(pearson_alignment_df = pd.DataFrame(), iterator = 0 , batch_size = 100, overlap_threshold = overlap_threshold, reduce_size_bins_gt_1000 = 50, reduce_size_bins_gt_100 = 25, reduce_size_bins_gt_10 = 10): 
    if iterator == 2:
        print("Final print")
        print(pearson_alignment_df)
        return pearson_alignment_df
    elif iterator == 0:
        reference_dict = dict(pearson_alignment_df.iloc[0,:])
        template_length = int(reference_dict["Length"])
        list_of_lengths = size_bin_selection(template_length,reduce_size_bins_gt_1000,reduce_size_bins_gt_100,reduce_size_bins_gt_10)
        #print(list_of_lengths)
        list_read_groups = []
        list_read_singles = []
        single_reads_df = pd.DataFrame()

        for length_segment in tqdm(list_of_lengths, total = len(list_of_lengths)):
            minimum_length = length_segment[1]
            maximum_length = length_segment[0]
            temp_alignment_df = pearson_alignment_df.loc[pearson_alignment_df["Length"] <= maximum_length,:]
            temp_alignment_df = temp_alignment_df.loc[temp_alignment_df["Length"] >= minimum_length,:].sort_values(by=['Refstart','Length'], ascending=[True,False]).reset_index(drop=True)
            
            if temp_alignment_df.shape[0] > 0:
                if not temp_alignment_df.shape[0] > batch_size:
                    temp_alignment_json = temp_alignment_df.to_dict(orient="records")
                    with open(f"{output}temp_parallel_comp_{length_segment[0]}_{length_segment[1]}_small.json","w") as json_file:
                        json.dump(temp_alignment_json, json_file)
                    out_list_read_groups,out_list_read_singles, out_single_reads_df = pearson_correlation(f"{output}temp_parallel_comp_{length_segment[0]}_{length_segment[1]}_small.json",reference_dict)
                    for index,read_group in enumerate(out_list_read_groups):
                        list_read_groups.append(read_group)
                    for single in out_list_read_singles:
                        list_read_singles.append(single)
                    if out_single_reads_df.shape[0] != 0:
                        single_reads_df = pd.concat([single_reads_df,out_single_reads_df], axis = 0)
                else:
                    length = math.floor(temp_alignment_df.shape[0] / batch_size)
                    bins = [k for k in range(0,length + 1)]
                    list_for_parallel_comp = []
                    for index,k in enumerate(bins):
                        if k < length:
                            k_1 = k*batch_size
                            k_2 = (k+1)*batch_size
                            df_for_parallel_comp = temp_alignment_df.iloc[k*batch_size:((k+1)*batch_size),:]
                            df_for_parallel_comp_json = df_for_parallel_comp.to_dict(orient="records")
                            with open(f"{output}temp_parallel_comp_{length_segment[0]}_{length_segment[1]}_{k_1}_{k_2}_{index}_big.json","w") as json_file:
                                json.dump(df_for_parallel_comp_json, json_file)
                        else:
                            k_1 = k*batch_size
                            k_2 = len(temp_alignment_df)
                            df_for_parallel_comp = temp_alignment_df.iloc[k*batch_size:len(temp_alignment_df),:]
                            df_for_parallel_comp_json = df_for_parallel_comp.to_dict(orient="records")
                            with open(f"{output}temp_parallel_comp_{length_segment[0]}_{length_segment[1]}_{k_1}_{k_2}_{index}_big.json","w") as json_file:
                                json.dump(df_for_parallel_comp_json, json_file)
                        list_for_parallel_comp.append(f"{output}temp_parallel_comp_{length_segment[0]}_{length_segment[1]}_{k_1}_{k_2}_{index}_big.json")
                    with Pool(cores) as p:
                        pool_output = p.starmap(pearson_correlation, zip(list_for_parallel_comp,repeat(reference_dict)))
                        list_for_parallel_comp = []
                    for index,objects in enumerate(pool_output):
                        for index2,read_group in enumerate(objects[0]):
                            list_read_groups.append(read_group)
                        single_reads_df = pd.concat([single_reads_df,objects[2]], axis = 0)
                    p.close()

        #single_reads_df = single_reads_df.drop_duplicates(subset=["ID"])
        single_reads_df = single_reads_df.reset_index(drop=True)
        
        # Fusion dataframes become summarized in this step. All Sequences become summarized to a consensus sequence and the proportions of A,C,G,T and - content will be stored. 
        # The consensus sequence is determined by the majority of base abundancies for a given fusion group. The start and end of references is defined by the first and last occuring base in the consensus sequence.
        # A consensus dataframe is created with all the new consensus sequences
        consensus_rows = []         
        with Pool(cores) as p:
            pool_output = p.starmap(fusion_read_groups, zip([path for path in list_read_groups],repeat(reference_dict)))
        for objects in pool_output:
            out_dict = objects[0]
            if len(out_dict) != 0:
                consensus_rows.append(out_dict)
            single_reads_df = pd.concat([single_reads_df,objects[1]], axis = 0)
        p.close()
        for path in list_read_groups:
            os.remove(path)
        print("Single reads list\n")
        print(single_reads_df)
        print(sum(single_reads_df["n_Reads"]))
        print("\n")

        for idx, row in tqdm(single_reads_df.iterrows(), total=single_reads_df.shape[0]):
            if row["ID"] != "RNA45SN1":
                consensus_rows.append({"ID" : row["ID"], "Sequence" : row["Sequence"], "Proportion_Sequence": [], "Length": row["Length"], "Refstart": row["Refstart"], "Refend": row["Refend"],"n_Reads": row["n_Reads"],"IDS": [row["ID"]]})
        new_ref = {"ID" : str(reference), "Sequence" : str(fasta_file.fetch(reference)), "Proportion_Sequence": [], "Length": len(fasta_file.fetch(reference)), "Refstart": 0, "Refend": len(fasta_file.fetch(reference)),"n_Reads": 1,"IDS": []}
        consensus_rows.append(new_ref)
        consensus_df = pd.DataFrame.from_dict(consensus_rows)
        consensus_df = consensus_df.sort_values(by=["Refstart","Length"], ascending=[True,False]).reset_index(drop = True)

        #print("Consensus DF")
        print(consensus_df)
        print(sum(consensus_df["n_Reads"]))

        if pearson_alignment_df.shape[0] <= 5000 or batch_size > 300:
            final_iterator = 1
            if batch_size > 300:
                unclustered_df = consensus_df.loc[consensus_df["n_Reads"] == 1,:]
                unclustered_df.to_csv(output + "unclustered_reads_df.csv", sep = ";")
        else:
            final_iterator = 0

        new_batch_size = batch_size + 100
        summary_df = pearson_correlation_fusion(pearson_alignment_df = consensus_df, iterator = final_iterator, batch_size = new_batch_size, overlap_threshold=overlap_threshold, reduce_size_bins_gt_1000 = 100, reduce_size_bins_gt_100 = 50, reduce_size_bins_gt_10 = 25)
        return summary_df

    # This filtering step filters out consensus sequences that to not fit the reference to a certain threshold in terms of sequence matches (default 90%)
    elif iterator == 1:
        reference_dict = dict(pearson_alignment_df.iloc[0,:])
        list_read_groups = []
        list_read_singles = []
        single_reads_df = pd.DataFrame()
        pearson_alignment_df = pearson_alignment_df.loc[((pearson_alignment_df["ID"] != str(reference)))]
        unclustered_df = pearson_alignment_df.loc[pearson_alignment_df["n_Reads"] == 1,:]
        unclustered_df.to_csv(output + "unclustered_reads_df.csv", sep = ";")
        if batch_size > 300:
            pearson_alignment_df = pearson_alignment_df.loc[pearson_alignment_df["n_Reads"] > 1,:]
        pearson_alignment_json = pearson_alignment_df.to_dict(orient="records")
        with open(f"{output}pearson_alignment_df_final.json","w") as json_file:
            json.dump(pearson_alignment_json, json_file)

        print("Last lap")
        out_list_read_groups,out_list_read_singles, out_single_reads_df = pearson_correlation(f"{output}pearson_alignment_df_final.json",reference_dict)
        for index,read_group in enumerate(out_list_read_groups):
            list_read_groups.append(read_group)
        single_reads_df = pd.concat([single_reads_df,out_single_reads_df], axis = 0)
        #single_reads_df = single_reads_df.drop_duplicates(subset=["ID"])
    
        final_consensus_rows = []
        with Pool(cores) as p:
            pool_output = p.starmap(fusion_read_groups, zip([path for path in list_read_groups],repeat(reference_dict)))
        for objects in pool_output:
            out_dict = objects[0]
            if len(out_dict) != 0:
                final_consensus_rows.append(out_dict)
            single_reads_df = pd.concat([single_reads_df,objects[1]], axis = 0)
        p.close()
        for path in list_read_groups:
            os.remove(path)

        print(single_reads_df)
        #print(single_reads_df.shape[0])
        print(sum(single_reads_df["n_Reads"]))

        for idx, row in tqdm(single_reads_df.iterrows(), total=single_reads_df.shape[0]):
            if row["ID"] != "RNA45SN1":
                final_consensus_rows.append({"ID" : row["ID"], "Sequence" : row["Sequence"], "Proportion_Sequence": [], "Length": row["Length"], "Refstart": row["Refstart"], "Refend": row["Refend"],"n_Reads": row["n_Reads"],"IDS": [row["ID"]]})
        final_consensus_df = pd.DataFrame.from_dict(final_consensus_rows)
        final_consensus_df = final_consensus_df.sort_values(by=['Refstart','Length'], ascending=[True,False]).reset_index(drop = True)

        print(final_consensus_df)
        print(sum(final_consensus_df["n_Reads"]))
        # This filtering step filters out consensus sequences that to not fit the reference to a certain threshold in terms of sequence matches (default 90%)
        percentage_of_fitting_bases = []
        
        for i,row_sequence,rowstart,rowend in tqdm(zip(range(final_consensus_df.shape[0]),final_consensus_df["Sequence"],final_consensus_df["Refstart"],final_consensus_df["Refend"]),total = final_consensus_df.shape[0]):
            row_sequence = "-" * rowstart + row_sequence + "-" * (reference_dict["Length"] - rowend)
            if i == 0:
                percentage_of_fitting_bases.append(1)
                continue
            n_of_bases = 0
            n_fitting_bases_to_reference = 0
            for index,val in enumerate(reference_dict["Sequence"]):
                if row_sequence[index] in ["A","C","G","T","U"]:
                    n_of_bases += 1
                    if row_sequence[index] == val:
                        n_fitting_bases_to_reference += 1
            percentage_of_fitting_bases.append(n_fitting_bases_to_reference / n_of_bases)


        final_consensus_df["Percentage_of_fits"] = percentage_of_fitting_bases
        final_consensus_df = final_consensus_df.loc[((final_consensus_df["Percentage_of_fits"])) > 0.9,:]
        final_consensus_df = final_consensus_df.sort_values(by=['Refstart','Length'], ascending=[True,False]).reset_index(drop = True)
        new_iterator = 2
        unclustered_df_last = final_consensus_df.loc[final_consensus_df["n_Reads"] == 1,:]
        unclustered_df_last.to_csv(output + "unclustered_reads_df_after_last_lap.csv", sep = ";")
        final_consensus_df = final_consensus_df.loc[final_consensus_df["n_Reads"] != 1,:]
        summary_df = pearson_correlation_fusion(pearson_alignment_df = final_consensus_df, iterator = new_iterator, batch_size = batch_size, overlap_threshold=overlap_threshold, reduce_size_bins_gt_1000 = 100, reduce_size_bins_gt_100 = 50, reduce_size_bins_gt_10 = 25)
        return summary_df

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

logger.info("-----------------------------------------------------------------------------------------")
logger.info("---                            Similarity driven analysis                             ---")


logger.info("---                                     Load Samfile                                  ---")
samfile = pysam.AlignmentFile(bamfile_name, "rb")
number_of_reads = samfile.count()


logger.info("---                            Load Reference in fasta format                         ---")
fasta_file = pysam.FastaFile(fasta_name)
reference = fasta_file.references[0]



#Create a list of dictionairies with reference data as first entry, the samples will follow in further steps 
temp_list = list()
ref_data = {"ID" : str(reference), "Sequence" : str(fasta_file.fetch(reference)),"Proportion_Sequence": [], "Length": len(fasta_file.fetch(reference)), "Refstart": 0, "Refend": len(fasta_file.fetch(reference)), "n_Reads": 1, "IDS": []}
temp_list.append(ref_data)


# Create a predictionairy for each read the Bam file contains 
# Read Sequence is primed with dashes in legnth of reference string
# Length of the sequence is primed by reference_end - reference_start
# Create alignment dataframe will primed placeholders 
logger.info("-----------------------------------------------------------------------------------------")
logger.info("---                Create alignment with sequence and cigarstring of bam file         ---")
cigartuples = []
read_ids = []
read_sequences = []
read_starts = []
read_ends = []
counter_forward = 0
counter_reverse = 0
for read in tqdm(samfile.fetch(),total=number_of_reads):
    read_data = {"ID" : read.query_name,"Sequence" : "-" * int(read.reference_end - read.reference_start),"Proportion_Sequence": [], "Length": int(read.reference_end - read.reference_start), "Refstart": read.reference_start, "Refend": read.reference_end, "n_Reads": 1,"IDS": [read.query_name]}
    temp_list.append(read_data)
    if read.is_forward:
        counter_forward += 1
        read_sequences.append(read.get_forward_sequence())
    else:
        counter_reverse += 1
        read_sequences.append(read.query_sequence)
    cigartuples.append(read.cigartuples)
    read_ids.append(read.query_name)
    
    read_starts.append(read.reference_start)
    read_ends.append(read.reference_end)
alignment_df = pd.DataFrame.from_dict(temp_list)
query_names = list()
# Create an alignment representation of each read to the RNA45S reference to locate Start, End and Length of the read
# Alignment is performed by following the 
# Sort read by start and readlength such that the read with the ealiest start point of longest length will occur first in the table
final_strings = []
for c_tuple_read, read_id, read_sequence, start_point_ref, end_point_ref in tqdm(zip(cigartuples,read_ids, read_sequences, read_starts, read_ends), total = len(read_ids)):
    string = ""
    position_read = 0
    for c_tuple in c_tuple_read: 
        operation = c_tuple[0]
        length = c_tuple[1]
        difference = 0
        if operation == 0:
            string += read_sequence[position_read:(position_read + length)]
            position_read += length
        elif operation == 1:
            position_read += length
        elif operation == 2:
            string += "-"*length
        elif operation == 3:
            string += read_sequence[position_read] * length
        elif operation == 4:
            position_read += length
    final_strings.append(string)

read_ids.insert(0,"RNA45SN1")
final_strings.insert(0,alignment_df.iloc[0]["Sequence"])
alignment_df["Sequence"] = final_strings 



indices_list = [index for index in range(alignment_df.shape[0]-1)]
intermediate_list = list(zip(alignment_df["Sequence"],alignment_df["Refstart"],alignment_df["Refend"]))
read_list_df = [intermediate_list[i] for i in tqdm(indices_list)]
reference_len = len(fasta_file.fetch(reference))
reference_seq = str(fasta_file.fetch(reference))



logger.info("---                Quality check for reference and query identity                     ---")
with Pool(cores) as p:
    pool_output = p.starmap(percentage_of_fits_parrallel, zip(indices_list,read_list_df,repeat(reference_len),repeat(reference_seq)))
selected_indices = []
for objects in pool_output:
    if objects[0] >= 0:
        selected_indices.append(objects[0])
read_list_df = []
p.close()

alignment_df = alignment_df.iloc[selected_indices,:]
alignment_df = alignment_df.sort_values(by=['Refstart','Length'], ascending=[True,False])
alignment_df.index = alignment_df["ID"]


logger.info("---                      Perform pearson correlation and clustering                     ---")
final_df = pearson_correlation_fusion(alignment_df,0,100)
final_df["rel_n_Reads"] = final_df["n_Reads"] / number_of_reads
final_df.to_csv(output + "fragment_df.csv", sep=";")



logger.info("---                                  Write results                                     ---")
simple_consensus_df = pd.DataFrame(final_df[["ID","Refstart","Refend","Length","n_Reads","rel_n_Reads","Percentage_of_fits"]])
simple_consensus_df.to_csv(output + "fragment_df_simple.csv", sep=";")
bed_consensus_df = final_df[["ID","Refstart","Refend","ID"]]
bed_consensus_df.iloc[:,0] = "RNA45SN1"
create_colored_bed(table_name = output + "fragment_df.csv", output_folder = output, output_file = "no_template.bed", sample_type = sample_type)
alignment_df.to_csv(output + "alignment_df.csv", sep = ";")

logger.info("---                                      Done                                          ---")
logger.info("------------------------------------------------------------------------------------------")


logger.info(f"### Metainformation ###")
logger.info(f"### Number of reads: {number_of_reads} ###")
logger.info(f"### Number of forward reads: {counter_forward} ###")
logger.info(f"### Number of reverse reads: {counter_reverse} ###\n")

logger.info("### Results after similarity clustering ###")
logger.info(simple_consensus_df)
logger.info("############################################\n")
