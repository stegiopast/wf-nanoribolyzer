from pathlib import Path
import numpy as np
import pysam
from tqdm import tqdm
from itertools import repeat
import polars as pl
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import plotly
import argparse
import json
from random import seed,uniform

# Takes bam file products of dorado when running it
# with modification detection flags.
# Extracts modification ratios. 
# Visualization of detected modifications as lineplot. 
# Works only on directRNA samples.

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

opt_parser.add_argument(
    "-b",
    "--bamfile",
    dest="bamfile_name",
    help="Insert a bamfile holding modification tags (basecalled with modification detection flags)",
    metavar="FILE",
)

opt_parser.add_argument(
    "-r",
    "--reference",
    dest="reference",
    help="Insert a path to the reference file",
    metavar="FILE",
)

opt_parser.add_argument(
    "-m",
    "--modification_list",
    dest="modification_list",
    help="Insert the list of literature based modifications",
    metavar="FILE",
)
opt_parser.add_argument(
    "-o",
    "--output_path",
    dest="output_path",
    help="Insert an output filepath to write to",
    metavar="FILE",
)



options = opt_parser.parse_args()

bamfile_path = options.bamfile_name#"/home/stefan/analysis_modifications_dRNA/IVPA_Nuc_filtered_pod5_rebasecalled_psU_m6A_aligned.bam"
reference_path = options.reference
literature_mod_df_path = options.modification_list
output_path = options.output_path


#####################################################################################################################################
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                            Script                                                                 #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#                                                                                                                                   #
#####################################################################################################################################

literature_mod_df = pd.read_csv(literature_mod_df_path ,sep="\t",header=None,index_col=None) #"/home/stefan/wf-nanoribolyzer/references/rRNA_modifications_conv.bed"
possible_columns = ["reference","start","end","modification","A","B","C","D","E","F","G","H"]
literature_mod_df.columns = possible_columns[0:len(literature_mod_df.columns)]



fasta_file = pysam.FastaFile(reference_path)
reference = fasta_file.references[0]
reference_sequence = str(fasta_file.fetch(reference))
reads_aligning = [0 for i in range(len(reference_sequence))]

###Adenosin###
print(literature_mod_df)
m6a_mod_df = literature_mod_df.loc[literature_mod_df["modification"] == "m6A"]
Am_mod_df = literature_mod_df.loc[literature_mod_df["modification"] == "Am"]
Ino_mod_df = literature_mod_df.loc[literature_mod_df["modification"] == "Ino"]
m62a_mod_df = literature_mod_df.loc[literature_mod_df["modification"] == "m62A"]

mod_positions_m6a = [0 for i in range(len(reference_sequence))]
mod_ids_m6a = [[] for i in range(len(reference_sequence))]
mod_positions_Am = [0 for i in range(len(reference_sequence))]
mod_ids_Am = [[] for i in range(len(reference_sequence))]
mod_positions_Ino = [0 for i in range(len(reference_sequence))]
mod_ids_Ino = [[] for i in range(len(reference_sequence))]

##Uridine###
psu_mod_df = literature_mod_df.loc[literature_mod_df["modification"] == "psu"]
Um_mod_df = literature_mod_df.loc[literature_mod_df["modification"] == "Um"]

mod_positions_pseU = [0 for i in range(len(reference_sequence))]
mod_ids_pseU = [[] for i in range(len(reference_sequence))]
mod_positions_Um = [0 for i in range(len(reference_sequence))]
mod_ids_Um = [[] for i in range(len(reference_sequence))]
number_of_basecalled_C = [0 for i in range(len(reference_sequence))]

###Guanosin###
Gm_mod_df = literature_mod_df.loc[literature_mod_df["modification"] == "Gm"]

mod_positions_Gm = [0 for i in range(len(reference_sequence))]
mod_ids_Gm = [[] for i in range(len(reference_sequence))]



###Cytosin###
Cm_mod_df = literature_mod_df.loc[literature_mod_df["modification"] == "Cm"]

mod_positions_Cm = [0 for i in range(len(reference_sequence))]
mod_ids_Cm = [[] for i in range(len(reference_sequence))]





bamfile = pysam.AlignmentFile(bamfile_path, mode="rb")
for i in tqdm(bamfile.fetch(until_eof=True)):
    if i.is_supplementary:
        continue
    start = i.reference_start
    end = i.reference_end
    for index in range(start,end):
        reads_aligning[index] += 1
    mod_obj = i.modified_bases
    if mod_obj != None:
        try:
            mod_m6a = list(mod_obj[('A', 0, 'a')])
        except KeyError:
            mod_m6a = None
        try:
            mod_Am = list(mod_obj[('A', 0, 69426)])
        except:
            mod_Am = None
        try:
            mod_Ino = list(mod_obj[('A', 0, 17596)])
        except:
            mod_Ino = None
        try:
            mod_pseU = list(mod_obj[('T', 0, 17802)])
        except KeyError:
            mod_pseU = None
        try:
            mod_Um = list(mod_obj[('T', 0, 19227)])
        except KeyError:
            mod_Um = None
        try:
            mod_Gm = list(mod_obj[('G', 0, 19229)])
        except KeyError:
            mod_Gm = None
        try:
            mod_Cm = list(mod_obj[('C', 0, 19228)])
        except KeyError:
            mod_Cm = None
            
        aligned_pairs = i.get_aligned_pairs(with_seq=True)
        alignment_dict = {}
        for pair_element in aligned_pairs:
            if None not in pair_element:
                alignment_dict[str(pair_element[0])] = {"index_query":pair_element[0],"index_reference":pair_element[1],"base_query": str(i.get_forward_sequence())[pair_element[0]],"base_reference": reference_sequence[pair_element[1]]}
        
        
        ###Adenin###
        if mod_m6a != None:
            for mod_base in mod_m6a:
                p = ((mod_base[1] + 1)/256)
                if p >= 0.95 and str(mod_base[0]) in alignment_dict:
                    mod_positions_m6a[alignment_dict[str(mod_base[0])]["index_reference"]] += 1
                    mod_ids_m6a[alignment_dict[str(mod_base[0])]["index_reference"]].append(i.query_name)
                    
        if mod_Am != None:
            for mod_base in mod_Am:
                p = ((mod_base[1] + 1)/256)
                if p >= 0.95 and str(mod_base[0]) in alignment_dict:
                    mod_positions_Am[alignment_dict[str(mod_base[0])]["index_reference"]] += 1
                    mod_ids_Am[alignment_dict[str(mod_base[0])]["index_reference"]].append(i.query_name)
        
        if mod_Ino != None:
            for mod_base in mod_Ino:
                p = ((mod_base[1] + 1)/256)
                if p >= 0.95 and str(mod_base[0]) in alignment_dict:
                    mod_positions_Ino[alignment_dict[str(mod_base[0])]["index_reference"]] += 1
                    mod_ids_Ino[alignment_dict[str(mod_base[0])]["index_reference"]].append(i.query_name)

        
        ###Uridine###
        if mod_pseU != None:
            for mod_base in mod_pseU:
                p = ((mod_base[1] + 1)/256)
                if p >= 0.95 and str(mod_base[0]) in alignment_dict:
                    mod_positions_pseU[alignment_dict[str(mod_base[0])]["index_reference"]] += 1
                    mod_ids_pseU[alignment_dict[str(mod_base[0])]["index_reference"]].append(i.query_name)
                    
        for pair_element in aligned_pairs:
            if None not in pair_element:
                if alignment_dict[str(pair_element[0])]["base_query"] == "C" and alignment_dict[str(pair_element[0])]["base_reference"] == "T":
                    number_of_basecalled_C[pair_element[1]] += 1
                    
        if mod_Um != None:
            for mod_base in mod_Um:
                p = ((mod_base[1] + 1)/256)
                if p >= 0.95 and str(mod_base[0]) in alignment_dict:
                    mod_positions_Um[alignment_dict[str(mod_base[0])]["index_reference"]] += 1
                    mod_ids_Um[alignment_dict[str(mod_base[0])]["index_reference"]].append(i.query_name)
                    
        ###Guanin###
        if mod_Gm != None:
            for mod_base in mod_Gm:
                p = ((mod_base[1] + 1)/256)
                if p >= 0.95 and str(mod_base[0]) in alignment_dict:
                    mod_positions_Gm[alignment_dict[str(mod_base[0])]["index_reference"]] += 1
                    mod_ids_Gm[alignment_dict[str(mod_base[0])]["index_reference"]].append(i.query_name)
        
        ###Cytosin###
        if mod_Cm != None:
            for mod_base in mod_Cm:
                p = ((mod_base[1] + 1)/256)
                if p >= 0.95 and str(mod_base[0]) in alignment_dict:
                    mod_positions_Cm[alignment_dict[str(mod_base[0])]["index_reference"]] += 1
                    mod_ids_Cm[alignment_dict[str(mod_base[0])]["index_reference"]].append(i.query_name)


reads_aligning = np.array(reads_aligning)
reads_aligning[reads_aligning == 0] = 1


positions = [i+1 for i in range(len(reads_aligning))]

mod_positions_m6a = np.array(mod_positions_m6a)
rel_mod_positions_m6a = mod_positions_m6a/reads_aligning

mod_positions_Am = np.array(mod_positions_Am)
rel_mod_positions_Am = mod_positions_Am/reads_aligning

mod_positions_Ino = np.array(mod_positions_Ino)
rel_mod_positions_Ino = mod_positions_Ino/reads_aligning

mod_positions_pseU = np.array(mod_positions_pseU)
rel_mod_positions_pseU = mod_positions_pseU/reads_aligning

mod_positions_Um = np.array(mod_positions_Um)
rel_mod_positions_Um = mod_positions_Um/reads_aligning

mod_positions_Gm = np.array(mod_positions_Gm)
rel_mod_positions_Gm = mod_positions_Gm/reads_aligning

mod_positions_Cm = np.array(mod_positions_Cm)
rel_mod_positions_Cm = mod_positions_Cm/reads_aligning

number_of_basecalled_C = np.array(number_of_basecalled_C)
rel_number_of_basecalled_C = number_of_basecalled_C/reads_aligning

seed(90)


dataset = {
    "position":positions,
    "reads_aligning":reads_aligning,
    "n_m6a":mod_positions_m6a,
    "n_Am":mod_positions_Am,
    "n_Ino":mod_positions_Ino,
    "n_pseU":mod_positions_pseU,
    "n_C":number_of_basecalled_C,
    "n_Um": mod_positions_Um,
    "n_Gm": mod_positions_Gm,
    "n_Cm": mod_positions_Cm,
    "rel_n_pseU":rel_mod_positions_pseU,
    "rel_n_C":number_of_basecalled_C / reads_aligning,
    "rel_n_C_and_PseU": rel_number_of_basecalled_C + rel_mod_positions_pseU,
    "rel_n_m6a":rel_mod_positions_m6a,
    "rel_n_Am":rel_mod_positions_Am,
    "rel_n_Ino": rel_mod_positions_Ino,
    "rel_n_Um":rel_mod_positions_Um,
    "rel_n_Gm":rel_mod_positions_Gm,
    "rel_n_Cm":rel_mod_positions_Cm,
}

modification_df = pd.DataFrame(dataset)
modification_df.to_csv(f"{output_path}/modification_quantification.csv",sep=";",header=True,index=None)

modified_ids_dataset = {
    "m6a": mod_ids_m6a,
    "Am": mod_ids_Am,
    "Ino": mod_ids_Ino,
    "pseU": mod_ids_pseU,
    "Um": mod_ids_Um,
    "Gm": mod_ids_Gm,
    "Cm": mod_ids_Cm    
}

json_dict = json.dumps(modified_ids_dataset)

"""
How to read json 
with open("modification_ids.json") as json_file:
   json_entry = json.load(json_file)
dictionairy = dict(eval(json_entry))
"""


with open(f"{output_path}/modification_ids.json", "w") as dict_file:
    json.dump(json_dict, dict_file)


layout = go.Layout(height = 800)
fig = go.Figure(layout=layout)


fig.add_trace(
    go.Scatter(
        x=[index for index in range(len(rel_mod_positions_pseU))],
        y=[value for index,value in enumerate(rel_mod_positions_pseU)],  # Use the same x position and y at the top of the line
        line_color="rgba(0,0,128,1)",
        showlegend=True,
        name="pseU freq."
    )
)


fig.add_trace(
        go.Scatter(
            x=[index for index in range(len(rel_number_of_basecalled_C))],
            y=[rel_mod_positions_pseU[index]+value for index,value in enumerate(rel_number_of_basecalled_C)],  # Use the same x position and y at the top of the line
            line_color="rgba(0,0,128,0.5)",
            showlegend=True,
            name="C/U missmatches freq."
        )
    )


fig.add_trace(
    go.Scatter(
        x=[index for index in range(len(rel_mod_positions_Um))],
        y=[value for index,value in enumerate(rel_mod_positions_Um)],  # Use the same x position and y at the top of the line
        line_color="rgba(128,0,0,1)",
        showlegend=True,
        name="Um freq."
    )
)




fig.add_trace(
    go.Scatter(
        x = [i for i in range(0,len(reference_sequence))],
        y = [0 for i in range(0,len(reference_sequence))],
        name= "",
        line_color = "white",
        )
    )

for index,T_position in psu_mod_df.iterrows():
        fig.add_shape(
            x0=T_position["end"]-1,
            x1=T_position["end"]-1,
            y0=0,
            y1=1,
            line=dict(
                color="rgba(0,0,128,0.5)",
                width=0.35,
                dash="dash"
                )
            )
    
for index,T_position in Um_mod_df.iterrows():
    fig.add_shape(
        x0=T_position["end"]-1,
        x1=T_position["end"]-1,
        y0=0,
        y1=1,
        line=dict(
            color="rgba(128,0,0,0.5)",
            width=0.35,
            dash="dash"
            )
        )

fig.add_trace(go.Scatter(
    x=[None], y=[None],  # Invisible point, used only for legend entry
    mode='lines',
    line=dict(
        color="rgba(0,0,128,0.5)",
        width=1,
        dash="dash"
    ),
    showlegend=True,
    name="known pseU"  # Legend entry name
))

fig.add_trace(go.Scatter(
    x=[None], y=[None],  # Invisible point, used only for legend entry
    mode='lines',
    line=dict(
        color="rgba(128,0,0,0.5)",
        width=1,
        dash="dash"
    ),
    showlegend=True,
    name="known Um"  # Legend entry name
))

    

fig.update_layout(
    title="Uridine based modifications",
    xaxis=dict(title="Position on reference",gridcolor = "white",tickformat="d"),
    yaxis=dict(title="Modification frequency",gridcolor = "white"),
    plot_bgcolor='rgba(0,0,0,0)',
)

plotly.io.write_html(fig, f"{output_path}/uridine_modification_abundance.html")








layout = go.Layout(height = 800)
fig2 = go.Figure(layout=layout)
fig2.add_trace(
    go.Scatter(
        x=[index for index in range(len(rel_mod_positions_m6a))],
        y=[value for index,value in enumerate(rel_mod_positions_m6a)],  # Use the same x position and y at the top of the line
        line_color="rgba(0, 0, 128, 1)",
        showlegend=True,
        name="m6a freq."
    )
)


fig2.add_trace(
    go.Scatter(
        x=[index for index in range(len(rel_mod_positions_Am))],
        y=[value for index,value in enumerate(rel_mod_positions_Am)],  # Use the same x position and y at the top of the line
        line_color="rgba(0, 128, 0, 1)",
        showlegend=True,
        name="Am freq."
    )
)


fig2.add_trace(
    go.Scatter(
        x=[index for index in range(len(rel_mod_positions_Ino))],
        y=[value for index,value in enumerate(rel_mod_positions_Ino)],  # Use the same x position and y at the top of the line
        line_color="rgba(128, 0, 0, 1)",
        showlegend=True,
        name="Ino freq."
    )
)

fig2.add_trace(
    go.Scatter(
        x = [i for i in range(0,len(reference_sequence))],
        y = [0 for i in range(0,len(reference_sequence))],
        name= "",
        line_color = "white",
        )
    )


for index,A_position in m6a_mod_df.iterrows():
    fig2.add_shape(
        x0=A_position["end"],
        x1=A_position["end"],
        y0=0,
        y1=1,
        line=dict(
            color="rgba(0,0,128,0.5)",
            width=0.35,
            dash="dash"
            )
        ) 
    
print(Am_mod_df)
for index,A_position in Am_mod_df.iterrows():
    fig2.add_shape(
        x0=A_position["end"]-1,
        x1=A_position["end"]-1,
        y0=0,
        y1=1,
        line=dict(
            color="rgba(0,128,0,0.5)",
            width=0.35,
            dash="dash"
            )
        )


fig2.add_trace(go.Scatter(
    x=[None], y=[None],  # Invisible point, used only for legend entry
    mode='lines',
    line=dict(
        color="rgba(0,0,128,0.5)",
        width=1,
        dash="dash"
    ),
    showlegend=True,
    name="known m6a"  # Legend entry name
))


fig2.add_trace(go.Scatter(
    x=[None], y=[None],  # Invisible point, used only for legend entry
    mode='lines',
    line=dict(
        color="rgba(0,128,0,0.5)",
        width=1,
        dash="dash"
    ),
    showlegend=True,
    name="known Am"  # Legend entry name
))

fig2.update_layout(
    title="Adenine based modifications",
    xaxis=dict(title="Position on reference",gridcolor = "white",tickformat="d"),
    yaxis=dict(title="Modification frequency",gridcolor = "white"),
    plot_bgcolor='rgba(0,0,0,0)'
)

plotly.io.write_html(fig2, f"{output_path}/adenine_modification_abundance.html")



layout = go.Layout(height = 800)
fig3 = go.Figure(layout=layout)
fig3.add_trace(
    go.Scatter(
        x=[index for index in range(len(rel_mod_positions_Gm))],
        y=[value for index,value in enumerate(rel_mod_positions_Gm)],  # Use the same x position and y at the top of the line
        line_color="rgba(0, 0, 128, 1)",
        showlegend=True,
        name="Gm freq."
    )
)


fig3.add_trace(
    go.Scatter(
        x = [i for i in range(0,len(reference_sequence))],
        y = [0 for i in range(0,len(reference_sequence))],
        name= "",
        line_color = "white",
        )
    )


for index,A_position in Gm_mod_df.iterrows():
    fig3.add_shape(
        x0=A_position["end"],
        x1=A_position["end"],
        y0=0,
        y1=1,
        line=dict(
            color="rgba(0,0,128,0.5)",
            width=0.2,
            dash="dash"
            )
        ) 
fig3.add_trace(go.Scatter(
    x=[None], y=[None],  # Invisible point, used only for legend entry
    mode='lines',
    line=dict(
        color="rgba(0,0,128,0.5)",
        width=1,
        dash="dash"
    ),
    showlegend=True,
    name="known Gm"  # Legend entry name
))
fig3.update_layout(
    title="Guanine based modifications",
    xaxis=dict(title="Position on reference",gridcolor = "white",tickformat="d"),
    yaxis=dict(title="Modification frequency",gridcolor = "white"),
    plot_bgcolor='rgba(0,0,0,0)'
)
plotly.io.write_html(fig3, f"{output_path}/guanine_modification_abundance.html")




layout = go.Layout(height = 800)
fig4 = go.Figure(layout=layout)
fig4.add_trace(
    go.Scatter(
        x=[index for index in range(len(rel_mod_positions_Cm))],
        y=[value for index,value in enumerate(rel_mod_positions_Cm)],  # Use the same x position and y at the top of the line
        line_color="rgba(0, 128, 0, 1)",
        showlegend=True,
        name="Cm freq."
    )
)
fig4.add_trace(
    go.Scatter(
        x = [i for i in range(0,len(reference_sequence))],
        y = [0 for i in range(0,len(reference_sequence))],
        name= "",
        line_color = "white",
        )
    )
for index,A_position in Cm_mod_df.iterrows():
    fig4.add_shape(
        x0=A_position["end"],
        x1=A_position["end"],
        y0=0,
        y1=1,
        line=dict(
            color="rgba(0,128,0,0.5)",
            width=0.2,
            dash="dash"
            )
        ) 
fig4.add_trace(go.Scatter(
    x=[None], y=[None],  # Invisible point, used only for legend entry
    mode='lines',
    line=dict(
        color="rgba(0,128,0,0.5)",
        width=1,
        dash="dash"
    ),
    showlegend=True,
    name="known Cm"  # Legend entry name
))
fig4.update_layout(
    title="Cytosine based modifications",
    xaxis=dict(title="Position on reference",gridcolor = "white",tickformat="d"),
    yaxis=dict(title="Modification frequency",gridcolor = "white"),
    plot_bgcolor='rgba(0,0,0,0)'
)
plotly.io.write_html(fig4, f"{output_path}/cytosine_modification_abundance.html")
