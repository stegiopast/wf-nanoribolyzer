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
from random import seed,uniform



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

literature_mod_df = pd.read_csv(literature_mod_df_path ,sep="\t",header=None,index_col=None) #"/home/stefan/wf-nanoribolyzer/references/rRNA_modifications_conv.bed"
literature_mod_df.columns = ["reference","start","end","modification","A","B","C"]

T_mod_df = literature_mod_df.loc[literature_mod_df["modification"] == "psu"]
A_mod_df = literature_mod_df.loc[literature_mod_df["modification"].isin(["Am","m62A","m6A"])]


fasta_file = pysam.FastaFile(reference_path)
reference = fasta_file.references[0]
reference_sequence = str(fasta_file.fetch(reference))

reads_aligning = [0 for i in range(len(reference_sequence))]
mod_positions_m6a = [0 for i in range(len(reference_sequence))]
mod_positions_pseU = [0 for i in range(len(reference_sequence))]
number_of_basecalled_C = [0 for i in range(len(reference_sequence))]

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
            mod_pseU = list(mod_obj[('T', 0, 17802)])
        except KeyError:
            mod_pseU = None
        aligned_pairs = i.get_aligned_pairs(with_seq=True)
        alignment_dict = {}
        for pair_element in aligned_pairs:
            if None not in pair_element:
                alignment_dict[str(pair_element[0])] = {"index_query":pair_element[0],"index_reference":pair_element[1],"base_query": str(i.get_forward_sequence())[pair_element[0]],"base_reference": reference_sequence[pair_element[1]]}
        if mod_m6a != None:
            for mod_base in mod_m6a:
                p = ((mod_base[1] + 1)/256)
                if p >= 0.95 and str(mod_base[0]) in alignment_dict:
                    mod_positions_m6a[alignment_dict[str(mod_base[0])]["index_reference"]] += 1

        if mod_pseU != None:
            for mod_base2 in mod_pseU:
                p = ((mod_base2[1] + 1)/256)
                if p >= 0.95 and str(mod_base2[0]) in alignment_dict:
                    mod_positions_pseU[alignment_dict[str(mod_base2[0])]["index_reference"]] += 1
                    
        for pair_element in aligned_pairs:
            if None not in pair_element:
                if alignment_dict[str(pair_element[0])]["base_query"] == "C" and alignment_dict[str(pair_element[0])]["base_reference"] == "T":
                    number_of_basecalled_C[pair_element[1]] += 1



reads_aligning = np.array(reads_aligning)
reads_aligning[reads_aligning == 0] = 1


positions = [i+1 for i in range(len(reads_aligning))]

mod_positions_m6a = np.array(mod_positions_m6a)
rel_mod_positions_m6a = mod_positions_m6a/reads_aligning

mod_positions_pseU = np.array(mod_positions_pseU)
rel_mod_positions_pseU = mod_positions_pseU/reads_aligning

number_of_basecalled_C = np.array(number_of_basecalled_C)
rel_number_of_basecalled_C = number_of_basecalled_C/reads_aligning

seed(90)


dataset = {
    "position":positions,
    "reads_aligning":reads_aligning,
    "n_pseU":mod_positions_pseU,
    "n_m6a":number_of_basecalled_C,
    "n_C":number_of_basecalled_C,
    "rel_n_pseU":rel_mod_positions_pseU,
    "rel_n_m6a":rel_mod_positions_m6a,
    "rel_n_C":rel_number_of_basecalled_C
}

modification_df = pd.DataFrame(dataset)
modification_df.to_csv(f"{output_path}/modification_quantification.csv",sep=";",header=True,index=None)

fig = go.Figure()


fig = go.Figure()


fig.add_trace(
    go.Scatter(
        x=[index for index in range(len(rel_mod_positions_pseU))],
        y=[value for index,value in enumerate(rel_mod_positions_pseU)],  # Use the same x position and y at the top of the line
        line_color="rgba(153,0,0,1)",
        showlegend=True,
        name="pseU freq."
    )
)



fig.add_trace(
        go.Scatter(
            x=[index for index in range(len(rel_number_of_basecalled_C))],
            y=[rel_mod_positions_pseU[index]+value for index,value in enumerate(rel_number_of_basecalled_C)],  # Use the same x position and y at the top of the line
            line_color="rgba(153,0,0,0.4)",
            showlegend=True,
            name="C/U missmatches freq."
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

for index,T_position in T_mod_df.iterrows():
    fig.add_shape(
        x0=T_position["end"],
        x1=T_position["end"],
        y0=0,
        y1=1,
        line=dict(
            color="rgba(255,106,106,0.5)",
            width=1,
            dash="dash"
            )
        )
    modification = T_position["modification"]
    fig.add_annotation(
        x=T_position["end"],
        y = uniform(0.8,0.99),
        text="pseU",
        arrowhead=1,
        showarrow=False,
        ax=0,  
        ay=-40,
        font=dict(size=10) 
        )
    

fig.update_layout(
    title=f"Modification basecalling",
    xaxis=dict(title="Position on reference",gridcolor = "white"),
    yaxis=dict(title="Modification frequency",gridcolor = "white"),
    plot_bgcolor='rgba(0,0,0,0)',
)

plotly.offline.plot(fig, filename=f"{output_path}/relative_pseU_modification_abundance.html")


fig2 = go.Figure()
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
        x = [i for i in range(0,len(reference_sequence))],
        y = [0 for i in range(0,len(reference_sequence))],
        name= "",
        line_color = "white",
        )
    )


for index,A_position in A_mod_df.iterrows():
    fig2.add_shape(
        x0=A_position["end"],
        x1=A_position["end"],
        y0=0,
        y1=1,
        line=dict(
            color="rgba(176,196,222,0.7)",
            width=2,
            dash="dash"
            )
        ) 
    modification = A_position["modification"]
    fig2.add_annotation(
        x=A_position["end"],
        y = uniform(0.8,0.99),
        text=f"{modification}",
        showarrow=False,
        arrowhead=1,
        ax=0,  
        ay=-40,
        font=dict(size=10)
        )


fig2.update_layout(
    title=f"Modification basecalling",
    xaxis=dict(title="Position on reference",gridcolor = "white"),
    yaxis=dict(title="Modification frequency",gridcolor = "white"),
    plot_bgcolor='rgba(0,0,0,0)'
)

plotly.offline.plot(fig2, filename=f"{output_path}/relative_m6A_modification_abundance.html")