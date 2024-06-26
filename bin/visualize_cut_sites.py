import pysam
from tqdm import tqdm
from pathlib import Path
import operator
import os
import pandas as pd
import json
import argparse
import plotly.graph_objects as go
import plotly
from random import seed,uniform
import plotly.io as pio

opt_parser = argparse.ArgumentParser()

opt_parser.add_argument("-s", "--start_sites_df", dest="start_sites_df", help="Insert a sample bam file", metavar="FILE")
opt_parser.add_argument("-e", "--end_sites_df",dest="end_sites_df", help="Insert a reference fasta file", metavar="FILE")
opt_parser.add_argument("-r", "--reference",dest="reference", help="Insert a path to the reference file", metavar="FILE")
opt_parser.add_argument("-o", "--output_path",dest="output_path", help="Insert an output path for the plot", metavar="FILE")

options = opt_parser.parse_args()

start_sites_df_path = options.start_sites_df
end_sites_df_path = options.end_sites_df
reference_path = options.reference
output_path = options.output_path


fasta_file = pysam.FastaFile(reference_path)
reference = fasta_file.references[0]
reference_sequence = str(fasta_file.fetch(reference))


start_sites_df = pd.read_csv(start_sites_df_path,sep="\t")
start_sites_df.columns = ["reference","start","end","type","score","strand","thikstart","thikend","color"]
start_sites_df["score"] = start_sites_df["score"]/1000000
end_sites_df = pd.read_csv(end_sites_df_path,sep="\t")
end_sites_df.columns = ["reference","start","end","type","score","strand","thikstart","thikend","color"]
end_sites_df["score"] = end_sites_df["score"]/1000000


fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x = [i for i in range(0,len(reference_sequence))],
        y = [0 for i in range(0,len(reference_sequence))],
        name= "Cut type",
        line_color = "rgba(256,256,256,0.01)",
        )
    )


fig.add_trace(
    go.Scatter(
        x = start_sites_df["start"],
        y = start_sites_df["score"],
        name="start site",
        line_color = "green",
        mode = "markers"
        )
    )

fig.add_trace(
    go.Scatter(
        x = end_sites_df["start"],
        y = end_sites_df["score"],
        name= "end site",
        line_color = "darkred",
        mode = "markers"
        )
    )


fig.update_layout(
    title=f"Cut-sites",
    xaxis=dict(
        title="Position on reference",
        gridcolor = "lightgray",
        tickvals = [i for i in range(0,len(reference_sequence),100)],
        ticktext = [str(i) for i in range(0,len(reference_sequence),100)]
        ),
    yaxis=dict(
        title="Cutsite frequency",
        gridcolor = "lightgray",
        range=[0,1]
        ),
    plot_bgcolor='rgba(0,0,0,0)'
)

plotly.offline.plot(fig, filename=f"{output_path}/cut_sites.html")





