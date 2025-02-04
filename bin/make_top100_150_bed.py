import pandas as pd
import argparse


# Argument Parser
opt_parser = argparse.ArgumentParser()

opt_parser.add_argument(
    "-i",
    "--input_file",
    dest="bedfile_name",
    help="Insert a sample bam file",
    metavar="FILE",
)

opt_parser.add_argument(
    "-o",
    "--output_path",
    dest="output",
    help="Insert an output directory to write to",
    metavar="FILE",
)

options = opt_parser.parse_args()
bedfile_name = options.bedfile_name
output = options.output


def top150_colored_bed(
    table_name: str = "",
    output_folder: str = output
):
    df = pd.read_csv(table_name, sep="\t", header=0)
    df.columns = [
            "ID",
            "Refstart",
            "Refend",
            "ID",
            "Score",
            "Strand",
            "Thikstart",
            "Thikend",
            "Colors"
            ]
    
    top100_df = df.sort_values(by="Score",ascending=False)[0:min(100,df.shape[0])]
    bed_df = top100_df
    bed_df.iloc[:, 0] = "RNA45SN1"
    with open(f"{output}top_100_no_template.bed", "w") as fp:
        fp.write(bed_df.to_csv(sep="\t", header=False, index=False))
        
    top150_df = df.sort_values(by="Score",ascending=False)[0:min(150,df.shape[0])]
    bed_df = top150_df
    bed_df.iloc[:, 0] = "RNA45SN1"
    with open(f"{output}top_150_no_template.bed", "w") as fp:
        fp.write(bed_df.to_csv(sep="\t", header=False, index=False))
        
top150_colored_bed(bedfile_name,output)