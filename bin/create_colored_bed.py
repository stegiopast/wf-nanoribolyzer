import argparse
import pandas as pd




opt_parser = argparse.ArgumentParser()

opt_parser.add_argument("-i", "--input_file", dest="csvfile_name", help="Insert a sample fragment csv file", metavar="FILE")
opt_parser.add_argument("-o", "--output_folder",dest="output", help="Insert an output directory to write to", metavar="FILE")
opt_parser.add_argument("-f", "--output_file",dest="output_file", help="Insert an output filename", metavar="FILE")
opt_parser.add_argument("-s", "--sample_type",dest="sample_type", help="Mention which type of sample you have", metavar="FILE")

options = opt_parser.parse_args()


csvfile_name = options.csvfile_name
output = options.output
output_file = options.output_file
sample_type = str(options.sample_type)




def create_colored_bed(table_name = "", output = output, output_file = output_file, sample_type = ""):
    df = pd.read_csv(table_name, sep="\t",header=0)
    df.columns = ["Fragment","Start","End","Fragment","rel_n_Reads","Strand","Thikstart","Thikend"]
    df["rel_n_Reads"] = df["rel_n_Reads"] / df["rel_n_Reads"].sum()
    print(df)    

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
    df["Strand"] = ["." for i in colors_list]
    df["rel_n_Reads"] = df["rel_n_Reads"] * 1000000
    
    with open(output + output_file, 'w') as fp:
        fp.write(df.to_csv(sep="\t",header=False,index=False)) 


    most_abundant_df = df[df["rel_n_Reads"] >= 1000]

    with open(output + "most_abundant_" + output_file, 'w') as fp:
        fp.write(most_abundant_df.to_csv(sep="\t",header=False,index=False))
    print(most_abundant_df) 

    most_abundant_df = df[df["rel_n_Reads"] >= 1000]

    with open(output + "most_abundant_" + output_file, 'w') as fp:
        fp.write(most_abundant_df.to_csv(sep="\t",header=False,index=False))
    print(most_abundant_df) 

    top100_df = df.sort_values(by="rel_n_Reads",ascending=False)[0:min(100,df.shape[0])]
    with open(f"{output}top_100_{output_file}", "w") as fp:
        fp.write(top100_df.to_csv(sep="\t", header=False, index=False))

create_colored_bed(csvfile_name, output, output_file, sample_type)