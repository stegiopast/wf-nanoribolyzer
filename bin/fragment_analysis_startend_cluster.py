#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from multiprocessing import Pool
from itertools import repeat
import seaborn as sb
import math
from sklearn.neighbors import NearestNeighbors
from matplotlib import pyplot as plt
import seaborn as sb
import math
from sklearn.cluster import HDBSCAN
import matplotlib
#%%
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




# opt_parser = argparse.ArgumentParser()

# opt_parser.add_argument("-i", "--input_file", dest="bedfile_name", help="Insert a sample bam file", metavar="FILE")
# opt_parser.add_argument("-o", "--output_path",dest="output", help="Insert an output directory to write to", metavar="FILE")
# opt_parser.add_argument("-c","--cores", dest= "cores", help="Insert number of cores",metavar="FILE")
# opt_parser.add_argument("-s", "--sample_type",dest="sample_type", help="Mention which type of sample you have", metavar="FILE")

# options = opt_parser.parse_args()

# bedfile_name = options.bedfile_name
# output = options.output
# sample_type = str(options.sample_type)
# cores = int(options.cores)

bedfile_name = "~/GiantDisk/Data_nano_ribolyzer/20230920_NucExp_IVPA_9bc_Cytoplasm1/basecalling_output/filtered.bed"

bed_dataframe = pd.read_csv(bedfile_name,sep="\t",header=None)
bed_dataframe.columns = ["Template","Start","End","Read_id","Score","Direction"]


print(bed_dataframe.Start[0:200])
print(bed_dataframe.End[0:200])


# %%
fig,ax = plt.subplots(1,1)
q = sb.scatterplot(data=bed_dataframe, x=bed_dataframe.Start,y=bed_dataframe.End, ax = ax, alpha = 0.1)
fig.show()

#%%
print(q.table)

# %%
bedfile_name = "~/GiantDisk/Data_nano_ribolyzer/20230920_NucExp_IVPA_9bc_Nucleus1/basecalling_output/filtered.bed"

bed_dataframe = pd.read_csv(bedfile_name,sep="\t",header=None)
bed_dataframe.columns = ["Template","Start","End","Read_id","Score","Direction"]


print(bed_dataframe.Start[0:200])
print(bed_dataframe.End[0:200])

# %%
fig,ax = plt.subplots(1,1)
q = sb.scatterplot(data=bed_dataframe, x=bed_dataframe.Start,y=bed_dataframe.End, ax = ax, alpha = 0.1)
#ax.set_xlim(4000,5500)
#ax.set_ylim(5050,5500)
fig.show()
# %%
#bedfile_name = "~/test_validation/validation.bed"
bedfile_name ="/home/stefan/GiantDisk/validation_datasets_NanoRibolyzer/dataset1/validation.bed"

bed_dataframe = pd.read_csv(bedfile_name,sep="\t",header=None)
bed_dataframe.columns = ["Template","Start","End","Read_id","Score","Direction"]


print(bed_dataframe.Start[0:200])
print(bed_dataframe.End[0:200])

# %%
#bedfile_name = "~/test_validation/validation.bed"
bedfile_name ="/home/stefan/GiantDisk/validation_datasets_NanoRibolyzer/dataset11/validation.bed"
validation_dataset_bedfile_name ="/home/stefan/GiantDisk/validation_datasets_NanoRibolyzer/dataset11/validation_dataset.bed"

bed_dataframe = pd.read_csv(bedfile_name,sep="\t",header=None)
bed_dataframe.columns = ["Template","Start","End","Read_id","Score","Direction"]

validation_bed_dataframe = pd.read_csv(validation_dataset_bedfile_name,sep="\t",header=None)
print(validation_bed_dataframe)
validation_bed_dataframe.columns = ["Template","Start","End","Read_id","Score","Direction","Thikstart","Thikend"]

print(bed_dataframe.Start[0:200])
print(bed_dataframe.End[0:200])


# %%
fig,ax = plt.subplots(1,1)
sb.scatterplot(data=bed_dataframe, x=bed_dataframe.Start,y=bed_dataframe.End, ax = ax, alpha = 0.1)#
sb.scatterplot(data=validation_bed_dataframe, x=validation_bed_dataframe.Start,y=validation_bed_dataframe.End, ax = ax, alpha = 1)
ax.set_xlim(0,13500)
ax.set_ylim(0,13500)
fig.show()



#%%
bed_dataframe

intensity_matrix = np.zeros((max(bed_dataframe.End)+1,max(bed_dataframe.End)+1),dtype=int)
id_dict = {}


intensity_matrix

for x,y,id in zip(bed_dataframe.Start,bed_dataframe.End,bed_dataframe.Read_id):
    if f"{x}:{y}" not in id_dict.keys():
        id_dict[f"{x}:{y}"] = []
    intensity_matrix[x][y] += 1
    id_dict[f"{x}:{y}"].append(id)
intensity_matrix.sum()
id_dict.keys()

intensity_matrix[intensity_matrix <= 1] = 0
#neighbors = NearestNeighbors(n_neighbors=2)
#neighbors_fit = neighbors.fit(intensity_matrix)
#distances, indices = neighbors_fit.kneighbors(intensity_matrix)

#ellbow = np.sort(distances, axis=0)
#ellbow = ellbow[:,1]
#print(np.mean(ellbow))
#mean_distance = np.mean(ellbow)

dbscan_dataframe = np.array([[x,y] for x,y in zip(*np.nonzero(intensity_matrix))])
dbscan_dataframe

dbscan_dataframe_intensity = pd.DataFrame(np.array([[x,y, intensity_matrix[x,y]] for x,y in zip(*np.nonzero(intensity_matrix))]),columns=["X","Y","Counts"])
print(dbscan_dataframe_intensity)
# %%

fig,ax = plt.subplots(1,1)
sb.scatterplot(x=[i[0] for i in dbscan_dataframe],y=[i[1] for i in dbscan_dataframe], ax = ax, alpha = 1)#
sb.scatterplot(data=validation_bed_dataframe, x=validation_bed_dataframe.Start,y=validation_bed_dataframe.End, ax = ax, alpha = 0.8)

#ax.set_xlim(0,12000)
#ax.set_ylim(12000,13000)
fig.show()

# %%
fig,ax = plt.subplots(1,1)
q = sb.scatterplot(data=dbscan_dataframe_intensity,x="X",y="Y",hue=dbscan_dataframe_intensity.Counts,size=dbscan_dataframe_intensity.Counts, ax = ax, alpha = 1)#
# ax.set_xlim(3000,7000)
# ax.set_ylim(3000,7000)
fig.show()

# %%
mean =dbscan_dataframe_intensity.Counts.mean()
std = dbscan_dataframe_intensity.Counts.std()
dbscan_dataframe_intensity_high = dbscan_dataframe_intensity[dbscan_dataframe_intensity.Counts > mean + std]
validation_bed_dataframe_high = validation_bed_dataframe[validation_bed_dataframe["Score"] > validation_bed_dataframe["Score"].mean() + validation_bed_dataframe["Score"].std()]
print(dbscan_dataframe_intensity_high)
print(validation_bed_dataframe)
fig,ax = plt.subplots(1,1)
sb.scatterplot(data=dbscan_dataframe_intensity_high,x="X",y="Y", ax = ax, alpha = 1)
sb.scatterplot(data=validation_bed_dataframe_high, x=validation_bed_dataframe_high.Start,y=validation_bed_dataframe_high.End, ax = ax, alpha = 0.5)
ax.set_xlim(0,13500)
ax.set_ylim(0,13500)
#ax.set_xlim(4000,6000)
#ax.set_ylim(4000,6000)
fig.show()
dbscan_dataframe_intensity_high.shape

#%%
print(dbscan_dataframe_intensity.Counts.mean())
print(dbscan_dataframe_intensity.Counts.std())
#%%
from sklearn.cluster import DBSCAN
result = DBSCAN(eps=1, min_samples=5).fit(dbscan_dataframe)
import matplotlib
labels = result.labels_ / np.max(result.labels_)
plt.scatter(dbscan_dataframe[:,0], dbscan_dataframe[:,1], c=[matplotlib.cm.twilight(float(i)) for i in labels])







# %%
from sklearn.cluster import HDBSCAN
result = HDBSCAN(min_cluster_size = 5,store_centers='centroid').fit(dbscan_dataframe)
import matplotlib
labels = result.labels_ / np.max(result.labels_)
plt.scatter(dbscan_dataframe[:,0], dbscan_dataframe[:,1], c=[matplotlib.cm.twilight(float(i)) for i in labels])
#ax.set_xlim(0,13500)
#ax.set_ylim(0,13500)
ax.set_xlim(4000,6000)
ax.set_ylim(4000,6000)



# %%
from sklearn.cluster import OPTICS
result = OPTICS(min_samples = 5).fit(dbscan_dataframe)
import matplotlib
labels = result.labels_ / np.max(result.labels_)

plt.scatter(dbscan_dataframe[:,0], dbscan_dataframe[:,1], c=[matplotlib.cm.twilight(float(i)) for i in labels])

# %%
fig,ax = plt.subplots(1,1)
ax.scatter(dbscan_dataframe[:,0], dbscan_dataframe[:,1], c=[matplotlib.cm.twilight(float(i)) for i in labels],alpha=0.2)
sb.scatterplot(data=validation_bed_dataframe_high, x=validation_bed_dataframe_high.Start,y=validation_bed_dataframe_high.End, ax = ax, alpha = 1)
#ax.set_xlim(0,13500)
#ax.set_ylim(0,13500)
ax.set_xlim(4000,6000)
ax.set_ylim(4000,6000)
# %%
fig,ax = plt.subplots(1,1)
ax.scatter(dbscan_dataframe[:,0], dbscan_dataframe[:,1], c=[matplotlib.cm.twilight(float(i)) for i in labels],alpha=0.2)
sb.scatterplot(data=validation_bed_dataframe_high, x=validation_bed_dataframe_high.Start,y=validation_bed_dataframe_high.End, ax = ax, alpha = 1)
#ax.set_xlim(0,13500)
#ax.set_ylim(0,13500)
ax.set_xlim(4000,6000)
ax.set_ylim(4000,6000)
# %%
import dask
# %%
import dask.dataframe as dd

# %%
fig,ax = plt.subplots(1,1)
q = sb.scatterplot(intensity_matrix, ax = ax, alpha = 0.1)#
#ax.set_xlim(0,12000)
#ax.set_ylim(12000,13000)
fig.show()
# %%
print(intensity_matrix)