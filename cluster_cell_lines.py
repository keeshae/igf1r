# clustering the IGF1R recruitment
# this code clusters the cell lines if you call cluster_dendrogram
# clusters cell lines and proteins if you call cluster_heatmap
# modified from https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/

import pandas as pd
import seaborn as sns; sns.set(color_codes=True)
from clusterDefs import cluster_heatmap

# read in the csv
df = pd.read_csv('endpt_recruitment_data.csv', sep=',')

# make a list of the cell lines
cell_lines = df['CellLine'].tolist()
df.drop(['CellLine'], axis=1, inplace=True)  # leave just numerical data in df (first column was cell line strings)

# make a list of the proteins
proteins = list(df.columns.values)

### display heatmap with dendrogram
cluster_heatmap(df, cell_lines)





