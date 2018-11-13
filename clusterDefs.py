from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import seaborn as sns; sns.set(color_codes=True)


###### Clustering with heatmap
# can change distance metric and linkage method here
def cluster_heatmap(dataframe, ylabels):
    g = sns.clustermap(dataframe, method='ward', metric='euclidean', yticklabels=ylabels, cmap="magma_r", linewidths=0.75, figsize=(7,10), row_cluster=True, col_cluster=False, cbar_kws={'label': 'Protein Bound : Total IGF1R'})
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.subplots_adjust(left=0.06, bottom=0.15, right=0.67, top=0.95, wspace=0, hspace=0)
    plt.show()
    #plt.savefig('clustermap.png', format='png', dpi=600)


##### hierarchical clustering (no heatmap shown)
# input data as numpy array
# methods to compute linkage => 'single', 'complete', 'average', 'ward'
# distance metrics => 'euclidean', 'cityblock', 'hamming', 'cosine'
def cluster_dendrogram(dataframe, leaf_labels):
    # convert pd df to numpy array
    D = dataframe.values

    Z = linkage(D, method='ward', metric='euclidean')

    # plot dendrogram
    plt.figure(figsize=(10, 6))
    ax = plt.subplot()
    plt.subplots_adjust(left=0.07, bottom=0.3, right=0.98, top=0.95, wspace=0, hspace=0)
    plt.xlabel('Cell Line')
    plt.ylabel('Distance')

    dendrogram(Z, leaf_rotation=90., leaf_font_size=10., labels=leaf_labels)

    plt.show()

