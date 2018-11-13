# Keesha Erickson, Nov 2018
# plot heatmaps from bngl files
# plot rankings of recruitment
# clustering by cell line and protein recruitment

from heatmaps_demo_def import heatmapsIGF1R, clusteringIGF1R, rankIGF1R

# location of gdat files
# there can only be gdat files in here (rm cdat and net)
loc = 'C:/Users/Keesha/PycharmProjects/IGF1R/NCI60/bnglout/'

# plot heatmaps of protein recruitment
heatmapsIGF1R(loc)

# rank analysis
rankIGF1R(loc)

# cluster cell lines
clusteringIGF1R(loc)

