# Pearson's correlation coefficients
# for indiv cell lines
# pairwise

from bngl_pop_editor import pearson_calc
from bngl_pop_editor import pearson_heatmap
from bngl_pop_editor import pop_calc

cell = 'A549'
path = 'C:/Users/Keesha/PycharmProjects/IGF1R/'+cell+'/'

# make matrix of correlations coefficients and write to csv
pearson_calc(path, cell)

# make heatmap.. might need to manually clean up csv file generated from pearson calc
#pearson_heatmap(cell)

# just get ss rankings
#pop_calc(path,cell)









