# Make "population" model files
# n=5000 copies of each cell line bngl with slight modifications in protein copy number

from bngl_pop_editor import pop_variants
import pandas as pd

modelfile = 'generic_model.bngl'
copynums = 'all_copy_nums.csv'  # proteomics data for all cell lines

model = open(modelfile, 'r').readlines()  # read in generic bngl file

cpc = pd.read_csv(copynums, sep='\t')  # reads tab-separated csv file that contains new copy numbers into pandas df

# number of variants to make for each cell line
n = 3

# sdev in copy number
std = 0.2

# Make population bngl files
pop_variants(model, cpc, std, n)

