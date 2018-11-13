# Make cell line specific bngl file for the IGF1R model

from bngl_editor import change_cpc

modelfile = 'generic_model.bngl'  # template
copynums = 'all_copy_nums.csv'  # proteomics data for all cell lines
loc = 'C:/Users/Keesha/PycharmProjects/IGF1R/NCI60/'  # where we want the model files

model = open(modelfile, 'r').readlines()
change_cpc(model, copynums, loc)
