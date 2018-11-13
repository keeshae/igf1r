# BNGL Editor will contain functions for reading, writing, and editing BNGL files

import pandas as pd
import os

# writes a new bngl file
def write_bngl(model, newfilename):
    with open(newfilename, 'w') as new:
        for item in model:
            new.write(item)


# changes the copy numbers of proteins in the bngl file according to newcpc
def change_cpc(model, newcpc, loc):

    cpc = pd.read_csv(newcpc, sep='\t')  # reads tab-separated csv file that contains new copy numbers into pandas df
    cell_lines = list(cpc.columns.values)  # makes a list of the cell lines

    # iterate over all cell lines
    for col in range(1, len(cpc.columns)):

        modmodel = list(model)  # makes a copy of model for modification
        cell = cell_lines[col]  # Cell line name

        # iterate over all proteins
        for index, row in cpc.iterrows():

            loi = row['PROTEIN']+'\tcopynumber\n'  # line of interest (string)

            # Now find the location of the line of interest in the model
            target = model.index(loi)

            # Edit this location in the modmodel file with the copy number from the corresponding cell line
            modmodel[target] = modmodel[target].replace("copynumber", str(row[cell])+'*f # [=] copies per cell')

        # Write the modified model to a new cell line-specific file
        os.makedirs(os.path.dirname(loc), exist_ok=True)
        write_bngl(modmodel, loc + cell + '.bngl')
