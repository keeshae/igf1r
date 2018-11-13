# BNGL Editor will contain functions for reading, writing, and editing BNGL files

import pandas as pd
import math
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from scipy.stats import pearsonr

__author__ = 'keerickson'

# Font sizes
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 14

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# writes a new bngl file
def write_bngl(model, newfilename):
    with open(newfilename, 'w') as new:
        for item in model:
            new.write(item)


# makes n variants on a cell line model
# model is the generic bngl file
# std is the desired standard deviation in cell-to-cell heterogeneity
# cpc is a pandas df with new copy numbers (data) and cell lines (columns). First column is protein names.
def pop_variants(model,cpc,std,n):
    cell_lines = list(cpc.columns.values)  # makes a list of the cell lines

    # iterate over all cell lines
    for col in range(1, len(cpc.columns)):

        # modmodel = list(model)  # makes a copy of model for modification
        cell = cell_lines[col]  # Cell line name

        if not os.path.exists(cell):  # make a directory to store these cell line's population files
            os.makedirs(cell)

        # make n variants
        for i in range(0, n):
            modmodel = list(model)  # makes a copy of model for modification

            # iterate over all proteins
            for index, row in cpc.iterrows():
                loi = row['PROTEIN'] + '\tcopynumber\n'  # line of interest (string)

                #  Now find the location of the line of interest in the model
                target = model.index(loi)

                # Alter the copy number for the protein
                if row[cell] > 0:  # for all nonzero proteins
                    mean = math.log(
                        row[cell])  # average cpc is from the csv - log of the avg is the mean of the distribution
                    sample = np.random.normal(mean, std)  # draw sample from normal distribution
                    new_cpc = int(round(math.exp(sample)))
                else:
                    new_cpc = row[cell]  # keep it 0

                # Edit this location in the modmodel file with the copy number from the corresponding cell line
                modmodel[target] = modmodel[target].replace("copynumber", str(new_cpc) + '*f # [=] copies per cell')

            # Write the modified model to a new cell line-specific file in a directory
            write_bngl(modmodel, cell + '/' + cell + '-' + str(i) + '.bngl')

def pearson_heatmap(cell):
    # Import csv of correlation coefficients
    df = pd.read_csv(cell+'.csv', sep=',')

    # Get y-axis labels
    proteins = df['proteins'].tolist()
    df.drop(['proteins'], axis=1, inplace=True)  # get rid of this column so just numerical data is left

    # Get x-axis labels
    #proteins = list(df.columns.values)

    # Get minimum value for heatmap scale
    mincor = np.nanmin(df.values)
    if mincor > -0.2:
        mincor = -0.2

    # Make heatmap fig
    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.2, left=0.2, top=0.92, right=0.92)

    heatmap = ax.pcolor(df, edgecolors='w', linewidths=0.5, cmap='RdBu', vmin=mincor, vmax=-mincor)
    cbar = plt.colorbar(heatmap, format='%.2f')

    # Set ticks in center of cells
    ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)

    # Set number of ticks to display on x-axis
    #plt.locator_params(axis='x', nbins=6)  # use nticks for log scale

    # Set tick labels
    ax.set_xticklabels(proteins, rotation='vertical')
    ax.set_yticklabels(proteins)

    # label colorbar
    cbar.ax.set_ylabel("Pearson's r")

    # save plot
    plt.savefig(cell+'_correlation.png', dpi=600)

    # display plot
    plt.show()

# path is location of gdat files
# cell is name of the cell line
def pearson_calc(path, cell):
    results = pd.DataFrame([])  # empty dataframe for proteins
    results2 = pd.DataFrame([]) # empty dataframe for phosphosites

    proteins = ['IRS1', 'ABL2', 'NCK2',
                'PIK3R1', 'PIK3R3', 'VAV2',
                'PLCG2', 'SRC', 'STAT1', 'YES1',
                'RASA1', 'PIK3R2', 'SHC1', 'SYK',
                'BLK', 'CRKL', 'ITK', 'ZAP70']

    for counter, file in enumerate(glob.glob(path + cell + '*.gdat')):
        namedf = pd.read_csv(file, delim_whitespace=True, skiprows=range(1, 101),
                             usecols=[7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24])
        results = results.append(namedf) # protein recruitment (absolute)

        sitesdf = pd.read_csv(file, delim_whitespace=True, skiprows=range(1, 101), usecols=[1,2,3,4,5,6]) # phosphosite data
        results2 = results2.append(sitesdf)

    # rename the columns to something more manageable
    results.columns = proteins

    # get estimate for IGF1R copy number from average phosphosite data
    results2['IGF1R_guess'] = results2.mean(axis=1)/0.84

    # normalize the data in results by corresponding IGF1R_guess in results 2
    results = results.div(results2['IGF1R_guess'], axis=0) # divides each row (axis = 0) by corresponding value in results2 to get normalized recruitment 

    # write to file (for testing)
    #results.to_csv(cell+'_ss_recruitment.csv')
    #results2.to_csv(cell+'_ss_phosphorylation.csv')

    pearsons = pd.DataFrame([])
    pearsons['proteins'] = proteins  # first column has protein names

    # generate all pairwise Pearson r
    # first iterate through vectors in results df (columns)
    for p in proteins:
        p1 = list(results[p])
        corr_coeffs = []  # will hold a vector of correlation coefficients

        # then iterate through names in col 1 of pearsons df
        for index, row in pearsons.iterrows():
            p2 = list(results[row['proteins']])  # row['proteins'] has name of a protein in pearsons df

            # p1 and p2 are now vectors of ss recruitment results from gdat files for two proteins
            # calculate the Pearson's r between p1 and p2
            r, pval = pearsonr(p1, p2)

            corr_coeffs.append(r)

        pearsons[p] = corr_coeffs  # add vector of correlation coefficients as a new column in df

    # write to file
    pearsons.to_csv(cell + '.csv')

# just accumulate ss recruitment values in a csv
def pop_calc(path, cell):
    results = pd.DataFrame([])  # empty dataframe for proteins
    results2 = pd.DataFrame([]) # empty dataframe for phosphosites

    proteins = ['IRS1', 'ABL2', 'NCK2',
                'PIK3R1', 'PIK3R3', 'VAV2',
                'PLCG2', 'SRC', 'STAT1', 'YES1',
                'RASA1', 'PIK3R2', 'SHC1', 'SYK',
                'BLK', 'CRKL', 'ITK', 'ZAP70']

    for counter, file in enumerate(glob.glob(path + cell + '*.gdat')):
        namedf = pd.read_csv(file, delim_whitespace=True, skiprows=range(1, 101),
                             usecols=[7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24])
        results = results.append(namedf) # protein recruitment (absolute)

        sitesdf = pd.read_csv(file, delim_whitespace=True, skiprows=range(1, 101), usecols=[1,2,3,4,5,6]) # phosphosite data
        results2 = results2.append(sitesdf)

    # rename the columns to something more manageable
    results.columns = proteins

    # get estimate for IGF1R copy number from average phosphosite data
    results2['IGF1R_guess'] = results2.mean(axis=1)/0.84

    # normalize the data in results by corresponding IGF1R_guess in results 2
    results = results.div(results2['IGF1R_guess'], axis=0) # divides each row (axis = 0) by corresponding value in results2 to get normalized recruitment

    # write to file
    results.to_csv(cell+'_ss_recruitment.csv')