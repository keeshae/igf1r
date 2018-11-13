# make cell line specific heatmaps
# cluster cell lines / proteins by end point recruitment profile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import fileinput
import os

# location of gdat files
# make sure there are only gdat files in here (rm cdat and net)
loc = 'C:/Users/Keesha/PycharmProjects/IGF1R/NCI60/bnglout/'

# empty dataframes for gathering endpt and rank data
endpt = pd.DataFrame(index=[])
rank = pd.DataFrame(index=[])

# execute for all files in the directory
for filename in os.listdir(loc):
    for line in fileinput.input(loc+filename, inplace=True):
        # Remove spaces from header line after # sign
        # without this step, all the columns in the df will be offset
        print(line.replace("#              ", ""), end='')

    # Import a gdat file
    # delimiter is whitespace
    df = pd.read_csv(loc+filename, sep='\s+')

    # Get cell name for labeling (remove ".gdat")
    cell = os.path.splitext(filename)[0]

    # Remove the columns that we do not need
    # keep only bound proteins, drop time column and bound phosphosites
    df = df.drop(df.columns[0:7], axis=1)

    # Get protein names from original order (this is before any sorting)
    # all columns are protein names
    proteins = df.columns.tolist()

    # Normalize proteins: present in ratio of protein bound to total IGF1R
    # IGF1R copy number is 14810 per cell for all models
    df = df/14810

    # Set x-axis labels. We'll force these.
    times = [0, 10, 20, 30, 40, 50, 60]

    # Make heatmap
    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.25, left=0.25)  # room for labels

    # want proteins on y axis and time on x axis, so will transpose the df
    df2 = df.transpose()

    # sort df2 by the values in the last column (end point)
    # so the most highly recruited proteins are on top
    df2sort = df2.sort_values(by=df2.columns[df2.columns.size-1])
    # get vector of sorted protein names
    sorted_proteins = df2sort.index.tolist()

    # Make heatmap
    heatmap = ax.pcolor(df2sort, edgecolors='w', linewidths=0.5)
    cbar = plt.colorbar(heatmap, format='%.2f')

    # Set ticks in center of cells
    ax.set_xticks(np.arange(df2sort.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(df2sort.shape[0]) + 0.5, minor=False)

    # Set number of ticks to display on x-axis
    plt.locator_params(axis='x', nbins=len(times))  # use nticks for log scale

    # Set tick labels
    ax.set_yticklabels(sorted_proteins)
    ax.set_xticklabels(times)

    # Set axis labels
    ax.set_xlabel('Time (min)')
    cbar.ax.set_ylabel('Protein Bound : Total IGF1R')

    # Save plot
    plt.savefig(os.getcwd()+'/NCI60/heatmaps/timecourse_' + cell + '.png', dpi=600)

    # Display plot
    #plt.show()

    # Close plot
    plt.close()

    #### Clustering analysis ####
    # Add last column from df2 into new dataframe for clustering analysis
    # use df2 here so that protein order is the same in every cell line (df2 is unsorted)
    # Last column in df2 contains normalized end point recruitment profile for each cell line
    # columns in endpt df are labeled with the cell name
    endpt[cell] = df2[df2.columns[df2.columns.size-1]]

    #### Rank analysis ####
    # Add new column to df2 for rank of protein recruitment,
    # obtained from end pt column (df2.columns.size -1)
    df2['rank'] = df2[df2.columns.size-1].rank(ascending=False)
    # Add rank to new df, label column with cell type
    rank[cell] = df2['rank']


# polish dataframe for input to clustering analysis
# I actually want proteins as the column names and cell lines as the rows
endpt = endpt.transpose()
# Name the first column 'CellLine'
endpt.index.rename('CellLine', inplace=True)

# write endptdf to csv
endpt.to_csv('endpt_recruitment_data.csv')

#### fun stuff with ranking binding partners
# Add a column with the average rank across all cell lines
rank['average'] = rank.mean(axis=1)
rank['stdev'] = rank.std(axis=1)
rank = rank.sort_values(by=['average'])

# plot the ranks
fig, ax = plt.subplots(figsize=(14,6))
fig.subplots_adjust(top=0.9, left=0.11)  # room for labels
plt.errorbar(rank.index.tolist(), rank['average'], yerr=rank['stdev'], fmt='o')
ax.set_ylabel('Average recruitment rank across cell lines')
ax.set_xlabel('IGF1R binding partner')
plt.locator_params(axis='y', nbins=rank.index.size)  # use nticks for log scale
plt.show()

# write ranks to csv
rank.to_csv('ranked_recruitment_data.csv')