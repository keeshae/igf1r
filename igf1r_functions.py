# make cell line specific heatmaps
# cluster cell lines by end point recruitment profile
# identify highly ranked binding partners

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import fileinput
import os
from clusterDefs import cluster_heatmap

def clean_gdat(loc, filename):
    for line in fileinput.input(loc + filename, inplace=True):
        # Remove spaces from header line after # sign
        # without this step, all the columns in the df will be offset
        print(line.replace("#              ", ""), end='')

    # Import one gdat file
    # delimiter is whitespace
    df = pd.read_csv(loc + filename, sep='\s+')

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

    # want proteins on y axis and time on x axis, so will transpose the df
    df2 = df.transpose()

    return df2, cell, proteins

def heatmapsIGF1R(loc):
    # execute specified actions for all files in the directory
    for filename in os.listdir(loc):
        # clean and normalize the gdat file into df
        # cell contains the cell line name
        # proteins is the unsorted list of proteins
        df2, cell, proteins = clean_gdat(loc, filename)

        # Set x-axis labels for heatmap. We'll force these.
        times = [0, 10, 20, 30, 40, 50, 60]

        # Make heatmap
        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.25, left=0.25)  # room for labels

        # sort df2 by the values in the last column (end point)
        # so the most highly recruited proteins are on top
        df2sort = df2.sort_values(by=df2.columns[df2.columns.size-1])
        # get vector of sorted protein names for heatmap label
        sorted_proteins = df2sort.index.tolist()

        # Make heatmap
        heatmap = ax.pcolor(df2sort, cmap='OrRd', edgecolors='w', linewidths=0.5)
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
        # Add cell line name as title
        plt.title(cell)
        # Save plot
        plt.savefig(os.getcwd()+'/NCI60/heatmaps/timecourse_' + cell + '.png', dpi=600)
        # Display plot
        #plt.show()

        # Close plot
        plt.close()

#### Clustering analysis ####
def clusteringIGF1R(loc):
    # empty dataframe for gathering end point data
    endpt = pd.DataFrame(index=[])

    # execute specified actions for all files in the directory
    for filename in os.listdir(loc):
        # clean and normalize the gdat file into df
        # cell contains the cell line name
        # proteins is the unsorted list of proteins
        df2, cell, proteins = clean_gdat(loc, filename)

        # Add last column from df2 into new dataframe for clustering analysis
        # use df2 here so that protein order is the same in every cell line (df2 is unsorted)
        # Last column in df2 contains normalized end point recruitment profile for each cell line
        # columns in endpt df are labeled with the cell name
        endpt[cell] = df2[df2.columns[df2.columns.size-1]]

    # polish dataframe for input to clustering analysis
    # I actually want proteins as the column names and cell lines as the rows
    endpt = endpt.transpose()
    # Name the first column 'CellLine'
    endpt.index.rename('CellLine', inplace=True)

    # write endptdf to csv
    endpt.to_csv('endpt_recruitment_data.csv')

    # make a list of the cell lines
    cell_lines = endpt.index.tolist()

    # cluster cell lines
    cluster_heatmap(endpt, cell_lines)

#### Rank analysis ####
def rankIGF1R(loc):
    # empty dataframe for gathering end point data
    rank = pd.DataFrame(index=[])

    # execute specified actions for all files in the directory
    for filename in os.listdir(loc):
        # clean and normalize the gdat file into df
        # cell contains the cell line name
        # proteins is the unsorted list of proteins
        df2, cell, proteins = clean_gdat(loc, filename)

        # Add new column to df2 for rank of protein recruitment,
        # obtained from end pt column (df2.columns.size -1)
        df2['rank'] = df2[df2.columns.size-1].rank(ascending=False)

        # Add rank to new df, label column with cell type
        rank[cell] = df2['rank']

    # Add a column to rank df with the average rank across all cell lines
    rank['average'] = rank.mean(axis=1)
    rank['stdev'] = rank.std(axis=1)
    rank = rank.sort_values(by=['average'])

    # plot the average/std of ranks
    fig, ax = plt.subplots(figsize=(14,6))
    fig.subplots_adjust(top=0.9, left=0.11)  # room for labels
    plt.errorbar(rank.index.tolist(), rank['average'], yerr=rank['stdev'], fmt='o')
    ax.set_ylabel('Average recruitment rank across cell lines')
    ax.set_xlabel('IGF1R binding partner')
    plt.locator_params(axis='y', nbins=rank.index.size)  # use nticks for log scale
    plt.show()

    # write ranks to csv
    rank.to_csv('ranked_recruitment_data.csv')