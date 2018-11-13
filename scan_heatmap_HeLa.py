import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Import gdat file
df = pd.read_csv("IGF1_scan_HeLa_S3_Feb27.csv", sep=',')

# Get y-axis labels
proteins = df['protein'].tolist()
df.drop(['protein'], axis=1, inplace = True) #get rid of this column so just numerical data is left

# Get x-axis labels
#times = list(df.columns.values)
conc = [r'$10^{-12}$', r'$10^{-11}$', r'$10^{-10}$', r'$10^{-9}$', r'$10^{-8}$', r'$10^{-7}$', r'$10^{-6}$']

# Make heatmap

fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.25, left=0.25) # room for labels

heatmap = ax.pcolor(df, edgecolors='w', linewidths=0.5)
cbar = plt.colorbar(heatmap, format='%.2f')

# Set ticks in center of cells
ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)
ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)

# Set number of ticks to display on x-axis
plt.locator_params(axis='x', nbins=7) # use nticks for log scale

# Set tick labels
ax.set_xticklabels(conc)
ax.set_yticklabels(proteins)

# Set axis labels
ax.set_xlabel('IGF1 Concentration (M)')
cbar.ax.set_ylabel('Protein Bound : Total IGF1R')


plt.savefig('scan_heatmap_Feb27.svg',dpi=600)
# Display plot
plt.show()
