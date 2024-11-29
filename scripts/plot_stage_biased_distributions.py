#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
import numpy as np

data = pd.read_csv('../data/significantly_stage_biased_chromosomal_coords.csv')
E_biased = data.loc[data['source_variable'] == 'E_biased']
L1_biased = data.loc[data['source_variable'] == 'L1_biased']
L2_biased = data.loc[data['source_variable'] == 'L2_biased']
L3_biased = data.loc[data['source_variable'] == 'L3_biased']
D_biased = data.loc[data['source_variable'] == 'D_biased']
L4_biased = data.loc[data['source_variable'] == 'L4_biased']
A_biased = data.loc[data['source_variable'] == 'A_biased']

sns.set_style('whitegrid')
fig, [ax1,ax2,ax3,ax5,ax6,ax7] = plt.subplots(1, 6, figsize=(15, 3), sharey=True)

# Define colors
cmap = plt.cm.get_cmap('viridis')
num_bins = 21
colors = [cmap(i / num_bins) for i in range(num_bins)]
colors = [mcolors.to_hex(color) for color in colors]
colors = [colors[0], colors[4], colors[5], colors[8], colors[12], colors[16], colors[20]]

sns.histplot(data=E_biased, x='density', ax=ax1, color=colors[0], binwidth=1, alpha=1, discrete=True, element="step")
sns.histplot(data=L1_biased, x='density', ax=ax2, color=colors[1], binwidth=1, alpha=1, discrete=True, element="step")
sns.histplot(data=L2_biased, x='density', ax=ax3, color=colors[2], alpha=1, discrete=True, element="step")
#there are no L3 biased regions
sns.histplot(data=D_biased, x='density', ax=ax5, color=colors[4], binwidth=1, alpha=1, discrete=True, element="step")
sns.histplot(data=L4_biased, x='density', ax=ax6, color=colors[5], binwidth=1, alpha=1, discrete=True, element="step")
sns.histplot(data=A_biased, x='density', ax=ax7, color=colors[6], alpha=1, bins=[0, 1, 2, 3, 4], discrete=True, element="step")

ax1.set_title('Embryo-biased', fontsize=15)
ax2.set_title('L1-biased', fontsize=15)
ax3.set_title('L2-biased', fontsize=15)
ax5.set_title('Dauer-biased', fontsize=15)
ax6.set_title('L4-biased', fontsize=15)
ax7.set_title('Adult-biased', fontsize=15)

for ax in [ax1,ax2,ax3,ax5,ax6,ax7]:
    ax.set_ylabel("Count", fontsize=15)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_xlabel('Gene density', fontsize=15)

ax1.set_xticks([0,5,10,15,20,25,30,35])
ax2.set_xticks([0,2,4,6,8,10])
ax3.set_xticks([0,1,2,3])
ax5.set_xticks([0,2,4,6,8,10])
ax6.set_xticks([0,2,4,6,8])
ax7.set_xticks([0,1,2,3,4])

plt.tight_layout()
plt.savefig('../figures/stage_biased_regions_distributions.pdf')