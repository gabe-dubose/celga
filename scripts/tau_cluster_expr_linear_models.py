#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.lines as mlines

#load data
data = pd.read_csv('../data/tau_cluster_data.csv')
#partition data
E_data = data.loc[data['stage_bias'] == 'Egg']
L1_data = data.loc[data['stage_bias'] == 'L1']
D_data = data.loc[data['stage_bias'] == 'Dauer']
L4_data = data.loc[data['stage_bias'] == 'L4']
A_data = data.loc[data['stage_bias'] == 'Adult']

# Define colors
cmap = plt.cm.get_cmap('viridis')
num_bins = 21
colors = [cmap(i / num_bins) for i in range(num_bins)]
colors = [mcolors.to_hex(color) for color in colors]
colors = [colors[0], colors[4], colors[12], colors[16], colors[20]]

# Define linestyles
linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]

#make plot
# Set the style
sns.set_style('whitegrid')

# Create the plot
fig, [ax1,ax2,ax3,ax4,ax5] = plt.subplots(1, 5, figsize=(25, 3))

#Embryo
sns.scatterplot(data=E_data.loc[E_data['embryo_clustered'] == 0], x='tau', y='max_expression',
                ax=ax1, color='lightgray', edgecolor=None, alpha=0.5)
sns.scatterplot(data=E_data.loc[E_data['embryo_clustered'] == 1], x='tau', y='max_expression',
                ax=ax1, color=colors[0], edgecolor=None, alpha=0.1)

sns.regplot(data=E_data.loc[E_data['embryo_clustered'] == 1], x='tau', y='max_expression',
            ax=ax1, scatter=False, color=colors[0], line_kws={'linestyle': '-'}, ci=None)
sns.regplot(data=E_data.loc[E_data['embryo_clustered'] == 0], x='tau', y='max_expression',
            ax=ax1, scatter=False, color='tab:gray', line_kws={'linestyle': '--'}, ci=None)
#ax1.set_ylim(0,1)

#L1
sns.scatterplot(data=L1_data.loc[L1_data['L1_clustered'] == 0], x='tau', y='max_expression',
                ax=ax2, color='lightgray', edgecolor=None, alpha=0.5)
sns.scatterplot(data=L1_data.loc[L1_data['L1_clustered'] == 1], x='tau', y='max_expression',
                ax=ax2, color=colors[1], edgecolor=None, alpha=0.5)

sns.regplot(data=L1_data.loc[L1_data['L1_clustered'] == 1], x='tau', y='max_expression',
            ax=ax2, scatter=False, color=colors[1], line_kws={'linestyle': '-'}, ci=None)
sns.regplot(data=L1_data.loc[L1_data['L1_clustered'] == 0], x='tau', y='max_expression',
            ax=ax2, scatter=False, color='tab:gray', line_kws={'linestyle': '--'}, ci=None)

#Dauer
sns.scatterplot(data=D_data.loc[D_data['Dauer_clustered'] == 0], x='tau', y='max_expression',
                ax=ax3, color='lightgray', edgecolor=None, alpha=0.5)
sns.scatterplot(data=D_data.loc[D_data['Dauer_clustered'] == 1], x='tau', y='max_expression',
                ax=ax3, color=colors[2], edgecolor=None, alpha=0.5)

sns.regplot(data=D_data.loc[D_data['Dauer_clustered'] == 1], x='tau', y='max_expression',
            ax=ax3, scatter=False, color=colors[2], line_kws={'linestyle': '-'}, ci=None)
sns.regplot(data=D_data.loc[D_data['Dauer_clustered'] == 0], x='tau', y='max_expression',
            ax=ax3, scatter=False, color='tab:gray', line_kws={'linestyle': '--'}, ci=None)

#L4
sns.scatterplot(data=L4_data.loc[L4_data['L4_clustered'] == 0], x='tau', y='max_expression',
                ax=ax4, color='lightgray', edgecolor=None, alpha=0.5)
sns.scatterplot(data=L4_data.loc[L4_data['L4_clustered'] == 1], x='tau', y='max_expression',
                ax=ax4, color=colors[3], edgecolor=None, alpha=0.5)

sns.regplot(data=L4_data.loc[L4_data['L4_clustered'] == 1], x='tau', y='max_expression',
            ax=ax4, scatter=False, color=colors[3], line_kws={'linestyle': '-'}, ci=None)
sns.regplot(data=L4_data.loc[L4_data['L4_clustered'] == 0], x='tau', y='max_expression',
            ax=ax4, scatter=False, color=colors[3], line_kws={'linestyle': '--'}, ci=None)

#Adult
sns.scatterplot(data=A_data.loc[A_data['Adult_clustered'] == 0], x='tau', y='max_expression',
                ax=ax5, color='lightgray', edgecolor=None, alpha=0.5)
sns.scatterplot(data=A_data.loc[A_data['Adult_clustered'] == 1], x='tau', y='max_expression',
                ax=ax5, color=colors[4], edgecolor=None, alpha=0.5)

sns.regplot(data=A_data.loc[A_data['Adult_clustered'] == 1], x='tau', y='max_expression',
            ax=ax5, scatter=False, color=colors[4], line_kws={'linestyle': '-'}, ci=None)
sns.regplot(data=A_data.loc[A_data['Adult_clustered'] == 0], x='tau', y='max_expression',
            ax=ax5, scatter=False, color=colors[4], line_kws={'linestyle': '--'}, ci=None)

plt.tight_layout()
plt.show()
#plt.savefig('../figures/tau_expr_clust_linear_models.pdf')