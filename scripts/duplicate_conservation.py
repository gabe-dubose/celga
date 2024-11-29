#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
from scipy.stats import ks_2samp

##### NOTES #####
#The dominance index is calcualted as (Nmax/N)
#################

#load data
data = pd.read_csv('../data/duplication/homologous_group_phylo_data_tau.csv')
dpl_data = pd.read_csv('../data/dpl_tau_homologous_group_data.csv')

#calculate standard deviations
data['group.tau.sd'] = np.sqrt(data['group.tau.variance'])

#compare cele and dpl stage-bias conservation
sns.set_style('whitegrid')
fig, [ax1,ax2] = plt.subplots(1, 2, figsize=(8, 4))

#cele
sns.ecdfplot(data=data, x='group.tau.sd', ax=ax1, 
             linewidth=2, color='black')
#dpl
sns.ecdfplot(data=dpl_data, x='group.tau.stdev', ax=ax1,
             linestyle='--', linewidth=2, color='tab:gray')

#add legend
#define legend
legend_lines = [mlines.Line2D([0], [0], color='black', label=r'${C. elegans}$', linewidth=2),
                mlines.Line2D([0], [0], color='tab:gray', label=r'${D. plexippus}$', linewidth=2, linestyle='--'),]

ax1.legend(handles=legend_lines, loc='lower right', fontsize=10)
ax1.set_ylabel('Proportion', fontsize=15)
ax1.set_xlabel(r'$\tau$ st.dev.', fontsize=15)
ax1.tick_params(axis='x', labelsize=15)
ax1.tick_params(axis='y', labelsize=15)

#cele
sns.ecdfplot(data=data, x='group.stage.bias.frequence', ax=ax2, 
             linewidth=2, color='black')
#dpl
sns.ecdfplot(data=dpl_data, x='group.stage.bias.frequency', ax=ax2,
             linestyle='--', linewidth=2, color='tab:gray')

#add legend
ax2.legend(handles=legend_lines, loc='lower right', fontsize=10)
ax2.set_ylabel('Proportion', fontsize=15)
ax2.set_xlabel(r'Stage-bias uniformity ($N_{max}/N$)', fontsize=15)
ax2.tick_params(axis='x', labelsize=15)
ax2.tick_params(axis='y', labelsize=15)

ax1.text(-0.15, 1.1, 'A', transform=ax1.transAxes, fontsize=15, fontweight='bold', va='top', ha='right')
ax2.text(-0.15, 1.1, 'B', transform=ax2.transAxes, fontsize=15, fontweight='bold', va='top', ha='right')

plt.tight_layout()
#plt.savefig('../figures/cele_vs_dplex_stage_specificity_gene_duplication.pdf')

##### STATISTICS #####
ks_statistic, p_value = ks_2samp(list(data['group.tau.sd']), list(dpl_data['group.tau.stdev']))
print(f"Standard deviation: D = {ks_statistic}, p = {p_value}")
#St.dev distributions: D = 0.30627316556586687, p = 2.3885159008465274e-10


ks_statistic, p_value = ks_2samp(list(data['group.stage.bias.frequence']), list(dpl_data['group.stage.bias.frequency']))
print(f"Stage-bias uniformity: D = {ks_statistic}, p = {p_value}")
#Stage-bias uniformity: D = 0.35197725942646935, p = 1.3800602968666975e-13

######################

'''
sns.set_style('whitegrid')
fig, [ax1,ax2] = plt.subplots(1, 2, figsize=(8, 4))

sns.histplot(data=data, x='group.tau.sd', ax=ax1,
    stat="density", element="step", color='tab:gray', alpha=1)
sns.histplot(data=data, x='group.stage.bias.frequence', ax=ax2,
    stat="density", element="step", color='tab:gray', alpha=1)

ax1.set_ylabel('Density', fontsize=15)
ax1.set_xlabel(r'$\tau$ st.dev.', fontsize=15)
ax1.tick_params(axis='x', labelsize=15)
ax1.tick_params(axis='y', labelsize=15)

ax2.set_ylabel('Density', fontsize=15)
ax2.set_xlabel('Proportion of genes', fontsize=15)
ax2.tick_params(axis='x', labelsize=15)
ax2.tick_params(axis='y', labelsize=15)

plt.tight_layout()
plt.show()

sns.set_style('whitegrid')
fig, ax1 = plt.subplots(1, 1, figsize=(4, 4))

sns.regplot(data=data, x='group.tau.mean', y='group.tau.sd', order=2, ax=ax1)

'''
