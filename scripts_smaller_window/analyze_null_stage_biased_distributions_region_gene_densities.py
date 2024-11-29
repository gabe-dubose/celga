#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
import warnings
warnings.filterwarnings("ignore")

#load observed stage biased chromosomal distributions
observed_data = pd.read_csv('../data/significantly_stage_biased_chromosomal_coords.csv')

#compute mean observed cluster sizes and frequencies for each stage

#partition stages
E_biased = observed_data.loc[observed_data['source_variable'] == 'E_biased']
L1_biased = observed_data.loc[observed_data['source_variable'] == 'L1_biased']
L2_biased = observed_data.loc[observed_data['source_variable'] == 'L2_biased']
#There are no L3-biased regions
D_biased = observed_data.loc[observed_data['source_variable'] == 'D_biased']
L4_biased = observed_data.loc[observed_data['source_variable'] == 'L4_biased']
A_biased = observed_data.loc[observed_data['source_variable'] == 'A_biased']

#get average cluster size and numebr of clusters for each stage
#Embryo
E_cluster_size_mean = E_biased['density'].mean()
E_n_clusters = len(E_biased)
#L1
L1_cluster_size_mean = L1_biased['density'].mean()
L1_n_clusters = len(L1_biased)
#L2
L2_cluster_size_mean = L2_biased['density'].mean()
L2_n_clusters = len(L2_biased)
#L3 - since there are no L3 regions
L3_cluster_size_mean = 0
L3_n_clusters = 0
#Dauer
D_cluster_size_mean = D_biased['density'].mean()
D_n_clusters = len(D_biased)
#L4
L4_cluster_size_mean = L4_biased['density'].mean()
L4_n_clusters = len(L4_biased)
#Adult
A_cluster_size_mean = A_biased['density'].mean()
A_n_clusters = len(A_biased)

#load null simulation data
null_data = pd.read_csv('../data/null_gene_cluster_distributions.csv')

#initialize dictionary to store null info
null_dict = {'E_biased' : {'mean_cluster_size' : [], 'n_clusters' : []},
             'L1_biased' : {'mean_cluster_size' : [], 'n_clusters' : []},
             'L2_biased' : {'mean_cluster_size' : [], 'n_clusters' : []},
             'L3_biased' : {'mean_cluster_size' : [], 'n_clusters' : []},
             'D_biased' : {'mean_cluster_size' : [], 'n_clusters' : []},
             'L4_biased' : {'mean_cluster_size' : [], 'n_clusters' : []},
             'A_biased' : {'mean_cluster_size' : [], 'n_clusters' : []}}


#parse null data
#define stage labels
stages = list(null_data['source_variable'].unique())

#iterate through each null simulation
for i in range(1000):
    iteration_data = null_data.loc[null_data['iteration'] == i]

    #iterate through each life stage
    for stage in stages:
        #get info
        stage_data = iteration_data.loc[iteration_data['source_variable'] == stage]
        #calculate mean size and n clusters
        stage_mean = stage_data['density'].mean()
        stage_n = len(stage_data)
        #add to dictionary
        if str(stage_mean) != 'nan':
            null_dict[stage]['mean_cluster_size'].append(stage_mean)
        null_dict[stage]['n_clusters'].append(stage_n)

#compare observed to null
#mean gene density
E_p = 0
for value in null_dict['E_biased']['mean_cluster_size']:
    if value >= E_cluster_size_mean:
        E_p += 1
null_count = len(null_dict['E_biased']['mean_cluster_size'])
print(f"E: Full null, P(N>=O) ={E_p/1000}")
print(f'E: Non 0 null, P(N>=O) = {E_p/null_count}')

L1_p = 0
for value in null_dict['L1_biased']['mean_cluster_size']:
    if value >= L1_cluster_size_mean:
        L1_p += 1
null_count = len(null_dict['L1_biased']['mean_cluster_size'])
print(f"L1: Full null, P(N>=O) ={L1_p/1000}")
print(f'L1: Non 0 null, P(N>=O) = {L1_p/null_count}')

L2_p = 0
for value in null_dict['L2_biased']['mean_cluster_size']:
    if value >= L2_cluster_size_mean:
        L2_p += 1
null_count = len(null_dict['L2_biased']['mean_cluster_size'])
print(f"L2: Full null, P(N>=O) ={L2_p/1000}")
print(f'L2: Non 0 null, P(N>=O) = {L2_p/null_count}')

L3_p = 0
for value in null_dict['L3_biased']['mean_cluster_size']:
    if value >= L3_cluster_size_mean:
        L3_p += 1
null_count = len(null_dict['L3_biased']['mean_cluster_size'])
print(f"L3: Full null, P(N>=O) ={L3_p/1000}")
#print(f'L3: Non 0 null, P(N>=O) = {L3_p/null_count}')

D_p = 0
for value in null_dict['D_biased']['mean_cluster_size']:
    if value >= D_cluster_size_mean:
        D_p += 1
null_count = len(null_dict['D_biased']['mean_cluster_size'])
print(f"D: Full null, P(N>=O) ={D_p/1000}")
print(f'D: Non 0 null, P(N>=O) = {D_p/null_count}')

L4_p = 0
for value in null_dict['L4_biased']['mean_cluster_size']:
    if value >= L4_cluster_size_mean:
        L4_p += 1
null_count = len(null_dict['L4_biased']['mean_cluster_size'])
print(f"L4: Full null, P(N>=O) ={L4_p/1000}")
print(f'L4: Non 0 null, P(N>=O) = {L4_p/null_count}')

A_p = 0
for value in null_dict['A_biased']['mean_cluster_size']:
    if value >= A_cluster_size_mean:
        A_p += 1
null_count = len(null_dict['A_biased']['mean_cluster_size'])
print(f"A: Full null, P(N>=O) ={A_p/1000}")
print(f'A: Non 0 null, P(N>=O) = {A_p/null_count}')

##### Null/Observed comparisons results ######
# I'm recording these here just for future reference
#E: Full null, P(N>=O) =0.0
#E: Non 0 null, P(N>=O) = 0.0
#L1: Full null, P(N>=O) =0.0
#L1: Non 0 null, P(N>=O) = 0.0
#L2: Full null, P(N>=O) =0.0
#L2: Non 0 null, P(N>=O) = 0.0
#L3: Full null, P(N>=O) =0.0
#D: Full null, P(N>=O) =0.03
#D: Non 0 null, P(N>=O) = 0.04178272980501393
#L4: Full null, P(N>=O) =0.0
#L4: Non 0 null, P(N>=O) = 0.0
#A: Full null, P(N>=O) =0.0
#A: Non 0 null, P(N>=O) = 0.0

##############################################

#plot null distributions
sns.set_style('whitegrid')
fig, [ax1,ax2,ax3,ax4,ax5,ax6,ax7] = plt.subplots(1, 7, figsize=(18, 3))

# Define colors
cmap = plt.cm.get_cmap('viridis')
num_bins = 21
colors = [cmap(i / num_bins) for i in range(num_bins)]
colors = [mcolors.to_hex(color) for color in colors]
colors = [colors[0], colors[4], colors[5], colors[8], colors[12], colors[16], colors[20]]

#Embryo biased
sns.kdeplot(null_dict['E_biased']['mean_cluster_size'], ax=ax1, color=colors[0], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax1.axvline(x=E_cluster_size_mean, color=colors[0], linestyle='--', linewidth=2)

#L1 biased
sns.kdeplot(null_dict['L1_biased']['mean_cluster_size'], ax=ax2, color=colors[1], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax2.axvline(x=L1_cluster_size_mean, color=colors[1], linestyle='--', linewidth=2)

#L2 biased
sns.kdeplot(null_dict['L2_biased']['mean_cluster_size'], ax=ax3, color=colors[2], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax3.axvline(x=L2_cluster_size_mean, color=colors[2], linestyle='--', linewidth=2)

#L3 biased
sns.kdeplot(null_dict['L3_biased']['mean_cluster_size'], ax=ax4, color=colors[3], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax4.axvline(x=0, color=colors[3], linestyle='--', linewidth=2)

#D biased
sns.kdeplot(null_dict['D_biased']['mean_cluster_size'], ax=ax5, color=colors[4], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax5.axvline(x=D_cluster_size_mean, color=colors[4], linestyle='--', linewidth=2)

#L4 biased
sns.kdeplot(null_dict['L4_biased']['mean_cluster_size'], ax=ax6, color=colors[5], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax6.axvline(x=L4_cluster_size_mean, color=colors[5], linestyle='--', linewidth=2)

#A biased
sns.kdeplot(null_dict['A_biased']['mean_cluster_size'], ax=ax7, color=colors[6], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax7.axvline(x=A_cluster_size_mean, color=colors[6], linestyle='--', linewidth=2)

for ax in [ax1,ax2,ax3,ax4,ax5,ax6,ax7]:
    if ax == ax1:
        ax.set_ylabel("Probability density", fontsize=15)
    else:
        ax.set_ylabel("")
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_xlabel('Mean gene density', fontsize=15)

ax1.set_title('Embryo: P(N$\geq$O) < 0.001', fontsize=14, multialignment='center')
ax2.set_title('L1: P(N$\geq$O) < 0.001', fontsize=14, multialignment='center')
ax3.set_title('L2: NA', fontsize=14, multialignment='center')
ax4.set_title('L3: NA', fontsize=14, multialignment='center')
ax5.set_title('Dauer: P(N$\geq$O) = 0.003', fontsize=14, multialignment='center')
ax6.set_title('L4: P(N$\geq$O) < 0.001', fontsize=14, multialignment='center')
ax7.set_title('Adult: P(N$\geq$O) < 0.001', fontsize=14, multialignment='center')

plt.tight_layout()
plt.savefig('../figures/cluster_regions_gene_densities_null_vs_obs.pdf')
