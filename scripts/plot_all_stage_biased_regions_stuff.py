#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
import numpy as np

##### DATA #####
##### initial data #####
#load data stage biased distributions data
data = pd.read_csv('../data/significantly_stage_biased_chromosomal_coords.csv')
E_biased = data.loc[data['source_variable'] == 'E_biased']
L1_biased = data.loc[data['source_variable'] == 'L1_biased']
L2_biased = data.loc[data['source_variable'] == 'L2_biased']
L3_biased = data.loc[data['source_variable'] == 'L3_biased']
D_biased = data.loc[data['source_variable'] == 'D_biased']
L4_biased = data.loc[data['source_variable'] == 'L4_biased']
A_biased = data.loc[data['source_variable'] == 'A_biased']
#########################

##### stage biased frequency data #####
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


#################


##### PLOT #####
sns.set_style('whitegrid')
fig, [[ax1,ax2,ax3,ax5,ax6,ax7], [ax8,ax9,ax10,ax11,ax12,ax13], [ax14,ax15,ax16,ax17,ax18,ax19]] = plt.subplots(3, 6, figsize=(15, 7))

# Define colors
cmap = plt.cm.get_cmap('viridis')
num_bins = 21
colors = [cmap(i / num_bins) for i in range(num_bins)]
colors = [mcolors.to_hex(color) for color in colors]
colors = [colors[0], colors[4], colors[5], colors[8], colors[12], colors[16], colors[20]]

##### plot observed data #####
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
    if ax == ax1:
        ax.set_ylabel("Count", fontsize=15)
    else:
        ax.set_ylabel("")
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_xlabel('Gene density', fontsize=15)
    ax.set_ylim([0,7.2])

ax1.set_xticks([0,10,20,30])
ax2.set_xticks([0,2,4,6,8,10])
ax3.set_xticks([0,1,2,3])
ax5.set_xticks([0,2,4,6,8,10])
ax6.set_xticks([0,2,4,6,8])
ax7.set_xticks([0,1,2,3,4])
##################

##### plot null frequency #####

#Embryo biased
sns.kdeplot(null_dict['E_biased']['n_clusters'], ax=ax8, color=colors[0], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax8.axvline(x=E_n_clusters, color=colors[0], linestyle='--', linewidth=2)

#L1 biased
sns.kdeplot(null_dict['L1_biased']['n_clusters'], ax=ax9, color=colors[1], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax9.axvline(x=L1_n_clusters, color=colors[1], linestyle='--', linewidth=2)

#L2 biased
sns.kdeplot(null_dict['L2_biased']['n_clusters'], ax=ax10, color=colors[2], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax10.axvline(x=L2_n_clusters, color=colors[2], linestyle='--', linewidth=2)

#L3 biased
#sns.kdeplot(null_dict['L3_biased']['n_clusters'], ax=ax4, color=colors[3], 
#            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
#ax4.axvline(x=0, color=colors[3], linestyle='--', linewidth=2)
#note that L3 is excluded because there were no L3-biased regions

#D biased
sns.kdeplot(null_dict['D_biased']['n_clusters'], ax=ax11, color=colors[4], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax11.axvline(x=D_n_clusters, color=colors[4], linestyle='--', linewidth=2)

#L4 biased
sns.kdeplot(null_dict['L4_biased']['n_clusters'], ax=ax12, color=colors[5], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax12.axvline(x=L4_n_clusters, color=colors[5], linestyle='--', linewidth=2)

#A biased
sns.kdeplot(null_dict['A_biased']['n_clusters'], ax=ax13, color=colors[6], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax13.axvline(x=A_n_clusters, color=colors[6], linestyle='--', linewidth=2)

for ax in [ax8,ax9,ax10,ax11,ax12,ax13]:
    if ax == ax8:
        ax.set_ylabel("Density", fontsize=15)
    else:
        ax.set_ylabel("")
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_xlabel('n Stage-biased regions', fontsize=15)

ax8.set_title('Embryo: P(N$\geq$O) < 0.001', fontsize=14, multialignment='center')
ax9.set_title('L1: P(N$\geq$O) < 0.001', fontsize=14, multialignment='center')
ax10.set_title('L2: P(N$\geq$O) = 0.051', fontsize=14, multialignment='center')
ax11.set_title('Dauer: P(N$\geq$O) = 0.047', fontsize=14, multialignment='center')
ax12.set_title('L4: P(N$\geq$O) = 0.573', fontsize=14, multialignment='center')
ax13.set_title('Adult: P(N$\geq$O) < 0.001', fontsize=14, multialignment='center')
#################

##### plot null gene density #####

#Embryo biased
sns.kdeplot(null_dict['E_biased']['mean_cluster_size'], ax=ax14, color=colors[0], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax14.axvline(x=E_cluster_size_mean, color=colors[0], linestyle='--', linewidth=2)

#L1 biased
sns.kdeplot(null_dict['L1_biased']['mean_cluster_size'], ax=ax15, color=colors[1], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax15.axvline(x=L1_cluster_size_mean, color=colors[1], linestyle='--', linewidth=2)

#L2 biased
sns.kdeplot(null_dict['L2_biased']['mean_cluster_size'], ax=ax16, color=colors[2], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax16.axvline(x=L2_cluster_size_mean, color=colors[2], linestyle='--', linewidth=2)

#L3 biased
#sns.kdeplot(null_dict['L3_biased']['mean_cluster_size'], ax=ax4, color=colors[3], 
#            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
#ax4.axvline(x=0, color=colors[3], linestyle='--', linewidth=2)
#note that L3 is excluded because there were no L3-biased regions

#D biased
sns.kdeplot(null_dict['D_biased']['mean_cluster_size'], ax=ax17, color=colors[4], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax17.axvline(x=D_cluster_size_mean, color=colors[4], linestyle='--', linewidth=2)

#L4 biased
sns.kdeplot(null_dict['L4_biased']['mean_cluster_size'], ax=ax18, color=colors[5], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax18.axvline(x=L4_cluster_size_mean, color=colors[5], linestyle='--', linewidth=2)

#A biased
sns.kdeplot(null_dict['A_biased']['mean_cluster_size'], ax=ax19, color=colors[6], 
            fill=True, alpha=0.5, clip=(0, None), common_norm=True)
ax19.axvline(x=A_cluster_size_mean, color=colors[6], linestyle='--', linewidth=2)

for ax in [ax14,ax15,ax16,ax17,ax18,ax19]:
    if ax == ax14:
        ax.set_ylabel("Density", fontsize=15)
    else:
        ax.set_ylabel("")
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_xlabel('Mean gene density', fontsize=15)

ax14.set_title('Embryo: P(N$\geq$O) = 0.003', fontsize=14, multialignment='center')
ax15.set_title('L1: P(N$\geq$O) = 0.005', fontsize=14, multialignment='center')
ax16.set_title('L2: P(N$\geq$O) = 1', fontsize=14, multialignment='center')
ax17.set_title('Dauer: P(N$\geq$O) < 0.001', fontsize=14, multialignment='center')
ax18.set_title('L4: P(N$\geq$O) = 0.004', fontsize=14, multialignment='center')
ax19.set_title('Adult: P(N$\geq$O) = 0.058', fontsize=14, multialignment='center')

##################################

#add labels
ax1.text(-0.255, 1.35, 'A', transform=ax1.transAxes, fontsize=15, fontweight='bold', va='top', ha='right')
ax8.text(-0.255, 1.35, 'B', transform=ax8.transAxes, fontsize=15, fontweight='bold', va='top', ha='right')
ax14.text(-0.255, 1.35, 'C', transform=ax14.transAxes, fontsize=15, fontweight='bold', va='top', ha='right')

plt.tight_layout()
plt.savefig('../figures/stage_biased_regions_totals.pdf', bbox_inches='tight')