#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.lines as mlines

#load data
tau_dataframe = pd.read_csv('../data/stage_specificity_data.csv')
#adjust stage bias label
tau_dataframe['stage_bias'] = tau_dataframe['stage_bias'].str.replace(r'^e.*', 'Egg', regex=True)
tau_dataframe['stage_bias'] = tau_dataframe['stage_bias'].str.replace(r'^dauer.*', 'Dauer', regex=True)
tau_dataframe['stage_bias'] = tau_dataframe['stage_bias'].str.replace(r'^adult.*', 'Adult', regex=True)

#load genomic position data
genome_data = pd.read_csv('../data/gtf_parsed.csv')

#remove mitochondrial (NC_001328.1)
genome_data = genome_data.loc[genome_data['chromosome'] != 'NC_001328.1']

#adjust gene id
genome_data['gene'] = genome_data['gene'].str.replace(r'^CELE_', '', regex=True)

#separate by chromosome
#chromosome I = NC_003279.8
chromosome_I = genome_data.loc[genome_data['chromosome'] == 'NC_003279.8']
chromosome_I = pd.merge(chromosome_I, tau_dataframe, left_on='gene', right_on='gene_id', how='inner')
chromosome_I = chromosome_I.sort_values(by='start', ascending=True)
chromosome_I['gene_position'] = range(1, len(chromosome_I) + 1)

#chromosome II = NC_003280.10
chromosome_II = genome_data.loc[genome_data['chromosome'] == 'NC_003280.10']
chromosome_II = pd.merge(chromosome_II, tau_dataframe, left_on='gene', right_on='gene_id', how='inner')
chromosome_II = chromosome_II.sort_values(by='start', ascending=True)
chromosome_II['gene_position'] = range(1, len(chromosome_II) + 1)

#chromosome III = NC_003281.10
chromosome_III = genome_data.loc[genome_data['chromosome'] == 'NC_003281.10']
chromosome_III = pd.merge(chromosome_III, tau_dataframe, left_on='gene', right_on='gene_id', how='inner')
chromosome_III = chromosome_III.sort_values(by='start', ascending=True)
chromosome_III['gene_position'] = range(1, len(chromosome_III) + 1)

#chromosome IV = NC_003282.8
chromosome_IV = genome_data.loc[genome_data['chromosome'] == 'NC_003282.8']
chromosome_IV = pd.merge(chromosome_IV, tau_dataframe, left_on='gene', right_on='gene_id', how='inner')
chromosome_IV = chromosome_IV.sort_values(by='start', ascending=True)
chromosome_IV['gene_position'] = range(1, len(chromosome_IV) + 1)

#chromosome V = NC_003283.11
chromosome_V = genome_data.loc[genome_data['chromosome'] == 'NC_003283.11']
chromosome_V = pd.merge(chromosome_V, tau_dataframe, left_on='gene', right_on='gene_id', how='inner')
chromosome_V = chromosome_V.sort_values(by='start', ascending=True)
chromosome_V['gene_position'] = range(1, len(chromosome_V) + 1)

#chromosome X = NC_003284.9
chromosome_X = genome_data.loc[genome_data['chromosome'] == 'NC_003284.9']
chromosome_X = pd.merge(chromosome_X, tau_dataframe, left_on='gene', right_on='gene_id', how='inner')
chromosome_X = chromosome_X.sort_values(by='start', ascending=True)
chromosome_X['gene_position'] = range(1, len(chromosome_X) + 1)

#A function to compute Cohen's d. This will be used in identifying regions on chromosomes where
#the density and degree of stage-biased expression significantly exceeds what is expected 
#of any given region randomly. 
# Function to compute Cohen's d
def cohens_d(x, y):
    diff_mean = np.mean(x) - np.mean(y)
    pooled_std = np.sqrt(((len(x) - 1) * np.var(x) + (len(y) - 1) * np.var(y)) / (len(x) + len(y) - 2))
    return diff_mean / pooled_std

# A function to run a sliding window that identifies regions where stage-biased genes
# are more densely clustered than expected randomly. 
# The algorithm is as follows:
# 1) Chromosome coordinates are iterated through, in a sliding window of some size
# 2) For each window, a random sample of 1000 other taue values is taken. This serves
#    as the baseline expectation for if tau values were organized randomly
# 3) Only regions where there are at least three genes are considered
# 4) The standardized effect size (Cohen's d) is calculated between the observed 
#    tau values for the windw and the random samples. This quantifies the deviation 
#    of the tau values in said window from the null expectation of random organization.

def sliding_window_expression_specificity(data, interval):
    start_min = 1
    start_max = data['start'].max()
    
    # Initialize an empty list to store the results
    results = []

    # Iterate through windows
    for start_pos in range(start_min, start_max + 1, interval):
        end_pos = start_pos + interval - 1
        # Select the window
        window = data[(data['start'] >= start_pos) & (data['start'] <= end_pos)]

        # Calculate average tau
        mean_tau = window['tau'].mean()

        # Calculate maximum expression
        mean_max_expression = window['max_expression'].mean()

        # Add density 
        density = len(window)
        
        # Randomly sample 1000 tau values from the entire dataset
        random_sample = data['tau'].sample(1000, replace=True)
        
        # Compute Cohen's d between the window tau values and the random sample
        if len(window) >= 3:  # Ensure there is at least 3
            observed_sample = window['tau']
            cohens_d_value = cohens_d(observed_sample, random_sample)
        else:
            cohens_d_value = np.nan  # Not enough data to compute Cohen's d
        
        # Append the results for this window
        results.append({
            'window_start': start_pos,
            'window_stop': end_pos,
            'average_tau': mean_tau,
            'max_expression': mean_max_expression,
            'density': density,
            'cohens_d': cohens_d_value
        })

    # Convert the results into a DataFrame
    result_df = pd.DataFrame(results)
    return result_df

#get sliding window data for egg
E_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'Egg'], 10000)
E_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'Egg'], 10000)
E_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'Egg'], 10000)
E_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'Egg'], 10000)
E_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'Egg'], 10000)
E_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'Egg'], 10000)

#get sliding window data for L1
L1_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'L1'], 10000)
L1_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'L1'], 10000)
L1_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'L1'], 10000)
L1_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'L1'], 10000)
L1_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'L1'], 10000)
L1_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'L1'], 10000)

#get sliding window data for L2
L2_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'L2'], 10000)
L2_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'L2'], 10000)
L2_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'L2'], 10000)
L2_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'L2'], 10000)
L2_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'L2'], 10000)
L2_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'L2'], 10000)

#get sliding window data for L3
L3_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'L3'], 10000)
L3_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'L3'], 10000)
L3_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'L3'], 10000)
L3_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'L3'], 10000)
L3_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'L3'], 10000)
L3_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'L3'], 10000)

#get sliding window data for dauer
D_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'Dauer'], 10000)
D_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'Dauer'], 10000)
D_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'Dauer'], 10000)
D_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'Dauer'], 10000)
D_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'Dauer'], 10000)
D_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'Dauer'], 10000)

#get sliding window data for L4
L4_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'L4'], 10000)
L4_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'L4'], 10000)
L4_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'L4'], 10000)
L4_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'L4'], 10000)
L4_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'L4'], 10000)
L4_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'L4'], 10000)

#get sliding window data for Adults
A_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'Adult'], 10000)
A_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'Adult'], 10000)
A_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'Adult'], 10000)
A_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'Adult'], 10000)
A_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'Adult'], 10000)
A_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'Adult'], 10000)

#combine dataframes
# Create a dictionary of dataframes with variable names as keys
dataframes = {
    'E_chrI_tau': E_chrI_tau, 'E_chrII_tau': E_chrII_tau, 'E_chrIII_tau': E_chrIII_tau, 
    'E_chrIV_tau': E_chrIV_tau, 'E_chrV_tau': E_chrV_tau, 'E_chrX_tau': E_chrX_tau,
    'L1_chrI_tau': L1_chrI_tau, 'L1_chrII_tau': L1_chrII_tau, 'L1_chrIII_tau': L1_chrIII_tau, 
    'L1_chrIV_tau': L1_chrIV_tau, 'L1_chrV_tau': L1_chrV_tau, 'L1_chrX_tau': L1_chrX_tau,
    'L2_chrI_tau': L2_chrI_tau, 'L2_chrII_tau': L2_chrII_tau, 'L2_chrIII_tau': L2_chrIII_tau, 
    'L2_chrIV_tau': L2_chrIV_tau, 'L2_chrV_tau': L2_chrV_tau, 'L2_chrX_tau': L2_chrX_tau,
    'L3_chrI_tau': L3_chrI_tau, 'L3_chrII_tau': L3_chrII_tau, 'L3_chrIII_tau': L3_chrIII_tau, 
    'L3_chrIV_tau': L3_chrIV_tau, 'L3_chrV_tau': L3_chrV_tau, 'L3_chrX_tau': L3_chrX_tau,
    'D_chrI_tau': D_chrI_tau, 'D_chrII_tau': D_chrII_tau, 'D_chrIII_tau': D_chrIII_tau, 
    'D_chrIV_tau': D_chrIV_tau, 'D_chrV_tau': D_chrV_tau, 'D_chrX_tau': D_chrX_tau,
    'L4_chrI_tau': L4_chrI_tau, 'L4_chrII_tau': L4_chrII_tau, 'L4_chrIII_tau': L4_chrIII_tau, 
    'L4_chrIV_tau': L4_chrIV_tau, 'L4_chrV_tau': L4_chrV_tau, 'L4_chrX_tau': L4_chrX_tau,
    'A_chrI_tau': A_chrI_tau, 'A_chrII_tau': A_chrII_tau, 'A_chrIII_tau': A_chrIII_tau, 
    'A_chrIV_tau': A_chrIV_tau, 'A_chrV_tau': A_chrV_tau, 'A_chrX_tau': A_chrX_tau
}

total_df = pd.concat(dataframes, names=['variable']).reset_index(level=0).rename(columns={'variable': 'source_variable'})
total_df.to_csv('../data/stage_biased_chromosomal_coords.csv', index=False)

#make plot
# Define colors
cmap = plt.cm.get_cmap('viridis')
num_bins = 21
colors = [cmap(i / num_bins) for i in range(num_bins)]
colors = [mcolors.to_hex(color) for color in colors]
colors = [colors[0], colors[4], colors[5], colors[8], colors[12], colors[16], colors[20]]
grey_color = 'tab:grey'

# Set seaborn style and color palette
sns.set_style('white')
sns.set_palette(colors)
fig, [ax1, ax2, ax3, ax4, ax5, ax6] = plt.subplots(6, 1, figsize=(12, 12), sharey=True)

def plot_colored_bars(ax, window_start, cohens_d, stage_index, width=10000, label=None):
    # Determine colors for each bar
    bar_colors = [grey_color if d <= 1 else colors[stage_index % len(colors)] for d in cohens_d]
    # Plot bars with colors
    ax.bar(window_start, cohens_d, width=width, color=bar_colors, label=label)

# Egg
plot_colored_bars(ax1, E_chrI_tau['window_start'], E_chrI_tau['cohens_d'], label='Embryo', stage_index=0)
plot_colored_bars(ax2, E_chrII_tau['window_start'], E_chrII_tau['cohens_d'], label='Embryo', stage_index=0)
plot_colored_bars(ax3, E_chrIII_tau['window_start'], E_chrIII_tau['cohens_d'], label='Embryo', stage_index=0)
plot_colored_bars(ax4, E_chrIV_tau['window_start'], E_chrIV_tau['cohens_d'], label='Embryo', stage_index=0)
plot_colored_bars(ax5, E_chrV_tau['window_start'], E_chrV_tau['cohens_d'], label='Embryo', stage_index=0)
plot_colored_bars(ax6, E_chrX_tau['window_start'], E_chrX_tau['cohens_d'], label='Embryo', stage_index=0)

# L1
plot_colored_bars(ax1, L1_chrI_tau['window_start'], L1_chrI_tau['cohens_d'], label='L1', stage_index=1)
plot_colored_bars(ax2, L1_chrII_tau['window_start'], L1_chrII_tau['cohens_d'], label='L1', stage_index=1)
plot_colored_bars(ax3, L1_chrIII_tau['window_start'], L1_chrIII_tau['cohens_d'], label='L1', stage_index=1)
plot_colored_bars(ax4, L1_chrIV_tau['window_start'], L1_chrIV_tau['cohens_d'], label='L1', stage_index=1)
plot_colored_bars(ax5, L1_chrV_tau['window_start'], L1_chrV_tau['cohens_d'], label='L1', stage_index=1)
plot_colored_bars(ax6, L1_chrX_tau['window_start'], L1_chrX_tau['cohens_d'], label='L1', stage_index=1)

# L2
plot_colored_bars(ax1, L2_chrI_tau['window_start'], L2_chrI_tau['cohens_d'], label='L2', stage_index=2)
plot_colored_bars(ax2, L2_chrII_tau['window_start'], L2_chrII_tau['cohens_d'], label='L2', stage_index=2)
plot_colored_bars(ax3, L2_chrIII_tau['window_start'], L2_chrIII_tau['cohens_d'], label='L2', stage_index=2)
plot_colored_bars(ax4, L2_chrIV_tau['window_start'], L2_chrIV_tau['cohens_d'], label='L2', stage_index=2)
plot_colored_bars(ax5, L2_chrV_tau['window_start'], L2_chrV_tau['cohens_d'], label='L2', stage_index=2)
plot_colored_bars(ax6, L2_chrX_tau['window_start'], L2_chrX_tau['cohens_d'], label='L2', stage_index=2)

# L3
plot_colored_bars(ax1, L3_chrI_tau['window_start'], L3_chrI_tau['cohens_d'], label='L3', stage_index=3)
plot_colored_bars(ax2, L3_chrII_tau['window_start'], L3_chrII_tau['cohens_d'], label='L3', stage_index=3)
plot_colored_bars(ax3, L3_chrIII_tau['window_start'], L3_chrIII_tau['cohens_d'], label='L3', stage_index=3)
plot_colored_bars(ax4, L3_chrIV_tau['window_start'], L3_chrIV_tau['cohens_d'], label='L3', stage_index=3)
plot_colored_bars(ax5, L3_chrV_tau['window_start'], L3_chrV_tau['cohens_d'], label='L3', stage_index=3)
plot_colored_bars(ax6, L3_chrX_tau['window_start'], L3_chrX_tau['cohens_d'], label='L3', stage_index=3)

# Dauer
plot_colored_bars(ax1, D_chrI_tau['window_start'], D_chrI_tau['cohens_d'], label='Dauer', stage_index=4)
plot_colored_bars(ax2, D_chrII_tau['window_start'], D_chrII_tau['cohens_d'], label='Dauer', stage_index=4)
plot_colored_bars(ax3, D_chrIII_tau['window_start'], D_chrIII_tau['cohens_d'], label='Dauer', stage_index=4)
plot_colored_bars(ax4, D_chrIV_tau['window_start'], D_chrIV_tau['cohens_d'], label='Dauer', stage_index=4)
plot_colored_bars(ax5, D_chrV_tau['window_start'], D_chrV_tau['cohens_d'], label='Dauer', stage_index=4)
plot_colored_bars(ax6, D_chrX_tau['window_start'], D_chrX_tau['cohens_d'], label='Dauer', stage_index=4)

# L4
plot_colored_bars(ax1, L4_chrI_tau['window_start'], L4_chrI_tau['cohens_d'], label='L4', stage_index=5)
plot_colored_bars(ax2, L4_chrII_tau['window_start'], L4_chrII_tau['cohens_d'], label='L4', stage_index=5)
plot_colored_bars(ax3, L4_chrIII_tau['window_start'], L4_chrIII_tau['cohens_d'], label='L4', stage_index=5)
plot_colored_bars(ax4, L4_chrIV_tau['window_start'], L4_chrIV_tau['cohens_d'], label='L4', stage_index=5)
plot_colored_bars(ax5, L4_chrV_tau['window_start'], L4_chrV_tau['cohens_d'], label='L4', stage_index=5)
plot_colored_bars(ax6, L4_chrX_tau['window_start'], L4_chrX_tau['cohens_d'], label='L4', stage_index=5)

# Adult
plot_colored_bars(ax1, A_chrI_tau['window_start'], A_chrI_tau['cohens_d'], label='Adult', stage_index=6)
plot_colored_bars(ax2, A_chrII_tau['window_start'], A_chrII_tau['cohens_d'], label='Adult', stage_index=6)
plot_colored_bars(ax3, A_chrIII_tau['window_start'], A_chrIII_tau['cohens_d'], label='Adult', stage_index=6)
plot_colored_bars(ax4, A_chrIV_tau['window_start'], A_chrIV_tau['cohens_d'], label='Adult', stage_index=6)
plot_colored_bars(ax5, A_chrV_tau['window_start'], A_chrV_tau['cohens_d'], label='Adult', stage_index=6)
plot_colored_bars(ax6, A_chrX_tau['window_start'], A_chrX_tau['cohens_d'], label='Adult', stage_index=6)

# Add legend
legend_labels = ['Embryo', 'L1', 'L2', 'L3', 'Dauer', 'L4', 'Adult']
legend_colors = colors
handles = [plt.Line2D([0], [0], marker='s', color='w', label=label, markersize=10, markerfacecolor=color) 
           for label, color in zip(legend_labels, legend_colors)]

ax1.legend(handles=handles, title='Stage Bias', loc='upper right', bbox_to_anchor=(1.1075, 1.05))

# Add horizontal lines
for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.axhline(y=1, color='tab:gray', linestyle='--')
    ax.axhline(y=-1, color='tab:gray', linestyle='--')

# Adjust plot text
ax6.set_xlabel('Genomic position (100 Kb windows)', fontsize=15)

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.set_ylabel("Cohen's d", fontsize=15)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

ax1.text(-700000,2, 'Chromosome I', fontsize=15)
ax2.text(-700000,2, 'Chromosome II', fontsize=15)
ax3.text(-650000,2, 'Chromosome III', fontsize=15)
ax4.text(-800000,2, 'Chromosome IV', fontsize=15)
ax5.text(-950000,2, 'Chromosome V', fontsize=15)
ax6.text(-800000,2, 'Chromosome X', fontsize=15)
    
plt.tight_layout()
plt.savefig('../figures/stage_biased_regions.pdf')

#get gene cluster size distributions
#manage data a bit
E_chrI_tau['chromosome'] = 'I'
E_chrII_tau['chromosome'] = 'II'
E_chrIII_tau['chromosome'] = 'III'
E_chrIV_tau['chromosome'] = 'IV'
E_chrV_tau['chromosome'] = 'V'
E_chrX_tau['chromosome'] = 'X'
E_biased = pd.concat([E_chrI_tau, E_chrII_tau, E_chrIII_tau, E_chrIV_tau, E_chrV_tau, E_chrX_tau], ignore_index=True)

L1_chrI_tau['chromosome'] = 'I'
L1_chrII_tau['chromosome'] = 'II'
L1_chrIII_tau['chromosome'] = 'III'
L1_chrIV_tau['chromosome'] = 'IV'
L1_chrV_tau['chromosome'] = 'V'
L1_chrX_tau['chromosome'] = 'X'
L1_biased = pd.concat([L1_chrI_tau, L1_chrII_tau, L1_chrIII_tau, L1_chrIV_tau, L1_chrV_tau, L1_chrX_tau], ignore_index=True)

L2_chrI_tau['chromosome'] = 'I'
L2_chrII_tau['chromosome'] = 'II'
L2_chrIII_tau['chromosome'] = 'III'
L2_chrIV_tau['chromosome'] = 'IV'
L2_chrV_tau['chromosome'] = 'V'
L2_chrX_tau['chromosome'] = 'X'
L2_biased = pd.concat([L2_chrI_tau, L2_chrII_tau, L2_chrIII_tau, L2_chrIV_tau, L2_chrV_tau, L2_chrX_tau], ignore_index=True)

L3_chrI_tau['chromosome'] = 'I'
L3_chrII_tau['chromosome'] = 'II'
L3_chrIII_tau['chromosome'] = 'III'
L3_chrIV_tau['chromosome'] = 'IV'
L3_chrV_tau['chromosome'] = 'V'
L3_chrX_tau['chromosome'] = 'X'
L3_biased = pd.concat([L3_chrI_tau, L3_chrII_tau, L3_chrIII_tau, L3_chrIV_tau, L3_chrV_tau, L3_chrX_tau], ignore_index=True)

D_chrI_tau['chromosome'] = 'I'
D_chrII_tau['chromosome'] = 'II'
D_chrIII_tau['chromosome'] = 'III'
D_chrIV_tau['chromosome'] = 'IV'
D_chrV_tau['chromosome'] = 'V'
D_chrX_tau['chromosome'] = 'X'
D_biased = pd.concat([D_chrI_tau, D_chrII_tau, D_chrIII_tau, D_chrIV_tau, D_chrV_tau, D_chrX_tau], ignore_index=True)

L4_chrI_tau['chromosome'] = 'I'
L4_chrII_tau['chromosome'] = 'II'
L4_chrIII_tau['chromosome'] = 'III'
L4_chrIV_tau['chromosome'] = 'IV'
L4_chrV_tau['chromosome'] = 'V'
L4_chrX_tau['chromosome'] = 'X'
L4_biased = pd.concat([L4_chrI_tau, L4_chrII_tau, L4_chrIII_tau, L4_chrIV_tau, L4_chrV_tau, L4_chrX_tau], ignore_index=True)

A_chrI_tau['chromosome'] = 'I'
A_chrII_tau['chromosome'] = 'II'
A_chrIII_tau['chromosome'] = 'III'
A_chrIV_tau['chromosome'] = 'IV'
A_chrV_tau['chromosome'] = 'V'
A_chrX_tau['chromosome'] = 'X'
A_biased = pd.concat([A_chrI_tau, A_chrII_tau, A_chrIII_tau, A_chrIV_tau, A_chrV_tau, A_chrX_tau], ignore_index=True)

#get windows that are significantly biased
E_biased = E_biased.loc[E_biased['cohens_d'] >=1]
L1_biased = L1_biased.loc[L1_biased['cohens_d'] >=1]
L2_biased = L2_biased.loc[L2_biased['cohens_d'] >=1]
L3_biased = L3_biased.loc[L3_biased['cohens_d'] >=1]
D_biased = D_biased.loc[D_biased['cohens_d'] >=1]
L4_biased = L4_biased.loc[L4_biased['cohens_d'] >=1]
A_biased = A_biased.loc[A_biased['cohens_d'] >=1]

#write to dataframe
dataframes = {
    'E_biased': E_biased, 'L1_biased': L1_biased, 'L2_biased' : L2_biased,
    'L3_biased' : L3_biased, 'D_biased' : D_biased, 'L4_biased' : L4_biased,
    'A_biased' : A_biased
}

total_biased_df = pd.concat(dataframes, names=['variable']).reset_index(level=0).rename(columns={'variable': 'source_variable'})
total_biased_df.to_csv('../data/significantly_stage_biased_chromosomal_coords.csv', index=False)