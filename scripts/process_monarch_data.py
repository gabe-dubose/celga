#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json

#load data
data = pd.read_csv('../data/dpl_tpm_counts_kallisto.csv', index_col=0)

#remove infected
data = data[~data.index.str[7].str.contains('i')]

#log transform
data = np.log1p(data)

#get average by stage
data['stage'] = data.index.str[5]
data = data.groupby('stage').mean()

#adjust gene ids
data.columns = [col.split('_')[2] for col in data.columns]

#transpose data
data = data.T

#A function to calculate specificity index (tau) for a row in a pandas data frame
def compute_tau(row):
    if max(row) != 0:
        #normalize to maximal component value
        normalized_row = [value / max(row) for value in row]
        row_to_sum = [1-value for value in normalized_row]
        row_sum = sum(row_to_sum)
        n_minus_1 = len(normalized_row) - 1
        tau = row_sum / n_minus_1
        return tau
    else:
        return 'NA'
    
ids = []
tau_values = []
max_expression_values = []
stage_max = []

for index, row in data.iterrows():
    #add gene id
    gene_id = str(index)
    ids.append(gene_id)
    
    #add tau
    values = list(row)
    tau = compute_tau(values)
    tau_values.append(tau)
    
    #add max
    max_expression = max(row)
    max_expression_values.append(max_expression)
    
    # add the stage with the max value
    max_stage = row.idxmax()
    stage_max.append(max_stage)
    
#compute tau values
tau_data = {'gene_id' : ids, 'tau' : tau_values, 'max_expression' : max_expression_values, 'stage_bias': stage_max}
tau_dataframe = pd.DataFrame(tau_data)
tau_dataframe = tau_dataframe.loc[tau_dataframe['tau'] != 'NA']

tau_dataframe['tau'] = pd.to_numeric(tau_dataframe['tau'])
tau_dataframe['max_expression'] = pd.to_numeric(tau_dataframe['max_expression'])

#write to file
tau_dataframe.to_csv('../data/dpl_stage_specificity_data.csv', index=False)

#now to calculate stats for homologous groups
with open('../data/dpl_psiblast_id-30_e-neg-10_cov-1_sequence_clusters.json') as file:
    homolog_data = json.load(file)

#initialize vectors to store ifo
group_id = []
group_size = []
group_tau_mean = []
group_tau_std = []
group_stage_biase_frequency = []

#iterate throgh each group
for group in homolog_data:
    gene_ids = homolog_data[group]
    #only consider groups with enough duplicates for comparisons
    if len(gene_ids) >= 4:
        group_data = tau_dataframe[tau_dataframe['gene_id'].isin(gene_ids)]
        #add data
        group_id.append(group)
        group_size.append(len(group_data))
        group_tau_std.append(group_data['tau'].std())
        group_tau_mean.append(group_data['tau'].mean())
        #calculate stage bias frequency
        stage_counts = group_data['stage_bias'].value_counts()
        most_frequent_stage = stage_counts.idxmax()
        most_frequent_count = stage_counts.max()
        total_entries = len(group_data)
        group_stage_biase_frequency.append(most_frequent_count / total_entries)

data = {
    'group.id': group_id,
    'group.size': group_size,
    'group.tau.mean': group_tau_mean,
    'group.tau.stdev': group_tau_std,
    'group.stage.bias.frequency': group_stage_biase_frequency
}

data = pd.DataFrame(data)

data.to_csv('../data/dpl_tau_homologous_group_data.csv', index=False)