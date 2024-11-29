#!/usr/bin/env python3

import pandas as pd

#load tau data
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

chromosome_ids = {'NC_003279.8': 'I', 'NC_003280.10': 'II', 'NC_003281.10': 'III',
                'NC_003282.8': 'IV', 'NC_003283.11': 'V', 'NC_003284.9': 'X'}

genome_data['chromosome'] = genome_data['chromosome'].map(chromosome_ids)

# Function to find the cluster value
def find_clustered_value(row, total_biased):
    # Filter total_biased for the same chromosome
    chromosome_matches = total_biased[total_biased['chromosome'] == row['chromosome']]
    
    # Find the window where the start value is between window_start and window_stop
    window = chromosome_matches[
        (chromosome_matches['window_start'] <= row['start']) & 
        (chromosome_matches['window_stop'] >= row['start'])
    ]
    
    # If a matching window is found
    if not window.empty:
        cohens_d = window.iloc[0]['cohens_d']
        # Return 1 if cohens_d is >= 1, else 0
        return 1 if cohens_d >= 1 else 0
    
    # Return 0 if no matching window is found
    return 0

#load chromosomal coords
biased_regions = pd.read_csv('../data/significantly_stage_biased_chromosomal_coords.csv')
#process into different dataframes by stage
E_biased = biased_regions.loc[biased_regions['source_variable'] == 'E_biased']
L1_biased = biased_regions.loc[biased_regions['source_variable'] == 'L1_biased']
L2_biased = biased_regions.loc[biased_regions['source_variable'] == 'L2_biased']
L3_biased = biased_regions.loc[biased_regions['source_variable'] == 'L3_biased']
D_biased = biased_regions.loc[biased_regions['source_variable'] == 'D_biased']
L4_biased = biased_regions.loc[biased_regions['source_variable'] == 'L4_biased']
A_biased = biased_regions.loc[biased_regions['source_variable'] == 'A_biased']

#add embryo biased clusters
genome_data['embryo_clustered'] = genome_data.apply(lambda row: find_clustered_value(row, E_biased), axis=1)

#add L1 biased clusters
genome_data['L1_clustered'] = genome_data.apply(lambda row: find_clustered_value(row, L1_biased), axis=1)

#add L2 biased clusters
genome_data['L2_clustered'] = genome_data.apply(lambda row: find_clustered_value(row, L2_biased), axis=1)

#Note: there's only 1 L2 biased cluster
#Note: there are no L3 biased clusters

#add Dauer biased clusters
genome_data['Dauer_clustered'] = genome_data.apply(lambda row: find_clustered_value(row, D_biased), axis=1)

#add L4 biased clusters
genome_data['L4_clustered'] = genome_data.apply(lambda row: find_clustered_value(row, L4_biased), axis=1)

#add Adult biased clusters
genome_data['Adult_clustered'] = genome_data.apply(lambda row: find_clustered_value(row, A_biased), axis=1)

#add tau data
genome_data = pd.merge(genome_data, tau_dataframe, how='left', left_on='gene', right_on='gene_id')

genome_data.to_csv('../data/tau_cluster_data.csv', index=False)