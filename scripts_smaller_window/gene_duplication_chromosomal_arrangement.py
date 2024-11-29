#!/usr/bin/env python3

import pandas as pd
import json
import uuid

#load stage biased regions data
stage_biased_regions = pd.read_csv('../data/significantly_stage_biased_chromosomal_coords.csv')

#load tau cluster data
tau_cluster_data = pd.read_csv('../data/tau_cluster_data.csv')

#load homologous group data
with open('../data/duplication/homologous_groups_strict_cutoffs.json', 'r') as file:
    homolog_data = json.load(file)

#reconfigure homolgous group data to be searchable by gene ids
#initialize new dictionary
homolog_data_by_gene = {}
for group in homolog_data:
    gene_ids = homolog_data[group]
    #go through each gene
    for gene_id in gene_ids:
        #reformat a little
        gene_id = gene_id.split('_')[1]
        #add to group
        homolog_data_by_gene[gene_id] = group

#initialize pandas dataframe to hold data
duplication_data = pd.DataFrame(columns=['chromosome', 'window_start', 'window_stop', 'stage_bias', 'non_homolog_ratio', 'n_genes'])

#define stages
stages = {'E_biased' : 'Egg', 'L1_biased' : 'L1', 'L2_biased' : 'L2',
          'L3_biased' : 'L3', 'D_biased' : 'Dauer', 'L4_biased' : 'L4',
          'A_biased' : 'Adult'}

#iterate through stage-biased regions
for index, row in stage_biased_regions.iterrows():
    #get relevant info
    chromosome = row["chromosome"]
    region_start = row["window_start"]
    region_stop = row["window_stop"]
    bias = row['source_variable']
    n_genes = row['density']

    print(f"{bias}\t{chromosome}\t{region_start}\t{region_stop}\t{n_genes}")

    #get the genes in this region from tau_cluster_data
    chromosome_data = tau_cluster_data.loc[tau_cluster_data['chromosome'] == chromosome]
    genes = chromosome_data[(chromosome_data['start'] >= region_start) & (chromosome_data['stop'] <= region_stop)]
    #iterate through each stage
    stage = stages[bias]
    stage_biased_data = genes.loc[genes['stage_bias'] == stage]
    #now get groups of stage biased genes
    gene_ids = list(stage_biased_data['gene_id'])
    #get what groups each gene is apart of
    groups = [homolog_data_by_gene[gene] if gene in homolog_data_by_gene else str(uuid.uuid4()) for gene in gene_ids]
    groups = list(set(groups))
    #define n_groups
    n_groups = len(groups)

    #if region is stage biased
    #calculate homolog / n_gene ratio
    non_homolog_ratio = n_groups / n_genes

    results = [chromosome, region_start, region_stop, stage, non_homolog_ratio, n_genes]
    duplication_data = pd.concat([duplication_data, pd.DataFrame([results], columns=duplication_data.columns)], ignore_index=True)

duplication_data.to_csv('../data/stage_biased_regions_duplication.csv', index=False)

#1. iterate through each stage-biased region for each stage
#2. get all of the genes in that region from tau cluster data
#3. get which (if any) homologous group each gene is in
#4. Then, for each region, calculate non-homolog ratio:
#       n homologous groups / n genes

#   for example:
#       1 homologous group / 10 genes = 1/10
#       1 homologous group / 6 genes = 1/6
#       3 homologous groups / 6 genes = 3/6
#       6 homologous groups / 6 genes = 1
#   Therefore, higher values (approaching 1) indicate the group is comprised largley of non-homologous genes
#       and that their clustering cannot entirely be explained by gene duplication
#   Conversely, lower numbers indicate that the cluster is comprised of more homologous genes,
#       which could contribute to their expression clustering. 
#   Even lower numbers indicate even larger regions that are comprised of a single homologous group

# Note that this approach does not count anything if a gene did not have any homologs. 
# For example, if I say added a "no_homologs" tag, and all genes in a group were from unique groups,
# this would register as a value of 1/6, since the algorithm counts only unique groups.
# Therefore, if a group consistents totally of genes with no detected homologs, a unique id is added.
# Therefore, each gene would get a unique id, and it would still tally up to 6/6 = 1