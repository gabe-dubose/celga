#!/usr/bin/env python3

##### NOTES #####
# Raw expression data is from:
#  Boeck et. al., 2016. The time-resolved transcriptome of C. elegans. Genome Research
# This data is the Supplemental_table_S2.tsv file. 
#################

import pandas as pd
import numpy as np

#load data
counts_data = pd.read_csv('../data/Supplemental_Table_S2.tsv', sep=" ")
counts_data.rename(columns={"WormbaseName": "gene_id_wormbase"}, inplace=True)
metadata = pd.read_csv('../data/metadata.csv')

#filter for specific samples
#assemble list of samples
samples = list(metadata['sample_id'])
samples = [sample + "_dcpm" for sample in samples]
samples.append('gene_id_wormbase')
#remove samples that aren't present
not_present = ['N2_EE_50-90_dcpm', 'N2_EE_50-120_dcpm', 'N2_EE_50-150_dcpm', 'N2_EE_50-180_dcpm']
samples = [sample for sample in samples if sample not in not_present]

#filter counts data to only include relevant samples
counts_data = counts_data[samples]
counts_data.set_index('gene_id_wormbase', inplace=True)
counts_data = counts_data.transpose()
#write to file
counts_data.to_csv('../data/counts_data.csv')
#log transform
counts_data_log = np.log(counts_data + 1)
counts_data_log.to_csv('../data/counts_data_log.csv')

#get average tables
metadata = metadata.set_index('sample_id')
counts_data_log.index = counts_data_log.index.str.replace('_dcpm', '', regex=False)
merged_data_log = counts_data_log.join(metadata, how='inner')

#get average expression by life stage
average_expression_by_stage = merged_data_log.groupby('stage').mean(numeric_only=True)
average_expression_by_stage.to_csv('../data/average_log_expression_by_stage.csv')