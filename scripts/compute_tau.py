#!/usr/bin/env python3

import pandas as pd

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
    
#load data
expression_data = pd.read_csv('../data/average_log_expression_by_stage.csv', index_col=0).rename_axis("genes")
expression_data = expression_data.transpose()

ids = []
tau_values = []
max_expression_values = []
stage_max = []

for index, row in expression_data.iterrows():
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
tau_dataframe.to_csv('../data/stage_specificity_data.csv', index=False)