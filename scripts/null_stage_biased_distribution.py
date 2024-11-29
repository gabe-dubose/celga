#!/usr/bin/env python3

##### NOTES #####
# The goal of this simulation is to estimate the frequency and size of 
# stage-biased gene clusters that can be expected if the arrangement
# of stage-biased patterns of expression was random. To accomplish this, 
# the null model algorithm will be a Monte-Carlo simulation that randomly distribute expression 
# patterns (which stage the expression is biased towards and its specificity to said stage).
# It will then calculate the distribution of stage biased regions in the same way that 
# they were quantified. This process will be repeated to generate a null distribution. 
#################

import pandas as pd
import numpy as np

#load stage-specificity data
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

#A function to randomly shuffel tau dataframe
def random_shuffle(tau_dataframe):
    random_df = tau_dataframe
    tau_shuffled = random_df['gene_id'].sample(frac=1).reset_index(drop=True)
    random_df['gene_id'] = tau_shuffled
    return random_df

#A function to run Monte-Carlo simulation
def null_simulation(tau_data, iterations):
    
    #initialize output df
    columns = ["source_variable", "window_start", "window_stop", 
    "average_tau", "max_expression", "density", 
    "cohens_d", "chromosome", "iteration"]

    full_null_df = pd.DataFrame(columns=columns)

    #start simulations
    for i in range(iterations):
        random_tau_distribution = random_shuffle(tau_data)
    
        #separate by chromosome
        #chromosome I = NC_003279.8
        chromosome_I = genome_data.loc[genome_data['chromosome'] == 'NC_003279.8']
        chromosome_I = pd.merge(chromosome_I, random_tau_distribution, left_on='gene', right_on='gene_id', how='inner')
        chromosome_I = chromosome_I.sort_values(by='start', ascending=True)
        chromosome_I['gene_position'] = range(1, len(chromosome_I) + 1)

        #chromosome II = NC_003280.10
        chromosome_II = genome_data.loc[genome_data['chromosome'] == 'NC_003280.10']
        chromosome_II = pd.merge(chromosome_II, random_tau_distribution, left_on='gene', right_on='gene_id', how='inner')
        chromosome_II = chromosome_II.sort_values(by='start', ascending=True)
        chromosome_II['gene_position'] = range(1, len(chromosome_II) + 1)

        #chromosome III = NC_003281.10
        chromosome_III = genome_data.loc[genome_data['chromosome'] == 'NC_003281.10']
        chromosome_III = pd.merge(chromosome_III, random_tau_distribution, left_on='gene', right_on='gene_id', how='inner')
        chromosome_III = chromosome_III.sort_values(by='start', ascending=True)
        chromosome_III['gene_position'] = range(1, len(chromosome_III) + 1)

        #chromosome IV = NC_003282.8
        chromosome_IV = genome_data.loc[genome_data['chromosome'] == 'NC_003282.8']
        chromosome_IV = pd.merge(chromosome_IV, random_tau_distribution, left_on='gene', right_on='gene_id', how='inner')
        chromosome_IV = chromosome_IV.sort_values(by='start', ascending=True)
        chromosome_IV['gene_position'] = range(1, len(chromosome_IV) + 1)

        #chromosome V = NC_003283.11
        chromosome_V = genome_data.loc[genome_data['chromosome'] == 'NC_003283.11']
        chromosome_V = pd.merge(chromosome_V, random_tau_distribution, left_on='gene', right_on='gene_id', how='inner')
        chromosome_V = chromosome_V.sort_values(by='start', ascending=True)
        chromosome_V['gene_position'] = range(1, len(chromosome_V) + 1)

        #chromosome X = NC_003284.9
        chromosome_X = genome_data.loc[genome_data['chromosome'] == 'NC_003284.9']
        chromosome_X = pd.merge(chromosome_X, random_tau_distribution, left_on='gene', right_on='gene_id', how='inner')
        chromosome_X = chromosome_X.sort_values(by='start', ascending=True)
        chromosome_X['gene_position'] = range(1, len(chromosome_X) + 1)

        #get sliding window data for egg
        E_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'Egg'], 100000)
        E_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'Egg'], 100000)
        E_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'Egg'], 100000)
        E_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'Egg'], 100000)
        E_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'Egg'], 100000)
        E_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'Egg'], 100000)

        #get sliding window data for L1
        L1_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'L1'], 100000)
        L1_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'L1'], 100000)
        L1_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'L1'], 100000)
        L1_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'L1'], 100000)
        L1_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'L1'], 100000)
        L1_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'L1'], 100000)

        #get sliding window data for L2
        L2_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'L2'], 100000)
        L2_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'L2'], 100000)
        L2_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'L2'], 100000)
        L2_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'L2'], 100000)
        L2_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'L2'], 100000)
        L2_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'L2'], 100000)

        #get sliding window data for L3
        L3_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'L3'], 100000)
        L3_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'L3'], 100000)
        L3_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'L3'], 100000)
        L3_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'L3'], 100000)
        L3_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'L3'], 100000)
        L3_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'L3'], 100000)

        #get sliding window data for dauer
        D_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'Dauer'], 100000)
        D_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'Dauer'], 100000)
        D_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'Dauer'], 100000)
        D_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'Dauer'], 100000)
        D_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'Dauer'], 100000)
        D_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'Dauer'], 100000)

        #get sliding window data for L4
        L4_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'L4'], 100000)
        L4_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'L4'], 100000)
        L4_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'L4'], 100000)
        L4_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'L4'], 100000)
        L4_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'L4'], 100000)
        L4_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'L4'], 100000)

        #get sliding window data for Adults
        A_chrI_tau = sliding_window_expression_specificity(chromosome_I.loc[chromosome_I['stage_bias'] == 'Adult'], 100000)
        A_chrII_tau = sliding_window_expression_specificity(chromosome_II.loc[chromosome_II['stage_bias'] == 'Adult'], 100000)
        A_chrIII_tau = sliding_window_expression_specificity(chromosome_III.loc[chromosome_III['stage_bias'] == 'Adult'], 100000)
        A_chrIV_tau = sliding_window_expression_specificity(chromosome_IV.loc[chromosome_IV['stage_bias'] == 'Adult'], 100000)
        A_chrV_tau = sliding_window_expression_specificity(chromosome_V.loc[chromosome_V['stage_bias'] == 'Adult'], 100000)
        A_chrX_tau = sliding_window_expression_specificity(chromosome_X.loc[chromosome_X['stage_bias'] == 'Adult'], 100000)

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

        #combine in
        dataframes = {
            'E_biased': E_biased, 'L1_biased': L1_biased, 'L2_biased' : L2_biased,
            'L3_biased' : L3_biased, 'D_biased' : D_biased, 'L4_biased' : L4_biased,
            'A_biased' : A_biased
        }

        null_df = pd.concat(dataframes, names=['variable']).reset_index(level=0).rename(columns={'variable': 'source_variable'})
        null_df['iteration'] = i

        #add null iteration to dataframe
        full_null_df = pd.concat([full_null_df, null_df], ignore_index=True)

        #update
        print(f"Iteration {i+1} of  {iterations} complete.")

    return full_null_df

#run simulations
null_distributions = null_simulation(tau_dataframe, 1000)

#write to file
null_distributions.to_csv('../data/null_gene_cluster_distributions.csv', index=False)