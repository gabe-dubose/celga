#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.lines as mlines

#load data
data = pd.read_csv('../data/tau_cluster_data.csv')

# Set the style
sns.set_style('whitegrid')

# Create the plot
fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))

# Define colors
cmap = plt.cm.get_cmap('viridis')
num_bins = 21
colors = [cmap(i / num_bins) for i in range(num_bins)]
colors = [mcolors.to_hex(color) for color in colors]
colors = [colors[0], colors[4], colors[12], colors[16], colors[20]]

# Define linestyles
linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]

# Plot each line with different colors and linestyles using line_kws
sns.regplot(data=data.loc[data['stage_bias'] == 'Egg'], x='tau', y='embryo_clustered', logistic=True, ax=ax1,
            scatter=False, color=colors[0], line_kws={'linestyle': linestyles[0]})

sns.regplot(data=data.loc[data['stage_bias'] == 'L1'], x='tau', y='L1_clustered', logistic=True, ax=ax1,
            scatter=False, color=colors[1], line_kws={'linestyle': linestyles[1]})

#note that relationship between expression specificity and probability of being cluster associated is not significant

sns.regplot(data=data.loc[data['stage_bias'] == 'Dauer'], x='tau', y='Dauer_clustered', logistic=True, ax=ax1,
            scatter=False, color=colors[2], line_kws={'linestyle': linestyles[2]})

sns.regplot(data=data.loc[data['stage_bias'] == 'L4'], x='tau', y='L4_clustered', logistic=True, ax=ax1,
            scatter=False, color=colors[3], line_kws={'linestyle': linestyles[3]})

sns.regplot(data=data.loc[data['stage_bias'] == 'Adult'], x='tau', y='Adult_clustered', logistic=True, ax=ax1,
            scatter=False, color=colors[4], line_kws={'linestyle': linestyles[4]})

ax1.set_ylim(0, 1)
ax1.set_ylabel('Proportion of genes in stage-biased region', fontsize=15)
ax1.set_xlabel('Stage-specificity (τ)', fontsize=15)
ax1.tick_params(axis='x', labelsize=15)
ax1.tick_params(axis='y', labelsize=15)

# Create legend handles with both color and linestyle
legend_labels = ['Embryo', 'L1', 'Dauer', 'L4', 'Adult']
handles = [plt.Line2D([0], [0], color=color, linestyle=linestyle, label=label) 
           for label, color, linestyle in zip(legend_labels, colors, linestyles)]

# Add legend
ax1.legend(handles=handles, title='Stage Bias', loc='upper left')

plt.tight_layout()

plt.savefig('../figures/clustered_by_tau_logistic.pdf')

'''
Notes on statistics

Embryo model:
Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept) -4.99178    0.09606  -51.97   <2e-16 ***
tau          5.22699    0.16005   32.66   <2e-16 ***

    Null deviance: 6237.6  on 13371  degrees of freedom
Residual deviance: 5032.1  on 13370  degrees of freedom
AIC: 5036.1

#L1 model
Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  -10.693      1.271  -8.414  < 2e-16 ***
tau           11.360      1.616   7.030 2.07e-12 ***

    Null deviance: 329.41  on 800  degrees of freedom
Residual deviance: 232.54  on 799  degrees of freedom
AIC: 236.54

#Dauer model
Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  -7.8019     0.6000  -13.00   <2e-16 ***
tau           6.5640     0.7899    8.31   <2e-16 ***

    Null deviance: 799.92  on 2973  degrees of freedom
Residual deviance: 694.40  on 2972  degrees of freedom
AIC: 698.4

#L4 model
Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  -15.848      3.202   -4.95 7.43e-07 ***
tau           14.386      3.727    3.86 0.000114 ***

    Null deviance: 154.06  on 1269  degrees of freedom
Residual deviance: 126.79  on 1268  degrees of freedom
AIC: 130.79

#Adult model
Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)   -9.808      2.048  -4.790 1.67e-06 ***
tau           10.023      2.482   4.038 5.39e-05 ***

    Null deviance: 146.758  on 341  degrees of freedom
Residual deviance:  99.503  on 340  degrees of freedom
AIC: 103.5

'''