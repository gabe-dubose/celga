#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

#load data
data = pd.read_csv('../data/stage_biased_regions_duplication.csv')
data = data.dropna(subset=['n_genes', 'non_homolog_ratio'])

sns.set_style('whitegrid')
fig, [ax1,ax2] = plt.subplots(1, 2, figsize=(10, 5))

#plot distribution
sns.ecdfplot(data=data, x='non_homolog_ratio', ax=ax1, linewidth=2, color='black')
#style
ax1.set_ylabel("Proportion", fontsize=15)
ax1.set_xlabel("Non-homologous ratio\n(n homologous groups/gene density)", fontsize=15)
ax1.tick_params(axis='x', labelsize=15)
ax1.tick_params(axis='y', labelsize=15)

# Fit lowess regression model
x = data['n_genes']
y = data['non_homolog_ratio']
lowess_result = lowess(y, x, frac=0.75)
lowess_x = lowess_result[:, 0]
lowess_y = lowess_result[:, 1]

# Bootstrap for confidence intervals
n_boot = 1000
y_boot = np.empty((n_boot, len(lowess_x)))

for i in range(n_boot):
    sample = data.sample(frac=1, replace=True)
    y_boot[i] = lowess(sample['non_homolog_ratio'], sample['n_genes'], frac=0.75)[:, 1]

ci_lower = np.percentile(y_boot, 2.5, axis=0)
ci_upper = np.percentile(y_boot, 97.5, axis=0)

#plot correlation with gene density
#sns.regplot(data=data, y='non_homolog_ratio', x='n_genes', ax=ax2, color = 'black', scatter=False)
sns.scatterplot(data=data, y='non_homolog_ratio', x='n_genes', ax=ax2, color = 'black', alpha=1, s=75, edgecolor=None)
ax2.plot(lowess_x, lowess_y, color='black')
ax2.plot(lowess_x, ci_lower, linestyle='--', color='black')
ax2.plot(lowess_x, ci_upper, linestyle='--', color='black')

#style
ax2.set_ylabel("Non-homologous ratio", fontsize=15)
ax2.set_xlabel("Gene density", fontsize=15)
ax2.tick_params(axis='x', labelsize=15)
ax2.tick_params(axis='y', labelsize=15)

ax1.text(-0.11, 1.1, 'A', transform=ax1.transAxes, fontsize=15, fontweight='bold', va='top', ha='right')
ax2.text(-0.11, 1.1, 'B', transform=ax2.transAxes, fontsize=15, fontweight='bold', va='top', ha='right')
ax2.text(22, 0.95, r'$\rho = -0.32$', fontsize=15)

plt.tight_layout()
plt.savefig('../figures/non-homologous-ratio.pdf')