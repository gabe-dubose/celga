#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
from scipy.stats import ks_2samp

##### NOTES #####
#The dominance index is calcualted as (Nmax/N)
#################

#load data
data = pd.read_csv('../data/cele_homologous_group_tau_var_orthogroup_only.csv')
dpl_data = pd.read_csv('../data/dplex_homologous_group_tau_var_orthogroup_only.csv')
dmel_data = pd.read_csv('../data/dmel_homologous_group_tau_var_orthogroup_only.csv')


# calculate standard deviations
data['group.tau.sd'] = np.sqrt(data['group.tau.variance....group.tau.variance'])
dmel_data['group.tau.sd'] = np.sqrt(dmel_data['group.tau.variance....group.tau.variance'])
dpl_data['group.tau.sd'] = np.sqrt(dpl_data['group.tau.variance....group.tau.variance'])

#compare cele and dpl stage-bias conservation
sns.set_style('whitegrid')
fig, [ax1,ax2] = plt.subplots(1, 2, figsize=(8, 4))

#cele
sns.ecdfplot(data=data, x='group.tau.sd', ax=ax1, 
             linewidth=2, color='black')
#dpl
sns.ecdfplot(data=dpl_data, x='group.tau.sd', ax=ax1,
             linestyle='--', linewidth=2, color='tab:gray')
#dmel
sns.ecdfplot(data=dmel_data, x='group.tau.sd', ax=ax1,
             linestyle=':', linewidth=2, color='tab:gray')

#add legend
#define legend
legend_lines = [mlines.Line2D([0], [0], color='black', label=r'${C. elegans}$', linewidth=2),
                mlines.Line2D([0], [0], color='tab:gray', label=r'${D. plexippus}$', linewidth=2, linestyle='--'),
                mlines.Line2D([0], [0], color='tab:gray', label=r'${D. melanogaster}$', linewidth=2, linestyle=':')]

ax1.legend(handles=legend_lines, loc='lower right', fontsize=10)
ax1.set_ylabel('Proportion', fontsize=15)
ax1.set_xlabel(r'$\tau$ st.dev.', fontsize=15)
ax1.tick_params(axis='x', labelsize=15)
ax1.tick_params(axis='y', labelsize=15)

#cele
sns.ecdfplot(data=data, x='group.stage.bias.frequency....group.stage.bias.frequency', ax=ax2, 
             linewidth=2, color='black')
#dpl
sns.ecdfplot(data=dpl_data, x='group.stage.bias.frequency....group.stage.bias.frequency', ax=ax2,
             linestyle='--', linewidth=2, color='tab:gray')
#dmel
sns.ecdfplot(data=dmel_data, x='group.stage.bias.frequency....group.stage.bias.frequency', ax=ax2,
             linestyle=':', linewidth=2, color='tab:gray')

#add legend
ax2.legend(handles=legend_lines, loc='lower right', fontsize=10)
ax2.set_ylabel('Proportion', fontsize=15)
ax2.set_xlabel(r'Stage-bias uniformity ($N_{max}/N$)', fontsize=15)
ax2.tick_params(axis='x', labelsize=15)
ax2.tick_params(axis='y', labelsize=15)

ax1.text(-0.15, 1.1, 'A', transform=ax1.transAxes, fontsize=15, fontweight='bold', va='top', ha='right')
ax2.text(-0.15, 1.1, 'B', transform=ax2.transAxes, fontsize=15, fontweight='bold', va='top', ha='right')

plt.tight_layout()

plt.savefig('../figures/cele_vs_dplex_stage_specificity_gene_duplication_orthogroups_only.pdf')
plt.savefig('../figures/cele_vs_dplex_stage_specificity_gene_duplication_orthogroups_only.png', dpi=600)

##### STATISTICS #####
#stdev
# cele vs dplex
ks_statistic, p_value = ks_2samp(list(data['group.tau.sd'].dropna()), list(dpl_data['group.tau.sd'].dropna()))
print(f"Standard deviation: D = {ks_statistic}, p = {p_value}")
#Standard deviation: D = 0.3905308464849354, p = 3.1016635845187326e-06

# cele vs dmel
ks_statistic, p_value = ks_2samp(list(data['group.tau.sd'].dropna()), list(dmel_data['group.tau.sd'].dropna()))
print(f"Standard deviation: D = {ks_statistic}, p = {p_value}")
#Standard deviation: D = 0.36949891067538126, p = 0.0050830744125562035

# stage bias uniformity
# cele vs dplex
ks_statistic, p_value = ks_2samp(list(data['group.stage.bias.frequency....group.stage.bias.frequency']), list(dpl_data['group.stage.bias.frequency....group.stage.bias.frequency']))
print(f"Stage-bias uniformity: D = {ks_statistic}, p = {p_value}")
#Stage-bias uniformity: D = 0.6283887931816121, p = 5.1715129109898295e-24

# cele vs dmel
ks_statistic, p_value = ks_2samp(list(data['group.stage.bias.frequency....group.stage.bias.frequency']), list(dmel_data['group.stage.bias.frequency....group.stage.bias.frequency']))
print(f"Stage-bias uniformity: D = {ks_statistic}, p = {p_value}")
#Stage-bias uniformity: D = 0.5351832182687202, p = 2.8511800206908317e-07
######################