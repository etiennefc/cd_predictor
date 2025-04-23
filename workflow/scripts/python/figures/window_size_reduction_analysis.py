#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np 

# Load dfs
half_fixed_len = int(int(snakemake.params.fixed_length) / 2)
positives_info = pd.read_csv(snakemake.input.positives_info, sep='\t')
positives_info['half_len'] = (positives_info['end'] - positives_info['start'] + 1) / 2
snoBIRD_N_preds = pd.read_csv(snakemake.input.snoBIRD_N_preds, sep='\t', names=['chr', 'start', 'end', 'gene_id', 'score', 'strand', 'len'])
snoBIRD_N_preds = snoBIRD_N_preds[['chr']]


# Create column for real_gene_id and number of Ns
snoBIRD_N_preds[['real_gene_id', 'number_of_N']] = snoBIRD_N_preds['chr'].str.split('_reduced_', expand=True)
snoBIRD_N_preds['number_of_N'] = snoBIRD_N_preds['number_of_N'].str.strip('N')
snoBIRD_N_preds['number_of_N'] = snoBIRD_N_preds['number_of_N'].astype(int)
snoBIRD_N_preds = snoBIRD_N_preds.drop(columns='chr')

# Drop predictions for which the number of N intrudes the snoRNA sequence (these predictions are not relevant as the snoRNA sequence is shortened)
max_n = dict(zip(positives_info.gene_id, positives_info.half_len))
positives_info = positives_info[positives_info['gene_id'].isin(
                        list(pd.unique(snoBIRD_N_preds['real_gene_id'])))]
positives_info = positives_info.sort_values('half_len', ascending=False)
min_half_sno_len = positives_info['half_len'].min()
snoBIRD_N_preds = snoBIRD_N_preds[snoBIRD_N_preds['number_of_N'] < snoBIRD_N_preds['real_gene_id'].map(max_n)]


# Create 1 column per N number (1: predicted as CD by snoBIRD; 0: not predicted)
df_pivot = snoBIRD_N_preds.pivot_table(index='real_gene_id', columns='number_of_N', aggfunc=lambda x: 1, fill_value=0)
df_pivot = df_pivot.reset_index()

# Add missing column up to the biggest N-mask range that include a snoRNA
missing_cols = [i for i in range(1, int(half_fixed_len - min_half_sno_len + 1)) if i not in df_pivot.columns]
for i in missing_cols:
    df_pivot[i] = 0

# Get the proportion of predicted C/D snoRNA per N content that were added 
# (from 0 to 97N, where 97N reaches the center of all windows from both sides (97*2 = 194nt))
proportion = dict(df_pivot.drop(columns=['real_gene_id']).mean())
proportion = {int(k):v*100 for k,v in proportion.items()}
proportion[0] = 100  # for no N added, all the selected examples were predicted as CD by SnoBIRD
proportion = dict(sorted(proportion.items()))
x = list(proportion.keys())
y = list(proportion.values())


# Get the proportion of snoRNAs that are included or not when reducing the window size
thresholds = np.arange(half_fixed_len, min_half_sno_len -1, -1)  # size thresholds (from 97 which includes all sno to the shortest test set snoRNA half-length)
cumul_prop = [100 * len(positives_info[positives_info['half_len'] < t]) / len(positives_info) for t in thresholds]


# Create the line plot
plt.rcParams['svg.fonttype'] = 'none'
rc = {'ytick.labelsize': 20, 'xtick.labelsize': 20}
plt.rcParams.update(**rc)
fig, ax = plt.subplots(1, 1, figsize=(12, 10))
ax.plot(x, y, color='lightblue')
ax.plot(x, cumul_prop, color='grey')
ax.set_title('Cumulative % of test set C/D snoRNAs that are predicted by SnoBIRD\n'+
            'depending on the number of N-masks on both sides of the window', fontdict={'fontsize': 22})
ax.set_xlabel('Number of N-masks on both sides of the window', fontdict={'fontsize': 22})
ax.set_ylabel('Proportion of test set C/D snoRNAs (%)', fontdict={'fontsize': 22})

legend_list = []
legend_element = mpatches.Patch(color='lightblue', label='% of predicted snoRNAs with N-mask')
legend_element2 = mpatches.Patch(color='grey', label='% of snoRNAs entirely\nincluded below masked length')
legend_list.append(legend_element)
legend_list.append(legend_element2)
plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0.3,0.1),
                fontsize=20)
plt.savefig(snakemake.output.lineplot, dpi=600, bbox_inches='tight')

