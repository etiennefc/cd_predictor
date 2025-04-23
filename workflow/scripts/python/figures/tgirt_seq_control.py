#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
from gtfparse import read_gtf
import sys
import numpy as np
import scipy.stats as stats

# Load dfs and define outputs
tpm_tgirt = snakemake.input.tpm_tgirt
gtfs = snakemake.input.gtf
human_tpm = pd.read_csv(snakemake.input.human_tpm, sep='\t')
output_scatter_tpm = snakemake.output.tpm_skov
output_scatter_tpm_rank = snakemake.output.tpm_skov_rank
output_pie = snakemake.output.pie_biotype
output_df_skov = snakemake.output.df_skov
output_df_chicken = snakemake.output.df_chicken
output_df_macaque = snakemake.output.df_macaque

for path in tpm_tgirt:
    if 'homo_sapiens' in path:
        df_human = pd.read_csv(path, sep='\t')
    elif 'gallus_gallus' in path:
        df_chicken = pd.read_csv(path, sep='\t')
        df_chicken['Avg_gallus_gallus'] = df_chicken[['Chicken_1', 'Chicken_2']].mean(axis=1)
        df_chicken = df_chicken[df_chicken['Avg_gallus_gallus'] > 0 ]
        gtf_chicken = read_gtf([g for g in gtfs if 'gallus_gallus' in g][0]).to_pandas()
        df_chicken = df_chicken.merge(gtf_chicken[['gene_id', 'gene_biotype']], on='gene_id', how='left')
    elif 'macaca_mulatta' in path:
        df_macaque = pd.read_csv(path, sep='\t')
        df_macaque['Avg_macaca_mulatta'] = df_macaque[['Monkey_1', 'Monkey_2']].mean(axis=1)
        df_macaque = df_macaque[df_macaque['Avg_macaca_mulatta'] > 0 ]
        gtf_macaque = read_gtf([g for g in gtfs if 'macaca_mulatta' in g][0]).to_pandas()
        df_macaque = df_macaque.merge(gtf_macaque[['gene_id', 'gene_biotype']], on='gene_id', how='left')



# Filter and merge dfs
human_tpm = human_tpm[['gene_id', 'gene_name', 'gene_biotype', 'gene_biotype2', 'SKOV_1', 'SKOV_2']]
human_tpm['SKOV_avg_frag'] = human_tpm[['SKOV_1', 'SKOV_2']].mean(axis=1)
human_tpm = human_tpm.drop(columns=['SKOV_1', 'SKOV_2'])
human_tpm = human_tpm.merge(df_human[['gene_id', 'SKOV_1']], on='gene_id', how='left')
human_tpm = human_tpm.rename(columns={'SKOV_1': 'SKOV_size_selection'})
human_df = human_tpm.copy()
human_tpm[['SKOV_avg_frag', 'SKOV_size_selection']] = human_tpm[['SKOV_avg_frag', 'SKOV_size_selection']].replace(0, 0.01)  # add pseudocount for log purposes
human_tpm['SKOV_size_selection_log'] = np.log10(human_tpm['SKOV_size_selection'])
human_tpm['SKOV_avg_frag_log'] = np.log10(human_tpm['SKOV_avg_frag'])


# Select only snoRNAs that are > 1 TPM in at least one of the two conditions
sno_human = human_tpm[human_tpm['gene_biotype2'] == 'snoRNA']
sno_human = sno_human[(sno_human[['SKOV_size_selection', 'SKOV_avg_frag_log']] >= 1).any(axis=1)]

# Create TPM scatter plot to compare snoRNA abundance in size-selected vs fragmented SKOV (avg)
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
plt.gca().set_aspect('equal', adjustable='box')
sns.scatterplot(data=sno_human, x='SKOV_size_selection_log', y='SKOV_avg_frag_log', ax=ax)
ax.set_xlabel('SnoRNA abundance in SKOV size-selection (log10(TPM))', fontdict={'fontsize': 12})
ax.set_ylabel('SnoRNA abundance in SKOV ribodepletion fragmented (log10(TPM))', fontdict={'fontsize': 12})
# Calculate r and R²
slope, intercept, r_value, p_value, std_err = stats.linregress(sno_human.SKOV_size_selection_log, sno_human.SKOV_avg_frag_log)
r_squared = r_value**2
plt.text(0.05, 0.95, f"Pearson's r = {r_value:.2f}", ha='left', va='top', 
         transform=plt.gca().transAxes, fontsize=12, color='black')
plt.text(0.05, 0.89, f"p-value = {p_value}", ha='left', va='top', 
         transform=plt.gca().transAxes, fontsize=12, color='black')
plt.text(0.05, 0.83, f"$R^2 = {r_squared:.2f}$", ha='left', va='top', 
         transform=plt.gca().transAxes, fontsize=12, color='black')
plt.margins(0.02)
plt.savefig(output_scatter_tpm, bbox_inches='tight', dpi=600)

print(sno_human[(sno_human['SKOV_avg_frag'] >= 1) & (sno_human['SKOV_size_selection'] < 1)][['gene_id', 'gene_name', 'SKOV_avg_frag', 'SKOV_size_selection']])
print(sno_human[(sno_human['SKOV_avg_frag'] < 1) & (sno_human['SKOV_size_selection'] >= 1)][['gene_id', 'gene_name', 'SKOV_avg_frag', 'SKOV_size_selection']])


# Create TPM rank scatter plot to compare snoRNA abundance in size-selected vs fragmented SKOV (avg)
sno_human = sno_human.sort_values('SKOV_avg_frag', ascending=False).reset_index(drop=True)
sno_human['rank_frag'] = sno_human.index + 1
sno_human = sno_human.sort_values('SKOV_size_selection', ascending=False).reset_index(drop=True)
sno_human['rank_size_selection'] = sno_human.index + 1
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
plt.gca().set_aspect('equal', adjustable='box')
sns.scatterplot(data=sno_human, x='rank_size_selection', y='rank_frag', ax=ax)
# Calculate r and R²
slope, intercept, r_value, p_value, std_err = stats.linregress(sno_human.rank_size_selection, sno_human.rank_frag)
r_squared = r_value**2

plt.text(0.05, 0.95, f"Pearson's r = {r_value:.2f}", ha='left', va='top', 
         transform=plt.gca().transAxes, fontsize=12, color='black')
plt.text(0.05, 0.89, f"p-value = {p_value}", ha='left', va='top', 
         transform=plt.gca().transAxes, fontsize=12, color='black')
plt.text(0.05, 0.83, f"$R^2 = {r_squared:.2f}$", ha='left', va='top', 
         transform=plt.gca().transAxes, fontsize=12, color='black')
plt.margins(0.02)
ax.set_xlabel('SnoRNA abundance (TPM) rank in SKOV size-selection', fontdict={'fontsize': 12})
ax.set_ylabel('SnoRNA abundance (TPM) rank in SKOV ribodepletion fragmented', fontdict={'fontsize': 12})
plt.savefig(output_scatter_tpm_rank, bbox_inches='tight', dpi=600)

# Create simplified biotype column and save dfs
human_df = human_df[human_df['SKOV_size_selection'] > 0]
human_df['gene_biotype2'] = human_df['gene_biotype'].fillna('other').replace({
                'Mt_tRNA': 'tRNA', 'processed_pseudogene': 'pseudogene', 'pre-tRNA': 'tRNA', 'tRNA_fragment': 'tRNA',
                'IG_V_pseudogene': 'pseudogene', 'polymorphic_pseudogene': 'pseudogene', 
                'transcribed_processed_pseudogene': 'pseudogene', 'transcribed_unitary_pseudogene': 'pseudogene', 
                'transcribed_unprocessed_pseudogene': 'pseudogene', 'translated_unprocessed_pseudogene': 'pseudogene',
                'unprocessed_pseudogene': 'pseudogene'})
human_df.loc[human_df['gene_biotype'].isin(['IG_C_gene', 'IG_V_gene', 'TR_C_gene', 'TR_V_gene']), 'gene_biotype2'] = 'protein_coding'
human_df.loc[human_df['gene_biotype'].isin(['Y_RNA', 'scRNA', 'TEC', 'vault_RNA', 'ribozyme']), 'gene_biotype2'] = 'other'
human_df.loc[human_df['gene_biotype'].isin(['rRNA_pseudogene', 'Mt_rRNA', 'ETS-RNA', 'ITS-RNA']), 'gene_biotype2'] = 'rRNA'
human_df = human_df.drop_duplicates('gene_id')
df_macaque['gene_biotype2'] = df_macaque['gene_biotype'].fillna('other').replace({
                'Mt_tRNA': 'tRNA', 'Mt_rRNA': 'rRNA', 'processed_pseudogene': 'pseudogene'})
df_macaque.loc[df_macaque['gene_biotype'].isin(['IG_C_gene', 'IG_V_gene', 'TR_C_gene', 'TR_V_gene']), 'gene_biotype2'] = 'protein_coding'
df_macaque.loc[df_macaque['gene_biotype'].isin(['Y_RNA', 'vault_RNA', 'ribozyme']), 'gene_biotype2'] = 'other'
df_macaque = df_macaque.drop_duplicates('gene_id')
df_chicken['gene_biotype2'] = df_chicken['gene_biotype'].fillna('other').replace({
                'Mt_tRNA': 'tRNA', 'Mt_rRNA': 'rRNA', 'processed_pseudogene': 'pseudogene'})
df_chicken.loc[df_chicken['gene_biotype'].isin(['Y_RNA', 'vault_RNA', 'ribozyme']), 'gene_biotype2'] = 'other'
df_chicken = df_chicken.drop_duplicates('gene_id')

human_df.to_csv(output_df_skov, sep='\t', index=False)
df_macaque.to_csv(output_df_macaque, sep='\t', index=False)
df_chicken.to_csv(output_df_chicken, sep='\t', index=False)

# Count the % of total TPM per biotype
tpm_col = ['SKOV_size_selection', 'Avg_macaca_mulatta', 'Avg_gallus_gallus']
biotypes = sorted(pd.unique(human_df['gene_biotype2']))
print(biotypes)
tpm_percent = []
for i, df in enumerate([human_df, df_macaque, df_chicken]):
    total_tpm = df[tpm_col[i]].sum()
    print(tpm_col[i])
    temp_percent = []
    for bio in biotypes:
        bio_sum = df[df['gene_biotype2'] == bio]
        if len(bio_sum) > 0:
            temp_percent.append(bio_sum[tpm_col[i]].sum() * 100 / total_tpm)
        else:
            temp_percent.append(0)
        if bio == 'rRNA':
            print(f'rRNA: {bio_sum[tpm_col[i]].sum() * 100 / total_tpm}%')
    tpm_percent.append(temp_percent)

# Create pie chart of TPM % per gene_biotype
colors = ['purple', 'green', 'pink', 'teal', 'lightgreen', 'orange', 'blue', 'red', 'black', 'darkred', 'lightblue', 'gold', 'grey']
ft.pie_multiple(1, 3, tpm_percent, biotypes, colors, 
    ['SKOV', 'Macaca mulatta', 'Gallus gallus'], 
    'TPM % per biotype across the size-selected TGIRT-Seq', 'gene_biotype', output_pie)


# Select snoRNAs in human






