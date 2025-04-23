#!/usr/bin/python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import functions as ft
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import collections as coll 

# Merge dfs and simplify species name
df = pd.read_csv(snakemake.input.cd_df, sep='\t')
shap_df = pd.read_csv(snakemake.input.shap_df, sep='\t')
shap_df = shap_df.merge(df[['gene_id', 'species_name', 'gene_biotype', 'sequence']])
shap_df['species_name'] = shap_df['species_name'].replace(snakemake.params.sp_name)
shap_df = shap_df.set_index('gene_id')


# Crate conf_val column based on the pseudogene class
shap_df.loc[(shap_df['predicted_label'] == 'snoRNA_pseudogene') & (shap_df['gene_biotype'] == 'snoRNA_pseudogene'), 
            'conf_val'] = 'TP'
shap_df.loc[(shap_df['predicted_label'] == 'snoRNA_pseudogene') & (shap_df['gene_biotype'] == 'expressed_CD_snoRNA'), 
            'conf_val'] = 'FP'
shap_df.loc[(shap_df['predicted_label'] == 'expressed_CD_snoRNA') & (shap_df['gene_biotype'] == 'expressed_CD_snoRNA'), 
            'conf_val'] = 'TN'
shap_df.loc[(shap_df['predicted_label'] == 'expressed_CD_snoRNA') & (shap_df['gene_biotype'] == 'snoRNA_pseudogene'), 
            'conf_val'] = 'FN'
shap_df = shap_df.drop(columns=['probability', 'predicted_label'])

shap_df2 = shap_df.copy()
shap_df3 = shap_df.copy()

# Create custom cmap
colors = ["#FF0B00", (0, 0, 0), (0, 1, 0)]  # Red, black, green
colors = ["#FF0B00", (1, 1, 1), (0, 1, 0)]  # Red, black, green

n_bins = 200  # Discretize the colormap into 100 bins
cmap_name = 'green_black_red'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)

# Create a clustered heatmap (on the y axis only) to show patterns of shap values across sequence order
# Add a vertical line in the middle of the heatmap to show the symmetry in shap signals for C and D boxes
shap_df = shap_df.drop(columns=['sequence'])
pred_label = shap_df.pop('species_name')
lut = snakemake.params.species_dict
row_colors = pred_label.map(lut)

gene_biotype = shap_df.pop('gene_biotype')
lut2 = {'snoRNA_pseudogene': '#66c2a5', 'expressed_CD_snoRNA': '#f46d43'}
row_colors2 = gene_biotype.map(lut2)

conf_val = shap_df.pop('conf_val')
lut3 = {'TP': 'lightblue', 'FP': 'green', 'TN': 'lightgreen', 'FN': 'blue'}
row_colors3 = conf_val.map(lut3)

clustermap = sns.clustermap(shap_df, col_cluster=False, row_colors=[row_colors, row_colors2, row_colors3], yticklabels=1, xticklabels=1, cmap=cm)
clustermap.ax_heatmap.tick_params(axis='x', labelsize=2)
clustermap.ax_heatmap.tick_params(axis='y', labelsize=2)
cbar = clustermap.ax_heatmap.collections[0].colorbar
cbar.set_label('Avg SHAP values', fontsize=8)
legend_list = []
for sp, color in snakemake.params.species_dict.items():
    legend_element = mpatches.Patch(color=color, label=sp)
    legend_list.append(legend_element)
for g, color in lut2.items():
    legend_element = mpatches.Patch(color=color, label=g)
    legend_list.append(legend_element)
for conf_val, color in lut3.items():
    legend_element = mpatches.Patch(color=color, label=conf_val)
    legend_list.append(legend_element)

plt.legend(handles=legend_list, bbox_to_anchor=(19,1),
                fontsize=6)
clustermap.ax_heatmap.set_title("C/D snoRNAs predicted as expressed or pseudogenes\nand clustered based on their SHAP values profiles", fontsize=12)
clustermap.ax_heatmap.axvline(x=98, color='black', linestyle='--')
plt.savefig(snakemake.output.heatmap_clustered, dpi=500)



# Create a NON-clustered heatmap to show patterns of shap values across sequences of expressed C/D and pseudogenes separately
# Sort the sequences in decreasing order of length to show that C and D motifs are identified in different ranges of length
# Add a vertical line in the middle of the heatmap to show the symmetry in shap signals for C and D boxes
shap_df2['len'] = shap_df2['sequence'].apply(lambda x: len(x))
pseudo = shap_df2[shap_df2['gene_biotype'] == 'snoRNA_pseudogene'].sort_values(by=['conf_val', 'len'], ascending=[False, False])
expressed = shap_df2[shap_df2['gene_biotype'] == 'expressed_CD_snoRNA'].sort_values(by=['conf_val', 'len'], ascending=[False, False])
shap_df2 = pd.concat([expressed, pseudo])

gene_biotype = shap_df2.pop('gene_biotype')
row_colors2 = gene_biotype.map(lut2)
#conf_val = shap_df2.pop('conf_val')
#lut3 = {'TP': 'lightblue', 'FP': 'green', 'TN': 'lightgreen', 'FN': 'blue'}
#row_colors3 = conf_val.map(lut3)

shap_df2 = shap_df2.drop(columns=['len', 'sequence', 'species_name', 'conf_val'])
clustermap = sns.clustermap(shap_df2, col_cluster=False, row_cluster=False, row_colors=[row_colors2], yticklabels=1, xticklabels=1, cmap=cm)
clustermap.ax_heatmap.tick_params(axis='x', labelsize=2, width=0.01)
clustermap.ax_heatmap.tick_params(axis='y', labelsize=2, width=0.01)
cbar = clustermap.ax_heatmap.collections[0].colorbar
cbar.set_label('Avg SHAP values', fontsize=8)
legend_list = []
#for conf_val, color in lut3.items():
#    legend_element = mpatches.Patch(color=color, label=conf_val)
#    legend_list.append(legend_element)
for g, color in lut2.items():
    legend_element = mpatches.Patch(color=color, label=g)
    legend_list.append(legend_element)
plt.legend(handles=legend_list, bbox_to_anchor=(19,1),
                fontsize=6)
clustermap.ax_heatmap.set_title("C/D snoRNAs predicted as expressed or pseudogenes\nand sorted by decreasing length per biotype and confusion value", fontsize=12)
clustermap.ax_heatmap.axvline(x=98, color='black', linestyle='--')
plt.savefig(snakemake.output.heatmap, dpi=500)


# Create a NON-clustered heatmap to show patterns of shap values across sequences of expressed C/D and pseudogenes separately
# Sort the sequences per species
# Add a vertical line in the middle of the heatmap to show the symmetry in shap signals for C and D boxes
pseudo = shap_df3[shap_df3['gene_biotype'] == 'snoRNA_pseudogene'].sort_values(by='species_name')
expressed = shap_df3[shap_df3['gene_biotype'] == 'expressed_CD_snoRNA'].sort_values(by='species_name')
shap_df3 = pd.concat([expressed, pseudo])
pred_label = shap_df3.pop('species_name')
row_colors = pred_label.map(lut)

gene_biotype = shap_df3.pop('gene_biotype')
row_colors2 = gene_biotype.map(lut2)
shap_df3 = shap_df3.drop(columns=['sequence', 'conf_val'])
clustermap = sns.clustermap(shap_df3, col_cluster=False, row_cluster=False, row_colors=[row_colors, row_colors2], yticklabels=1, xticklabels=1, cmap=cm)
clustermap.ax_heatmap.tick_params(axis='x', labelsize=2, width=0.01)
clustermap.ax_heatmap.tick_params(axis='y', labelsize=2, width=0.01)
cbar = clustermap.ax_heatmap.collections[0].colorbar
cbar.set_label('Avg SHAP values', fontsize=8)
legend_list = []
for sp, color in snakemake.params.species_dict.items():
    legend_element = mpatches.Patch(color=color, label=sp)
    legend_list.append(legend_element)
for g, color in lut2.items():
    legend_element = mpatches.Patch(color=color, label=g)
    legend_list.append(legend_element)
plt.legend(handles=legend_list, bbox_to_anchor=(19,1),
                fontsize=6)
clustermap.ax_heatmap.set_title("C/D snoRNAs predicted as expressed or pseudogenes\nand sorted by species and per gene biotype", fontsize=12)
clustermap.ax_heatmap.axvline(x=98, color='black', linestyle='--')
plt.savefig(snakemake.output.heatmap_species, dpi=500)


# TO create SHAP line plot of all examples
'''
df = pd.read_csv('data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_194nt.tsv', sep='\t')
shap_df2 = shap_df2.merge(df[['gene_id', 'extended_194nt_sequence']], on='gene_id', how='left')

other = shap_df2[shap_df2['predicted_label'] == 'CD_snoRNA']
seq_dict = dict(zip(other.gene_id, other.extended_194nt_sequence))
for sno, seq_ in seq_dict.items():
    xticklabels = ['CLS'] + [i for i in seq_] + ['SEP']
    avg_shap_per_nt = list(other[other['gene_id'] == sno].values.tolist()[0][3:-1])
    ft.shap_lineplot(range(len(avg_shap_per_nt)), avg_shap_per_nt, xticklabels, 'Sequence tokens', 
                    'Average SHAP values', f'{sno} predicted as CD_snoRNA', f'results/figures/lineplot/shap/snoBIRD/CD_snoRNA_{sno}.svg')
'''
