#!/usr/bin/python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from upsetplot import UpSet
import functions as ft

""" Generate an upset plot per confusion value (true positive, false positive,
    false negative or true negative) to see the intersection in the snoRNAs
    (mis)classified by all models (3 existing predictors and GRU NN)."""

sp_color_dict = snakemake.params.species_colors
biot_color_dict = snakemake.params.biotype_colors
model_colors = snakemake.params.models_colors
model_colors['else'] = 'black'

# Load predictions for each tool
snoreport = pd.read_csv([p for p in snakemake.input.other_tool if 'snoreport' in p][0], sep='\t')
snoscan = pd.read_csv([p for p in snakemake.input.other_tool if 'snoscan' in p][0], sep='\t')
infernal_rfam = pd.read_csv([p for p in snakemake.input.other_tool if 'infernal' in p][0], sep='\t')
infernal_rfam['infernal_rfam_prediction'] = infernal_rfam['infernal_rfam_prediction'].replace('expressed_CD_snoRNA', 'CD_snoRNA')
snoscan['snoscan_prediction'] = snoscan['snoscan_prediction'].replace('expressed_CD_snoRNA', 'CD_snoRNA')
snoreport['snoreport2_prediction'] = snoreport['snoreport2_prediction'].replace('expressed_CD_snoRNA', 'CD_snoRNA')

snoBIRD = pd.read_csv(snakemake.input.snoBIRD, sep='\t').rename(columns={'y_true': 'target', 'y_pred':'snoBIRD_prediction'})
snoBIRD['target'] = snoBIRD['target'].replace({0: 'other', 1: 'CD_snoRNA'})
snoBIRD['snoBIRD_prediction'] = snoBIRD['snoBIRD_prediction'].replace({0: 'other', 1: 'CD_snoRNA'})


# Merged dfs
merged_df = snoreport[['gene_id', 'gene_biotype', 
                    'species_name', 'snoreport2_prediction']].merge(snoscan[['gene_id', 'snoscan_prediction']], 
                    how='left', on='gene_id')
merged_df = merged_df.merge(infernal_rfam[['gene_id', 'infernal_rfam_prediction']], how='left', on='gene_id')
merged_df = merged_df.merge(snoBIRD[['gene_id', 'snoBIRD_prediction', 'target']], how='left', on='gene_id')
print(merged_df)





# Create conf_value column
mods = ['snoBIRD', 'snoreport2', 'snoscan', 'infernal_rfam']
for model in mods:
    merged_df.loc[(merged_df['target'] == 'CD_snoRNA') & (merged_df[f'{model}_prediction'] == 'CD_snoRNA'), f'conf_val_{model}'] = 'TP'
    merged_df.loc[(merged_df['target'] == 'CD_snoRNA') & (merged_df[f'{model}_prediction'] != 'CD_snoRNA'), f'conf_val_{model}'] = 'FN'
    merged_df.loc[(merged_df['target'] == 'other') & (merged_df[f'{model}_prediction'] != 'other'), f'conf_val_{model}'] = 'FP'
    merged_df.loc[(merged_df['target'] == 'other') & (merged_df[f'{model}_prediction'] == 'other'), f'conf_val_{model}'] = 'TN'


positives = merged_df[merged_df['target'] == 'CD_snoRNA']
negatives = merged_df[merged_df['target'] == 'other']

# 2 fake FP (predicted by all models including snoBIRD), they're miRNA that overlap with a C/D snoRNA
#print(negatives[(negatives['infernal_rfam_prediction'] == 'CD_snoRNA') & (negatives['target'] == 'other') & (negatives['infernal_rfam_prediction'] == 'CD_snoRNA') & (negatives['snoreport2_prediction'] == 'CD_snoRNA') & (negatives['snoscan_prediction'] == 'CD_snoRNA')])


# Set index for TP and create upset
TP_FN = positives.set_index(positives.conf_val_snoBIRD == 'TP').set_index(
                positives.conf_val_snoreport2 == 'TP', append=True).set_index(
                positives.conf_val_snoscan == 'TP', append=True).set_index(
                positives.conf_val_infernal_rfam == 'TP', append=True)

plt.rcParams['svg.fonttype'] = 'none'
upset = UpSet(TP_FN, show_counts=True, sort_by='-degree', sort_categories_by='-cardinality')
upset.plot()
plt.savefig(snakemake.output.TP_FN_upset, dpi=600, bbox_inches='tight')

# Set index for TN and create upset
TN_FP = negatives.set_index(negatives.conf_val_snoBIRD == 'TN').set_index(
                negatives.conf_val_snoreport2 == 'TN', append=True).set_index(
                negatives.conf_val_snoscan == 'TN', append=True).set_index(
                negatives.conf_val_infernal_rfam == 'TN', append=True)

plt.rcParams['svg.fonttype'] = 'none'
upset = UpSet(TN_FP, show_counts=True, sort_by='-degree', sort_categories_by='-cardinality')
upset.plot()
plt.savefig(snakemake.output.TN_FP_upset, dpi=600, bbox_inches='tight')


# Create donut plots for all tools across confusion value with a hue of species or gene_biotype
ax_titles, sp_counts, biot_counts = [], [], []
for model in mods:
    for conf_val in ['TP', 'FN', 'TN', 'FP']:
        temp_df = merged_df[merged_df[f'conf_val_{model}'] == conf_val]
        temp_sp, temp_biot = [], []
        for sp in sp_color_dict.keys():
            temp_sp.append(len(temp_df[temp_df['species_name'] == sp]))
        for biot in biot_color_dict.keys():
            temp_biot.append(len(temp_df[temp_df['gene_biotype'] == biot]))
        sp_counts.append(temp_sp)
        biot_counts.append(temp_biot)
        ax_titles.append(f'{conf_val}_{model}')



ft.pie_multiple(4, 4, sp_counts, list(sp_color_dict.keys()), list(sp_color_dict.values()), ax_titles, 
        'Species distribution across prediction and tools', 'Species', snakemake.output.donut_conf_value_species)
ft.pie_multiple(4, 4, biot_counts, list(biot_color_dict.keys()), list(biot_color_dict.values()), ax_titles, 
        'Biotype distribution across prediction and tools', 'Species', snakemake.output.donut_conf_value_biotype)


