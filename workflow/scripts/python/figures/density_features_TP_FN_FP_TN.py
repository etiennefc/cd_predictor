#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft
import collections as coll

df = pd.read_csv(snakemake.input.feature_df, sep='\t')
pred_df = pd.read_csv(snakemake.input.preds, sep='\t')
feature_cols = ['score_c', 'score_d', 'score_c_prime', 'score_d_prime', 'box_score', 
                'len', 'sno_stability', 'normalized_sno_stability', 'terminal_stem_stability', 'terminal_stem_length_score', 'terminal_stem_combined']
color_dict = snakemake.params.color_conf_val

# Keep only 1 example per sno (no data augmentation)
filtered_pred_df = pred_df.copy()
filtered_pred_df = filtered_pred_df[~filtered_pred_df['gene_id'].str.contains('_[AB][0-9]*$')]

# Merge dfs
df_pred = filtered_pred_df.merge(df, how='left', on='gene_id')

# Add confusion value columns based on pseudogene class
df_pred.loc[(df_pred['y_true'] == 0) & (df_pred['y_pred'] == 0), 'conf_val_pseudo'] = 'TP'
df_pred.loc[(df_pred['y_true'] == 0) & (df_pred['y_pred'] == 1), 'conf_val_pseudo'] = 'FN'
df_pred.loc[(df_pred['y_true'] == 1) & (df_pred['y_pred'] == 0), 'conf_val_pseudo'] = 'FP'
df_pred.loc[(df_pred['y_true'] == 1) & (df_pred['y_pred'] == 1), 'conf_val_pseudo'] = 'TN'

# Create terminal stem combined (terminal stem mfe * terminal length score)
df_pred['terminal_stem_combined'] = df_pred['terminal_stem_stability'] * df_pred['terminal_stem_length_score']


df_pred.loc[df_pred['terminal_stem_combined'] < -100, 'new_pred'] = 1
df_pred.loc[(df_pred['box_score'] <= 2) & (df_pred['terminal_stem_combined'] >= -100), 'new_pred'] = 1
df_pred.loc[(df_pred['normalized_sno_stability'] <= -0.3) & (df_pred['terminal_stem_combined'] >= -100) & (df_pred['box_score'] > 2), 'new_pred'] = 1
df_pred['new_pred'] = df_pred['new_pred'].fillna(0)
df_pred.loc[(df_pred['y_true'] == 0) & (df_pred['new_pred'] == 0), 'new_conf_val_pseudo'] = 'TP'
df_pred.loc[(df_pred['y_true'] == 0) & (df_pred['new_pred'] == 1), 'new_conf_val_pseudo'] = 'FN'
df_pred.loc[(df_pred['y_true'] == 1) & (df_pred['new_pred'] == 0), 'new_conf_val_pseudo'] = 'FP'
df_pred.loc[(df_pred['y_true'] == 1) & (df_pred['new_pred'] == 1), 'new_conf_val_pseudo'] = 'TN'
print(coll.Counter(df_pred.new_conf_val_pseudo))


# Create density
species_tgirt = [['H_sapiens', 'M_musculus', 'D_melanogaster'], 'H_sapiens', 'M_musculus', 'D_melanogaster']
row_titles = species_tgirt.copy()
row_titles[0] = 'All species'
title = 'Feature distribution based on the prediction types of the 2nd model of snoBIRD for the snoRNA pseudogene class across subset of species'
ft.density_multiple(len(species_tgirt), len(feature_cols), df_pred, 'species_name', species_tgirt, 
                        'conf_val_pseudo', feature_cols, color_dict, row_titles, title, snakemake.output.density, 
                        warn_singular=False, hue_order=list(color_dict.keys()))


print(coll.Counter(df_pred['conf_val_pseudo']))
for i in ['H_sapiens', 'M_musculus', 'D_melanogaster']:
        print(i)
        print(coll.Counter(df_pred[df_pred['species_name'] == i]['conf_val_pseudo']))

df_pred[df_pred['conf_val_pseudo'] == 'FP'].to_csv('test_pseudo.tsv', index=False, sep='\t')