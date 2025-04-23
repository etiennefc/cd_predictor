#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft
import subprocess as sp

shap_df = pd.read_csv(snakemake.input.shap_df, sep='\t')
cd_df = pd.read_csv(snakemake.input.cd_df, sep='\t')
fixed_length = snakemake.params.fixed_length

# Merge target value column to shap df
shap_df = shap_df.merge(cd_df[['gene_id', 'gene_biotype', f'extended_{fixed_length}nt_sequence']], how='left', on='gene_id')


# Add confusion value columns based on pseudogene class
shap_df.loc[(shap_df['gene_biotype'] == 'snoRNA_pseudogene') & (shap_df['predicted_label'] == 'snoRNA_pseudogene'), 
					'conf_val_pseudo'] = 'TP'
shap_df.loc[(shap_df['gene_biotype'] == 'snoRNA_pseudogene') & (shap_df['predicted_label'] == 'expressed_CD_snoRNA'), 
					'conf_val_pseudo'] = 'FN'
shap_df.loc[(shap_df['gene_biotype'] == 'expressed_CD_snoRNA') & (shap_df['predicted_label'] == 'snoRNA_pseudogene'), 
					'conf_val_pseudo'] = 'FP'
shap_df.loc[(shap_df['gene_biotype'] == 'expressed_CD_snoRNA') & (shap_df['predicted_label'] == 'expressed_CD_snoRNA'), 
					'conf_val_pseudo'] = 'TN'

# Create line plot for all examples of a given confusion value

for v in ['TP', 'TN', 'FP', 'FN']:
	path = [p for p in snakemake.output if v in p][0]
	df = shap_df[shap_df['conf_val_pseudo'] == v]
	seq_dict = dict(zip(df.gene_id, df[f'extended_{fixed_length}nt_sequence']))
	sp.call(f'mkdir -p {path}', shell=True)
	for sno, seq_ in seq_dict.items():
		print(v, sno)
		xticklabels = ['CLS'] + [i for i in seq_] + ['SEP']
		avg_shap_per_nt = list(df[df['gene_id'] == sno].values.tolist()[0][3:-4])
		ft.shap_lineplot(range(len(avg_shap_per_nt)), avg_shap_per_nt, xticklabels, 'Sequence tokens', 
				'Average SHAP values', f'{sno} predicted as a {v}', f'{path}/{sno}.svg')