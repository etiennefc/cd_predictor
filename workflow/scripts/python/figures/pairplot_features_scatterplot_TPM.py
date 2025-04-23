#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft

df = pd.read_csv(snakemake.input.feature_df, sep='\t')
human_sno_tpm_df = pd.read_csv(snakemake.input.human_sno_tpm_df, sep='\t')
human_pseudosno_tpm_df = pd.read_csv(snakemake.input.human_pseudosno_tpm_df, sep='\t')
mouse_sno_tpm_df = pd.read_csv(snakemake.input.mouse_sno_tpm_df, sep='\t')
mouse_pseudosno_tpm_df = pd.read_csv(snakemake.input.mouse_pseudosno_tpm_df, sep='\t')
feature_cols = ['score_c', 'score_d', 'score_c_prime', 'score_d_prime', 'box_score', 'gene_biotype', 
                'len', 'sno_stability', 'normalized_sno_stability', 'terminal_stem_stability', 'terminal_stem_length_score']
color_dict = snakemake.params.color_biotype
color_dict = {k: v for k,v in color_dict.items() if k != 'negatives'}

# Create avg tpm column
cols_human = [i for i in human_sno_tpm_df.columns 
            if pd.Series(i).str.contains('_[0-9]$').any()]
cols_mouse = [i for i in mouse_sno_tpm_df.columns 
            if pd.Series(i).str.contains('_[0-9]$').any()]

human_sno_tpm_df['avg_tpm'] = human_sno_tpm_df[cols_human].mean(axis=1)
human_pseudosno_tpm_df['avg_tpm'] = human_pseudosno_tpm_df[cols_human].mean(axis=1)
mouse_sno_tpm_df['avg_tpm'] = mouse_sno_tpm_df[cols_mouse].mean(axis=1)
mouse_pseudosno_tpm_df['avg_tpm'] = mouse_pseudosno_tpm_df[cols_mouse].mean(axis=1)



# Create common df with all features for the pairplot
human_df = pd.concat([human_sno_tpm_df, human_pseudosno_tpm_df])
human_df = human_df.drop(columns='gene_biotype').merge(df[['gene_id'] + feature_cols], how='left', on='gene_id')
mouse_df = pd.concat([mouse_sno_tpm_df, mouse_pseudosno_tpm_df])
mouse_df = mouse_df.drop(columns='gene_biotype').merge(df[['gene_id'] + feature_cols], how='left', on='gene_id')
droso_df = df[df['species_name'] == 'D_melanogaster']

plt.rcParams['svg.fonttype'] = 'none'
sns.pairplot(human_df[feature_cols], hue='gene_biotype', palette=color_dict)
plt.savefig(snakemake.output.pairplot_human)

sns.pairplot(mouse_df[feature_cols], hue='gene_biotype', palette=color_dict)
plt.savefig(snakemake.output.pairplot_mouse)

sns.pairplot(droso_df[feature_cols], hue='gene_biotype', palette=color_dict)
plt.savefig(snakemake.output.pairplot_droso)


# Create scatter plot per biotype and species
ft.scatter_multiple(2, 10, human_df, 'gene_biotype', list(color_dict.keys()), 
        [c for c in feature_cols if c != 'gene_biotype'], 'avg_tpm', list(color_dict.values()), 
        'C/D SnoRNA feature distribution according to average TPM in human', snakemake.output.scatter_tpm_human)

ft.scatter_multiple(2, 10, mouse_df, 'gene_biotype', list(color_dict.keys()), 
        [c for c in feature_cols if c != 'gene_biotype'], 'avg_tpm', list(color_dict.values()), 
        'C/D SnoRNA feature distribution according to average TPM in mouse', snakemake.output.scatter_tpm_mouse)

