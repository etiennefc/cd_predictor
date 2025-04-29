#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
from glob import glob 

colors = snakemake.params.color_dict
colors = {k:v for k,v in colors.items() if k != 'arabidopsis_thaliana'}
pred_path = snakemake.input.snoBIRD_preds
all_spe = ['schizosaccharomyces_pombe', 'drosophila_melanogaster', 'tetrahymena_thermophila', 
            'danio_rerio', 'gallus_gallus', 'plasmodium_falciparum', 
            'macaca_mulatta',  'homo_sapiens']

# Load prediction dfs
dfs = []
for p in pred_path:
    df = pd.read_csv(glob(p)[0], sep='\t')
    species_ = p.split('/')[4]
    df['species'] = species_
    dfs.append(df[['gene_id', 'predicted_label', 'species']])
    

all_cd = pd.concat(dfs)
all_cd['predicted_label'] = all_cd['predicted_label'].replace('CD_snoRNA_pseudogene', 'snoRNA_pseudogene')
print(all_cd)


# Nb of predicted expressed snoRNA and snoRNA pseudogene per species
exp_status = list(pd.unique(all_cd.predicted_label))
counts, sno_nb = [], []
for sp in all_spe:
    c = []
    sp_df = all_cd[all_cd['species'] == sp]
    sno_nb.append(len(sp_df))
    for status in exp_status:
        c.append(len(sp_df[sp_df['predicted_label'] == status]))
    counts.append(c)
percent_counts = ft.percent_count(counts)

# Change long species name for shorter one
short_name_dict = snakemake.params.species_short_name
short_names = [k for v in all_spe for k, val in short_name_dict.items() if val == v]
# Create bar chart
ft.stacked_bar2(percent_counts, short_names,
                exp_status, 'SnoBIRD predictions across species', 'Species name', 
                'Proportion of snoBIRD predictions (%)', colors, 0, 105, 
                [f'{i}' for i in sno_nb], snakemake.output.bar) 
