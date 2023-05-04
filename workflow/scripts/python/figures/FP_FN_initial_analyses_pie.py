#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 


snoreport = pd.read_csv(snakemake.input.snoreport, sep='\t')
snoscan = pd.read_csv(snakemake.input.snoscan, sep='\t')
infernal_rfam = pd.read_csv(snakemake.input.infernal_rfam, sep='\t')
pie_species = snakemake.output.pie_species
pie_neg_type = snakemake.output.pie_neg_type
error = snakemake.wildcards.error

# Count the FP that are part of a given species/gene_biotype
mods = ['snoreport2', 'snoscan', 'infernal_rfam']
species = sorted(list(pd.unique(snoreport.species_name)))
biotype = sorted(list(pd.unique(snoreport.gene_biotype)))
ax_titles = []
species_count, neg_type_count = [], []
for i, df in enumerate([snoreport, snoscan, infernal_rfam]):
    if error == 'FP':
        fp = df[(df['target'] == 'other') & (df[f'{mods[i]}_prediction'] == 'expressed_CD_snoRNA')]
        err = 'false positives'
    elif error == 'FN':
        fp = df[(df['target'] == 'expressed_CD_snoRNA') & (df[f'{mods[i]}_prediction'] == 'other')]
        err = 'false negatives'
    print(fp[['species_name', 'gene_biotype', 'gene_id']])
    ax_titles.append(f'{mods[i]} ({len(fp)})')
    temp_l, temp_l2 = [], []
    for sp in species:
        nb = len(fp[fp['species_name'] == sp])
        temp_l.append(nb)
    species_count.append(temp_l)
    for ty in biotype:
        nb2 = len(fp[fp['gene_biotype'] == ty])
        temp_l2.append(nb2)
    neg_type_count.append(temp_l2)

# Create pie chart for species distribution in FP
sp_colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
            '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
            'lightgrey','grey','black']
ft.pie_multiple(1, len(mods), species_count, species, sp_colors, ax_titles, 
                f'Species distribution across the {err}', 
                '', pie_species)
print(biotype)
# Create pie chart for gene_biotype distribution in FP
ty_colors = ['#d73027','orange','pink','lightgrey','grey',
            'black','green','blue','#4575b4']
ft.pie_multiple(1, len(mods), neg_type_count, biotype, ty_colors, ax_titles, 
                f'Negative type distribution across the {err}', 
                '', pie_neg_type)