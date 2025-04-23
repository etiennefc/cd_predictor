#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 

# Load dfs
dfs, dfs2 = [], []
for p in snakemake.input.sets_first_model:
    dfs.append(pd.read_csv(p, sep='\t'))
for p in snakemake.input.sets_second_model:
    dfs2.append(pd.read_csv(p, sep='\t'))

df1 = pd.concat(dfs)
df2 = pd.concat(dfs2)

# Select only the real examples, not the data augmented
positives = df1[~df1['gene_id'].str.contains('_[AB][0-9]+$')]
negatives = df1[df1['target'] == 'other']


# Convert species short name to long name
sp_name = snakemake.params.species_name
positives['species_name'] = positives['species_name'].replace(sp_name)
negatives['species_name'] = negatives['species_name'].replace(sp_name)

exp_cd = positives[positives['target'] == 'expressed_CD_snoRNA']
snoRNA_pseudo = positives[positives['target'] == 'snoRNA_pseudogene']

pie_species = snakemake.output.pie_species
pie_neg_type = snakemake.output.pie_neg_type
pie_neg_species = snakemake.output.pie_neg_species
pie_pseudo_species = snakemake.output.pie_pseudo_species
sp_colors = snakemake.params.species_colors
ty_colors = snakemake.params.biotype_colors



# Count the number of examples that are part of a given species/gene_biotype
species = ['homo_sapiens', 'mus_musculus', 'drosophila_melanogaster', 'ornithorhynchus_anatinus', 
            'caenorhabditis_elegans', 'macaca_mulatta', 'gallus_gallus', 
            'tetrahymena_thermophila', 'dictyostelium_discoideum', 'giardia_lamblia', 'leishmania_major', 
            'neurospora_crassa', 'saccharomyces_cerevisiae', 'candida_albicans', 'aspergillus_fumigatus',
            'arabidopsis_thaliana', 'oryza_sativa', 'ostreococcus_tauri']
biotype = ['random_exonic_region', 'random_intronic_region', 'random_intergenic_region',
            'tRNA', 'HACA_snoRNA', 'snRNA', 'pre_miRNA', 'shuffled_expressed_CD_snoRNA']
print(species)
print(biotype)
species_colors = [sp_colors[sp] for sp in species]
biotype_colors = [ty_colors[biot] for biot in biotype]
species_count, neg_type_count = [], []
sp_neg, pseudo_species_count = [], []
for sp in species:
    nb = len(exp_cd[exp_cd['species_name'] == sp])
    nb_sp_neg = len(negatives[negatives['species_name'] == sp])
    pseudo_nb = len(snoRNA_pseudo[snoRNA_pseudo['species_name'] == sp])
    species_count.append(nb)
    sp_neg.append(nb_sp_neg)
    pseudo_species_count.append(pseudo_nb)

for ty in biotype:
    nb2 = len(negatives[negatives['gene_biotype'] == ty])
    neg_type_count.append(nb2)

# Create pie chart for species distribution in expressed CD
ft.donut(species_count, species, species_colors, '', 
                '', '', pie_species)
# Create pie chart for species distribution in snoRNA_pseudogene
ft.donut(pseudo_species_count, species, species_colors, '', 
                '', '', pie_pseudo_species)                
# Create pie chart for gene_biotype distribution in negatives
ft.donut(neg_type_count, biotype, biotype_colors, '', 
                '', '', pie_neg_type)
# Create pie chart for species distribution in negatives
ft.donut(sp_neg, species, species_colors, '', 
                '', '', pie_neg_species)

print(exp_cd)
print(snoRNA_pseudo)
print(negatives)

