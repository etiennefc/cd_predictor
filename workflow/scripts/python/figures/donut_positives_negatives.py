#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 


positives = pd.read_csv(snakemake.input.positives, sep='\t')
negatives = pd.concat([pd.read_csv(path, sep='\t') for path in snakemake.input.negatives])
pie_species = snakemake.output.pie_species
pie_neg_type = snakemake.output.pie_neg_type
sp_colors = snakemake.params.species_colors
ty_colors = snakemake.params.biotype_colors

# Convert species short name to long name
sp_name = snakemake.params.species_name
positives['species_name'] = positives['species_name'].replace(sp_name)

# Remove snoRNA pseudogenes from negatives
negatives = negatives[negatives['gene_biotype'] != 'snoRNA_pseudogene']

# Count the number of examples that are part of a given species/gene_biotype
species = ['gallus_gallus', 'macaca_mulatta', 'drosophila_melanogaster', 'caenorhabditis_elegans', 
            'ornithorhynchus_anatinus', 'mus_musculus', 'homo_sapiens',
            'tetrahymena_thermophila', 'dictyostelium_discoideum', 'giardia_lamblia', 'leishmania_major', 
            'saccharomyces_cerevisiae', 'arabidopsis_thaliana', 
            'oryza_sativa', 'ostreococcus_tauri']
biotype = ['random_exonic_region', 'random_intergenic_region', 'random_intronic_region',
            'shuffled_expressed_CD_snoRNA', 'HACA_snoRNA', 'pre_miRNA', 'snRNA', 'tRNA']
print(species)
print(biotype)
species_colors = [sp_colors[sp] for sp in species]
biotype_colors = [ty_colors[biot] for biot in biotype]
species_count, neg_type_count = [], []

for sp in species:
    nb = len(positives[positives['species_name'] == sp])
    species_count.append(nb)
for ty in biotype:
    nb2 = len(negatives[negatives['gene_biotype'] == ty])
    neg_type_count.append(nb2)

# Create pie chart for species distribution in FP
ft.donut(species_count, species, species_colors, '', 
                '', '', pie_species)
# Create pie chart for gene_biotype distribution in FP
ft.donut(neg_type_count, biotype, biotype_colors, '', 
                '', '', pie_neg_type)