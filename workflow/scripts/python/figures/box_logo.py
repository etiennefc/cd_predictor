#!/usr/bin/python3
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from math import log2
import numpy as np

""" Create logo of C and D boxes from fasta of either expressed or
    not expressed C/D box snoRNAs."""

# Filter per species
df = pd.read_csv(snakemake.input.box_df, sep='\t')
species_name_dict = snakemake.params.species_short_name
df['species_long'] = df['species_name'].apply(lambda x: species_name_dict[x])
species = snakemake.wildcards.species
df = df[df['species_long'] == species]
if species in ['mus_musculus', 'homo_sapiens', 'drosophila_melanogaster']:
    df = df[df['gene_biotype'] == 'snoRNA_pseudogene']
else:
    df = df[df['gene_biotype'] == 'expressed_CD_snoRNA']

logo_outputs = [snakemake.output.logo_c, snakemake.output.logo_d, 
                snakemake.output.logo_c_prime, snakemake.output.logo_d_prime]
boxes = ['C', 'D', 'C_PRIME', 'D_PRIME']

# Get all box sequences (not sno_id) in a list
# Get the name of group of C/D to redirect figure to correct output
for j, box in enumerate(boxes):
    seqs = list(df[box+'_MOTIF'])

    #Get a count and probability matrix to create the logo
    counts_matrix = logomaker.alignment_to_matrix(seqs)
    prob_matrix = logomaker.transform_matrix(counts_matrix, from_type='counts',
                                            to_type='probability')

    # Create logo wo blanks
    rc = {'ytick.labelsize': 32}
    plt.rcParams.update(**rc)
    plt.rcParams['svg.fonttype'] = 'none'
    logo = logomaker.Logo(prob_matrix, color_scheme='classic')
    logo.ax.set_ylabel("Frequency", fontsize=35)
    logo.ax.set_xlabel(f"Position in {box} box", fontsize=35)
    plt.title(f"{box} box in {species}", fontsize=35)
    plt.savefig(logo_outputs[j], bbox_inches='tight', dpi=600)



