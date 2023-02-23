#!/usr/bin/python3
import pandas as pd

output = snakemake.output.random_ncRNA
tRNA_df = pd.read_csv(snakemake.input.tRNA, sep='\t')

species = ['Hsapi38', 'Mmusc39', 'Scere3', 'Osati7', 
            'Athal10', 'Lmajo2', 'Dmela6', 'Celeg11']  

# O. tauri, giardia, discoideum are not there (G. gallus, T. thermophila, M. mulatta & O. anatinus are there but not the latest version)

# Keep only relevant species and remove header line for each species
tRNA_df = tRNA_df[(tRNA_df['GenomeID'].isin(species)) & (tRNA_df['Score'] != 'pseudo')]

print(tRNA_df)
import collections as coll
print(coll.Counter(tRNA_df.GenomeID))
print(tRNA_df.columns)