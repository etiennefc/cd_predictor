#!/usr/bin/python3
import pandas as pd
import numpy as np
import itertools as it
import subprocess as sp

if snakemake.wildcards.species == 'plasmodium_falciparum':
    pd.DataFrame(columns=
                ['gene_id', 'gene_name', 'chr', 'strand', 'start', 'end']
                ).to_csv(snakemake.output.df, sep='\t', index=False)
else:
    cols = ['match', 'mismatch', 'rep_match', 'Ns', 'query_gap_count', 'query_gap_bases',
            'target_gap_count', 'target_gap_bases', 'strand', 'gene_name', 'query_size',
            'query_start', 'query_end', 'chr', 'target_size', 'start', 'end', 'block_count', 
            'block_sizes', 'query_starts', 'target_starts']
    # where query is a snoRNA and target is a location in the species genome
    df = pd.read_csv(snakemake.input.blat, sep='\t', skiprows=5, names=cols)

    # Custom filter to remove snoRNAs not on real assembled chromosomes **TO UPDATE FOR OTHER SPECIES
    df['chr'] = df['chr'].astype(str)
    df = df[~df['chr'].str.contains('QNV|RZJ')]

    # Remove alignments that have more than 5 nt gaps in the target sequence
    df = df[df['target_gap_bases'] <= 5]

    # Keep only the match with the highest matched nt per gene_name 
    df = df.loc[df.groupby('gene_name').match.idxmax()]

    # Remove duplicate snoRNAs with the same start/end (keep only the one with highest match)
    df = df.loc[df.groupby(['chr', 'strand','start', 'end']).match.idxmax()]
    df = df.drop_duplicates(['chr', 'strand', 'start', 'end'])
    df = df[['gene_name', 'chr', 'strand', 'start', 'end']].reset_index(drop=True)

    final_df = df.copy()
    final_df['gene_id'] = final_df['gene_name']

    # Add genome_version, species name and classification, validation method and article source columns
    final_df = final_df[['gene_id', 'gene_name', 'chr', 'strand', 'start', 'end']]
    final_df['start'] = final_df['start'] + 1
    print(final_df)
    final_df.to_csv(snakemake.output.df, sep='\t', index=False)
