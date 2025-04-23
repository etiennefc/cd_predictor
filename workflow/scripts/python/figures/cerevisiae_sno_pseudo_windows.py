#!/usr/bin/python3
import pandas as pd
import os
import collections as coll
import seaborn as sns
import matplotlib.pyplot as plt
import functions as ft
from pybedtools import BedTool
import subprocess as sp

window_size = int(snakemake.params.fixed_length)

# Load pos_blocks and create strand dict
pos_blocks = BedTool(snakemake.input.pos_blocks)
pos_blocks_df = pos_blocks.to_dataframe()
strand = dict(zip(pos_blocks_df['name'], pos_blocks_df['strand']))


# Create 190 nt windows 
windows = pos_blocks.window_maker(w=window_size, s=1, b=pos_blocks, i='src').to_dataframe()

# Add strand and keep only 190 nt window (not smaller ones created at the end of each interval)
windows['len'] = windows['end'] - windows['start']
windows = windows[windows['len'] == 190]
windows['strand'] = windows['name'].map(strand)
windows['len'] = windows['len'].replace(190, '.')


# Modify block_id so that windwos in the same block have an added suffix _1, _2, ...
windows['gene_id'] = windows.groupby('name').cumcount() + 1
windows['gene_id'] = windows['name'] + '_' + windows['gene_id'].astype(str)
windows = windows[['chrom', 'start', 'end', 'gene_id', 'len', 'strand']]

windows.to_csv('temp_window.bed', sep='\t', index=False, header=False)
filtered_windows = BedTool('temp_window.bed')

# Get sequence of these filtered windows
fasta = filtered_windows.sequence(fi=snakemake.input.genome, nameOnly=True, s=True)
seq_dict = {}
with open(fasta.seqfn, 'r') as fasta_file:
    for line in fasta_file:
        if '>' in line:
            block_name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
        else:
            seq = line.strip('\n')
            seq_dict[block_name] = seq

windows['seq'] = windows['gene_id'].map(seq_dict)
windows.to_csv(snakemake.output.df, sep='\t', index=False)


sp.call('rm temp_window.bed', shell=True)