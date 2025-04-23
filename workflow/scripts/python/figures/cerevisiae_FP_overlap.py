#!/usr/bin/python3
import pandas as pd
import os
import collections as coll
import seaborn as sns
import matplotlib.pyplot as plt
import functions as ft
from pybedtools import BedTool
import subprocess as sp
import re

# Load bed of FP
FP_df = pd.read_csv(snakemake.input.pos_bed, sep='\t', names=['chr', 
                'start', 'end', 'block_id', 'score', 'strand', 'len'])
print(len(FP_df))
FP_bed = BedTool(snakemake.input.pos_bed)

# Convert gtf to bed
gtf = pd.read_csv(snakemake.input.gtf, sep='\t', skiprows=5, 
                    names=['chr', 'source', 'feature', 'start', 
                            'end', 'score', 'strand', 'frame', 
                            'attributes'])
gtf = gtf[gtf['feature'] == 'gene']

gtf['gene_biotype'] = gtf['attributes'].str.extract(r'gene_biotype "(.+?)";')
gtf['gene_id'] = gtf['attributes'].str.extract(r'gene_id "(.+?)";')
gtf = gtf[['chr', 'start', 'end', 'gene_id', 'score', 'strand', 'gene_biotype']]
gtf.to_csv('gtf_bed.bed', index=False, header=False, sep='\t')
gtf_bed = BedTool('gtf_bed.bed')


# If positive overlaps with any genomic element, return as positives, otherwise intergenic
# enforce strandedness and make sure 100% (f=1) of the positive is included in the genomic element
cols = ['chr', 'start', 'end', 'block_id', 'score', 'strand', 'len', 'chr2', 'start2', 'end2', 'gene_id', 'score2', 'strand2', 'gene_biotype']
pos = FP_bed.intersect(gtf_bed, wb=True, s=True).to_dataframe(names=cols)

FP_df = FP_df[['chr', 'start', 'end', 'block_id', 'score']].merge(pos[['block_id', 'gene_biotype']], how='left', on='block_id')
FP_df['gene_biotype'] = FP_df['gene_biotype'].fillna('intergenic')

#print(FP_df[FP_df['gene_biotype'] == 'transposable_element'][['chr', 'start', 'end']])
print(coll.Counter(FP_df.gene_biotype))

#sp.call('rm gtf_bed.bed fptest.bed', shell=True)


