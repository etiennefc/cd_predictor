#!/usr/bin/python3
import pandas as pd
import numpy as np
import itertools as it
from pybedtools import BedTool
import subprocess as sp

cd_id_rnacentral = pd.read_csv(snakemake.input.cd_id_rnacentral, sep='\t', 
        names=['rnacentral_id'])  # downlaoded from RNAcentral directly
cols = ['chr', 'start', 'end', 'rnacentral_id', 'score', 'strand', 'start2', 
        'end2', 'other', 'num', 'other2', 'other3', 'dot', 'gene_biotype', 'source']
ncRNA_bed = pd.read_csv(snakemake.input.ncRNA_rnacentral_bed, sep='\t', names=cols)

# Filter ncRNA bed from RNAcentral to keep only C/D snoRNAs
sno_bed = ncRNA_bed[ncRNA_bed['gene_biotype'] == 'snoRNA']
cd_bed = sno_bed[sno_bed['rnacentral_id'].isin(cd_id_rnacentral.rnacentral_id)].sort_values('rnacentral_id')
cd_bed['sno_type'] = 'C/D'
cd_bed['gene_id'] = cd_bed.groupby('rnacentral_id').cumcount() + 1
cd_bed['gene_id'] = cd_bed['rnacentral_id'] +'_'+ cd_bed['gene_id'].astype(str).where(cd_bed['gene_id'] > 1, '')
cd_bed['gene_id'] = cd_bed['gene_id'].str.replace('_$', '', regex=True)

cd_bed = cd_bed[['rnacentral_id', 'chr', 'start', 'end', 'strand', 'gene_id', 'sno_type']]

cd_bed.to_csv(snakemake.output.df, sep='\t', index=False)

print(cd_bed)
