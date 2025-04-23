#!/usr/bin/python3
import pandas as pd
import subprocess as sp 
from pybedtools import BedTool

df = snakemake.input.df
genome_path = snakemake.input.genome_ensembl
sp.call(f'''sed 's/chr//g' {df}'''+''' | awk -v OFS='\t' 'NR>1 {print $2,$3,$4,$6,".",$5}' > temp_ttherm.bed''', shell=True)
bed = BedTool('temp_ttherm.bed')

# get sequence of snoRNAs
fasta = bed.sequence(fi=genome_path, s=True, nameOnly=True)
seq_dict = {}
with open(fasta.seqfn, 'r') as fasta_file:
    for line in fasta_file:
        if '>' in line:
            sno_name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
        else:
            seq = line.strip('\n')
            seq_dict[sno_name] = seq
with open(snakemake.output.fa, 'w') as f:
    for k,v in seq_dict.items():
        f.write(f'>{k}\n{v}\n')
sp.call('rm temp_ttherm.bed', shell=True)
