#!/usr/bin/python3
import pandas as pd

sno_literature = pd.read_csv(snakemake.input.sno_literature, sep='\t')
sno_literature = sno_literature[['gene_id', 'chr', 'strand', 'start', 'end', 'sequence']]

# Load human, mouse and S. cerevisiae dfs (C/D expressed in TGIRT-Seq)
sno_dfs = [sno_literature]
for path in snakemake.input.sno_tgirt:
    df = pd.read_csv(path, sep='\t')
    df = df[['gene_id', 'chr', 'strand', 'start', 'end', 'sequence']]
    sno_dfs.append(df)

# Create fasta of all expressed C/D sno sequences 
all_cd_df = pd.concat(sno_dfs)
all_cd_df.to_csv(snakemake.output.df, sep='\t', index=False)
d = dict(zip(all_cd_df.gene_id, all_cd_df.sequence))

with open(snakemake.output.fa, 'w') as f:
    for id, sequence in d.items():
        f.write(f'>{id}\n')
        f.write(f'{sequence}\n')
