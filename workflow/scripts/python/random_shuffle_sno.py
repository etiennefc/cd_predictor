#!/usr/bin/python3
import pandas as pd
from sklearn.utils import shuffle

rs = snakemake.params.random_state
output = snakemake.output.shuffled_sno_df
sno_df = pd.read_csv(snakemake.input.sno_sequences, sep='\t')

# Iterate over all expressed C/D and add shuffled sequence as column
seq_dict, shuffled_seq_dict = dict(zip(sno_df.gene_id, sno_df.sequence)), {}
for sno_id, seq in seq_dict.items():
    shuf = ''.join(shuffle(shuffle(seq, random_state=rs), random_state=rs))  # double shuffle
    shuffled_seq_dict[sno_id] = shuf

sno_df['shuffled_seq'] = sno_df.gene_id.map(shuffled_seq_dict)

sno_df.to_csv(output, sep='\t', index=False)