#!/usr/bin/python3
import pandas as pd

# Load test df and keep only the accurately predicted snoRNAs by snoBIRD
fixed_length = int(snakemake.params.fixed_length)
test_df = pd.read_csv(snakemake.input.test_set, sep='\t')
snoBIRD = pd.read_csv(snakemake.input.snoBIRD_test_set, sep='\t')

test_df = test_df[test_df['gene_id'].isin(list(snoBIRD['chr']))]

# Create masked windows from the test set positive examples by adding 
# 1 N nucleotide to replace each first-last, second-second last, etc. pairs up 
# until all the initial sequence is replaced by Ns
def N_replace(seq_, x_):
    # Replace the sequence by x_ N nt from the start and end of the sequence
    # ex: X_ = 2; seq_ = 'ATCGGAA' --> 'NNCGGNN'
    final_seq = 'N' * x_ + seq_[x_:-x_] + 'N' * x_ 

    return final_seq

sequences = {}
for i in range(1, int(fixed_length/2 + 1)):
    test_df[f'window_reduced_{i}N'] = test_df[f'extended_{fixed_length}nt_sequence'].apply(
                                                    lambda x: N_replace(x, i))
    for j in range(len(test_df)):
        id_ = test_df['gene_id'].iloc[j] + f'_reduced_{i}N'
        sequences[id_] = test_df[f'window_reduced_{i}N'].iloc[j]

# Create output fasta out of it 
with open(snakemake.output.windows, 'w') as file_:
    for k,v in sequences.items():
        file_.write(f'>{k}\n{v}\n') 
    
    