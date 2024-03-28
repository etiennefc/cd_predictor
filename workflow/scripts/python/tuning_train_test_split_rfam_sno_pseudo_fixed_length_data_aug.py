#!/usr/bin/python3
import pandas as pd
import collections as coll
import subprocess as sp 
from sklearn.utils import shuffle
from pybedtools import BedTool
import glob

seed = snakemake.params.random_seed
fixed_length = int(snakemake.wildcards.fixed_length)  # fixed window length with 15 nt up/down stream of snoRNA
species_dict = snakemake.params.species_dict
data_aug_num = snakemake.params.data_aug_num
bed_cols = ['chr', 'start', 'end', 'gene_id', 'score', 'strand', 'gene_biotype']
outputs = [snakemake.output.tuning, snakemake.output.training, snakemake.output.test]
# Load dfs
tuning = pd.read_csv(snakemake.input.tuning, sep='\t')
training = pd.read_csv(snakemake.input.training, sep='\t')
test = pd.read_csv(snakemake.input.test, sep='\t')
all_species = list(pd.unique(pd.concat([tuning, training, test])['species_name']))

# For each set and each species, create a bed file of snoRNAs
# Extend the sequence to fixed_length + 5 nt on each side
ext_seq_dict = {}
# 5 nucleotides Before and 5 nt After
suffixes = ['_B5', '_B4', '_B3', '_B2', '_B1', '', '_A1', '_A2', '_A3', '_A4', '_A5']
all_pos = []
for j, df in enumerate([tuning, training, test]):
    df['score'] = '.'
    df['length'] = df.end - df.start
    df['diff'] = fixed_length - df.length  # remaining nt to extend and split betwen the sno 5'/3' ends
    for spe in all_species:
        species_long = species_dict[spe]
        sp_df = df[df['species_name'] == spe]
        chr_size = glob.glob(snakemake.params.chr_size_dir+species_long+'*_chr_size.tsv')[0]
        genome_path = [i for i in snakemake.input.genomes if species_long in i][0]
        for row in sp_df.iterrows():
            row = pd.DataFrame(row[1]).T.reset_index(drop=True) # convert row to dataframe
            if row.loc[0, 'diff'] % 2 == 0: # even number: split the remaining nt equally each side of the sno
                l_extension = int(row['diff'] / 2)
                r_extension = int(l_extension)
            else: # odd number: split the remaining nt almost equally each side of the sno (1 nt more on the 3'end)
                l_extension = int((row['diff'] - 1) / 2)
                r_extension = int(l_extension + 1)
            bed = row[bed_cols]
            bed.to_csv(f'data_aug_bed_{species_long}.tsv', sep='\t', index=False, header=False)
            sno_bed = BedTool(f'data_aug_bed_{species_long}.tsv')
            extended_sno_bed = sno_bed.slop(g=chr_size, l=l_extension + data_aug_num, r=r_extension + data_aug_num)  # extend sequence to obtain fixed_length + 5nt on each side
            extended_sequence_fa = extended_sno_bed.sequence(fi=genome_path, s=True, nameOnly=True)

            with open(extended_sequence_fa.seqfn, 'r') as fa:
                for line in fa:
                    if '>' in line:
                        sno_name = line.strip('\n').replace('>', '').replace('(+)', '').replace('(-)', '')
                    else:
                        seq = line.strip('\n')
                        # Get the 11 sequences for the same snoRNA (i.e. the original fixed_length window 
                        # and 5 windows before and five windows after with a step of 1 nt)
                        data_aug = [seq[0+i:fixed_length+i] for i in range(0,2*data_aug_num+1)]
                        for i, suff in enumerate(suffixes):
                            ext_seq_dict[sno_name+suff] = data_aug[i]
    
    # Duplicate ten times the original snoRNAs and add corresponding adjusted
    # sequence (5 before the original window, the original window and 5 after)
    df_dups_l = []
    for i, suff in enumerate(suffixes):
        temp_df = df.copy()
        temp_df['gene_id'] = temp_df['gene_id'] + suff
        df_dups_l.append(temp_df)
    df_dups = pd.concat(df_dups_l)
    df_dups[f'extended_{fixed_length}nt_sequence'] = df_dups['gene_id'].map(ext_seq_dict)
    df_dups = df_dups.sort_values(by='gene_id').drop(columns=['score', 'length', 'diff']).drop_duplicates(subset=f'extended_{fixed_length}nt_sequence')
    df_dups.to_csv(outputs[j], sep='\t', index=False)
    all_pos.append(df_dups)

    sp.call('rm data_aug*', shell=True)

# Concat all positives
pd.concat(all_pos).to_csv(snakemake.output.all_positives, sep='\t', index=False)
