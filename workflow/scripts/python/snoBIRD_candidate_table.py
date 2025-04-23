#!/usr/bin/python3
import pandas as pd
from glob import glob
import re

candidates = snakemake.input.new_cd_df

# Load dfs
snoBIRD = {}
for p in snakemake.input.snoBIRD:
    sp = p.split('/')[4]
    df = pd.read_csv(glob(p+'/*tsv')[0], sep='\t', dtype = {'chr': 'str'})
    df['species'] = sp
    snoBIRD[sp] = df

snoBIRD_cols = ['gene_id', 'species', 'chr', 'start', 'end', 'strand', 'probability_CD', 
                'box_score', 'C_MOTIF', 'C_START', 'C_END', 'D_MOTIF', 'D_START', 
                'D_END', 'C_PRIME_MOTIF', 'C_PRIME_START', 'C_PRIME_END', 
                'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END', 'terminal_stem_score', 
                'normalized_sno_stability', 'probability_expressed_pseudogene', 
                'predicted_label', 'predicted_sequence']


final_dfs = []
#for p in glob('results/predictions/snoBIRD/bedgraph_overlap/multiple_filters/*_potential_new_cd_snoBIRD_preds.tsv'):
for p in candidates:
    if not re.search(r'drosophila|plasmodium|arabidopsis', p):
        sp_ = p.split('/')[-1].split('_potential_')[0]
        candidate_species = pd.read_csv(p, sep='\t')
        candidate_species = candidate_species[['gene_id', 'avg_coverage_samples']]
        snoBIRD_info_df = snoBIRD[sp_]
        candidate_species = candidate_species.merge(
                                                snoBIRD_info_df[snoBIRD_cols], 
                                                how='left', on='gene_id')
        final_dfs.append(candidate_species)

# Save final df
final_cols = ['prediction_id', 'species', 'mean_coverage', 'chr', 'start', 'end',
       'strand', 'probability_CD', 'box_score', 'C_MOTIF', 'C_START', 'C_END',
       'D_MOTIF', 'D_START', 'D_END', 'C_PRIME_MOTIF', 'C_PRIME_START',
       'C_PRIME_END', 'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END',
       'terminal_stem_score', 'normalized_sno_stability',
       'probability_expressed_pseudogene', 'predicted_label',
       'predicted_sequence']
final_df = pd.concat(final_dfs)
final_df = final_df.rename(columns={'gene_id': 'prediction_id', 'avg_coverage_samples': 'mean_coverage'})
final_df[['CD', 'number']] = final_df['prediction_id'].str.split('_', expand=True)  # sort by species and pred id
final_df['number'] = final_df['number'].astype(int)
final_df = final_df.sort_values(['species', 'number'])
final_df = final_df[final_cols].reset_index(drop=True)
final_df.to_csv(snakemake.output.final_table, sep='\t', index=False)