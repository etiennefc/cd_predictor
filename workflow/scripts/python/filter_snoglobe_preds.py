#!/usr/bin/python3
import pandas as pd

# Load dfs
sno_df = pd.read_csv(snakemake.input.sno_table, sep='\t')

snoglobe_preds = []
for p in snakemake.input.snoglobe:
    species = p.split('predictions_')[1].replace('.tsv', '')
    df = pd.read_csv(p, sep='\t')
    df['species'] = species

    #Filter to have min 3 consecutive windows and a mean score of >=0.98
    df = df[(df['count'] >= 3) & (df['mean_score'] >= 0.98)]
    snoglobe_preds.append(df)

snoglobe = pd.concat(snoglobe_preds)


# Merge with snoBIRD pred table to filter based on D and D' positions 
snoglobe = snoglobe.merge(sno_df[['prediction_id', 'species', 'D_MOTIF', 
                        'mean_coverage', 'chr', 'start', 'end', 'strand',
                        'probability_CD', 'box_score', 'C_MOTIF', 'C_PRIME_MOTIF',
                        'C_START', 'C_END', 'C_PRIME_START', 'C_PRIME_END',
        'D_START', 'D_END', 'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END', 
        'terminal_stem_score', 'normalized_sno_stability', 
        'probability_expressed_pseudogene', 'predicted_label', 'predicted_sequence']], 
        left_on=['species', 'sno_id'], right_on=['species', 'prediction_id'], 
        how='left')

# Keep only the highest probability interaction per snoRNA
idx = snoglobe.groupby(['species', 'sno_id'])['mean_score'].idxmax()
result_df = snoglobe.loc[idx].reset_index(drop=True)

# Filter these interactions so that they are located either between the 
# C and D' boxes or C' and D boxes (they can overlap by 1 nt with the D and D' box)
result_df = result_df[((result_df['sno_window_start'] > result_df['C_END']) & 
                (result_df['sno_window_end'] <= result_df['D_PRIME_START'])) | 
                ((result_df['sno_window_start'] > result_df['C_PRIME_END']) & 
                (result_df['sno_window_end'] <= result_df['D_START']))]

result_df[['prediction_id', 'species', 'mean_coverage', 'chr', 'start', 'end', 
        'strand', 'probability_CD', 'box_score', 'C_MOTIF', 'C_START', 'C_END', 
        'D_MOTIF', 'D_START', 'D_END', 'C_PRIME_MOTIF', 'C_PRIME_START', 
        'C_PRIME_END', 'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END',
        'terminal_stem_score', 'normalized_sno_stability', 
        'probability_expressed_pseudogene', 'predicted_label', 
        'predicted_sequence', 'target_id', 'target_window_start', 
        'target_window_end', 'mean_score', 'target_seq', 'sno_seq' 
        ]].to_csv(snakemake.output.filtered_preds, sep='\t', index=False)




