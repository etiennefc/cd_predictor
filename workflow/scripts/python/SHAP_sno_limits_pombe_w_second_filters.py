#!/usr/bin/python3
import sys
import subprocess as sp
import pandas as pd
import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from math import ceil
import utils as ut

# Define model path and other variables
shap_pombe = pd.read_csv(sys.argv[1], sep='\t')
genome = sys.argv[2]
fixed_length = int(sys.argv[6])
pred_seqs = pd.read_csv(sys.argv[3], sep='\t', names=['chr_window', 
                        'start_window', 'end_window', 'gene_id', 
                        'prob_first_model', 'strand_window', 'block_id', 
                        f'extended_{fixed_length}nt_sequence'])
second_model_preds = pd.read_csv(sys.argv[4], sep='\t').rename(
                    columns={'y_pred': 'second_model_prediction'})                        
pretrained_model = sys.argv[5]

half_window = int((fixed_length + 2) / 2)
batch_size = int(sys.argv[7])
output_df = sys.argv[8]
num_labels = 2
len_c_box, len_d_box = 7, 4

# The smallest reliably annotated C/D snoRNA is 50 nt long (NR_145814, a
# C/D pseudogene in human. To not miss any snoRNA, we define the minimal length
# to find the C or D box to be +-15 nt from the center of the predicted window
min_box_dist = 15
# Flanking nt extending after the snoRNA start/end are minimally of 15 nt
# (+ 5 nt to get to the C_start or D_end), so
# no C or D box should be found in the first and last 20 nt of the window
flanking_nt = 20

# Thus, the ranges in which a C and D boxes can be searched are defined below
## + 1 to account for the [CLS] ##
C_range = range(flanking_nt + 1, half_window-min_box_dist)
D_range = range(half_window+min_box_dist, fixed_length - flanking_nt + 1)

# Select only examples predicted as C/D from the first model
# This is an additional filter so that the centered window is also 
# predicted positively
shap_pombe = shap_pombe[shap_pombe['predicted_label'] == 'CD_snoRNA']

# Add extended sequence
shap_pombe = shap_pombe.merge(pred_seqs[
                        ['gene_id', f'extended_{fixed_length}nt_sequence', 
                        'chr_window', 'start_window', 'end_window', 
                        'strand_window']], how='left', on='gene_id')
shap_pombe = shap_pombe.rename(columns={'probability': 'probability_CD'})

shap_cols = [i for i in shap_pombe.columns if i.startswith('SHAP_')]


# Based on SHAP values, find the C and D box
# Then find C' and D' boxes and the overall box score
pombe_sno_df = ut.find_all_boxes(shap_pombe, fixed_length, shap_cols, C_range,
                            D_range, flanking_nt)

# Get snoRNA sequence, location, length and normalized_mfe
pombe_sno_df = pombe_sno_df.merge(pred_seqs[['gene_id', 'chr_window', 
                'start_window', 'end_window', 'strand_window']], 
                how='left', on='gene_id')
pombe_sno_df['predicted_sequence'] =  pombe_sno_df.apply(ut.get_seq, axis=1)
# find extended seq (15 nt flanking the snoRNA) to compute 
# the terminal stem stability
pombe_sno_df['predicted_extended_sequence'] = pombe_sno_df.apply(
                            lambda row: ut.get_seq(row, extended=True), axis=1)
pombe_sno_df[['chr', 'start', 'end', 'strand']] = pombe_sno_df.apply(
                            ut.get_sno_location, axis=1
)
pombe_sno_df['sno_length'] = pombe_sno_df['predicted_sequence'].apply(
                                lambda x: len(x))

pombe_sno_df = ut.sno_mfe(pombe_sno_df, 'predicted_sequence',
                    'sno_length')


# Get terminal stem stability/score and combined
pombe_sno_df = ut.terminal_stem(pombe_sno_df, 
                f'extended_{fixed_length}nt_sequence', extended_col=True, 
                extended_col_name='predicted_extended_sequence')

# Merge with second model predictions
pombe_sno_df = pombe_sno_df.merge(second_model_preds[
                    ['gene_id', 'block_id', 'second_model_probability', 
                    'second_model_prediction']], how='left', on='gene_id')

pombe_sno_df['normalized_sno_stability'] = pombe_sno_df[
                                        'normalized_sno_stability'].round(2)
pombe_sno_df['terminal_combined'] = pombe_sno_df['terminal_combined'].round(2)


# Filter predictions using the feature values
pombe_sno_df = ut.feature_filters(pombe_sno_df, 'second_model_probability', 
            'second_model_prediction', 'second_model_filtered_prediction')


# Correct the start/end of boxes based on the snoRNA seq and not 
# the extended window sequence
boxes = ['C', 'D', 'C_PRIME', 'D_PRIME']
pombe_sno_df = pombe_sno_df.apply(ut.correct_box_pos, axis=1, motifs=boxes)

# Save final df
final_cols = ['gene_id', 'block_id', 'chr', 'start', 'end', 'strand', 
            'probability_CD', 'C_MOTIF', 'C_START', 'C_END',
            'D_MOTIF', 'D_START', 'D_END', 'C_PRIME_MOTIF', 'C_PRIME_START', 
            'C_PRIME_END', 'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END', 
            'box_score', 'terminal_combined', 'normalized_sno_stability', 
            'second_model_probability', 'second_model_prediction', 
            'second_model_filtered_prediction', 'sno_length', 
            'predicted_sequence']

#print(pombe_sno_df[final_cols])
import collections as coll
print(coll.Counter(pombe_sno_df['second_model_filtered_prediction']))
print(pombe_sno_df[pombe_sno_df['gene_id'].isin(['CD_442', 'CD_720', 'CD_113', 'CD_367', 'CD_130', 'CD_387', 'CD_644', 'CD_704', 'CD_705', 'CD_865'])]['second_model_filtered_prediction'])
pombe_sno_df[final_cols].to_csv(output_df, sep='\t', index=False)

