#!/usr/bin/python3
import pandas as pd
import collections as coll
import re
import regex
import numpy as np
from scipy.signal import find_peaks, savgol_filter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import utils as ut

species_short = snakemake.params.species_short_name
species_short = {v:k for k,v in species_short.items()}
species_colors = snakemake.params.species_colors
species_colors = {species_short[k]:v for k,v in species_colors.items()}
biotype_colors = {'expressed_CD_snoRNA': '#a8ddb5', 'snoRNA_pseudogene': '#88419d'}

len_c_box, len_d_box = 7, 4
shap_df = pd.read_csv(snakemake.input.shap_df, sep='\t').rename(columns={'probability': 'probability_CD'})
FN = shap_df[shap_df['predicted_label'] == 'Other']
# We will deal will TP first using the shap values to determine the snoRNA C and D boxes
# The FN (false negatives) will be dealt using the annotated sequence of the snoRNA to not 
# miss these examples in our second model tune/train/test sets
shap_df = shap_df[shap_df['predicted_label'] != 'Other']
cd_df = pd.read_csv(snakemake.input.cd_df, sep='\t')
fixed_length = int(snakemake.input.cd_df.split('_')[-1].split('nt')[0])
half_window = int((fixed_length + 2) / 2)  # the +2 is to add for the [CLS] and 
                                            # [SEP] tokens at the start and end of the sequence
shap_df = shap_df.merge(cd_df[['gene_id', 'sequence', 
                        f'extended_{fixed_length}nt_sequence', 
                        'species_name', 'gene_biotype']], how='left', on='gene_id')
FN = FN.merge(cd_df[['gene_id', 'sequence', 
                        f'extended_{fixed_length}nt_sequence', 
                        'species_name', 'gene_biotype', 'chr', 'strand', 'start', 'end']], 
                        how='left', on='gene_id')
shap_cols = [i for i in shap_df.columns if i.startswith('SHAP_')]


# The smallest reliably annotated C/D snoRNA is 50 nt long (NR_145814, a C/D pseudogene in human)
# To not miss any snoRNA, we define the minimal length to find the C or D box to be +-15 nt from 
# the center of the predicted window
min_box_dist = 15
# Flanking nt extending after the snoRNA start/end are minimally of 15 nt, so no C or D box 
# should be found in the first and last 15 nt of the window
flanking_nt = 20

# Thus, the ranges in which a C and D boxes can be searched are defined below
## + 1 to account for the [CLS] ##
C_range = range(flanking_nt + 1, half_window-min_box_dist)
D_range = range(half_window+min_box_dist, fixed_length - flanking_nt + 1)


# Deal firstly with the TP examples using their SHAP values
# Based on SHAP values, find the C and D box
# Then find C' and D' boxes and the overall box score
TP_df = ut.find_all_boxes(shap_df, fixed_length, shap_cols, C_range,
                            D_range, flanking_nt)

# Get actual predicted sequence of snoRNAs
cols = ['C_MOTIF', 'C_START', 'C_END', 'D_MOTIF', 'D_START', 'D_END', 
        'C_PRIME_MOTIF', 'C_PRIME_START', 'C_PRIME_END', 'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END',
        'gene_id', 'gene_biotype', 'score_c', 'score_d', 'score_c_prime', 'score_d_prime', 'box_score', 
        'species_name', 'box_provenance']
TP_df['box_provenance'] = 'SHAP_defined'
TP_df = TP_df.merge(cd_df[['gene_id', 'species_name', 'gene_biotype']], how='left', on='gene_id')
TP_df = TP_df[cols]

# Deal with the FN examples to find their C and D box scores based on their sequence only
def find_shift(row):
    return row[f'extended_{fixed_length}nt_sequence'].find(row['sequence'])

FN['left_extension'] = FN.apply(find_shift, axis=1)
fn_examples = []
for i in range(len(FN)):
    sequence = FN.loc[i, 'sequence']
    gene_id = FN.loc[i, 'gene_id']
    gene_bio = FN.loc[i, 'gene_biotype']
    species_name = FN.loc[i, 'species_name']
    left_ext = FN.loc[i, 'left_extension']
    # Find C and D boxes in the 25 first/last nt respectively
    c_motif_FN, c_start_FN, c_end_FN = ut.find_c_box(sequence[0:25], 0, 0)
    d_motif_FN, d_start_FN, d_end_FN = ut.find_d_box(sequence[-25:], len(sequence) - 25, 0)
    
    
    # Find box score
    score_c_FN = ut.hamming(c_motif_FN)
    score_d_FN = ut.hamming(d_motif_FN, consensus='CTGA')

    # Only 2 examples where only the C motif is not found
    # Change to default position of C box starting after the fifth nt
    # and D box ending before the 5th-to-last nt
    if c_motif_FN == 'NNNNNNN':
        c_start_FN = 6
        c_end_FN = c_start_FN + len_c_box
    if d_motif_FN == 'NNNN':
        d_start_FN = len(sequence) - 5
        d_end_FN = d_start_FN + len_d_box
    
    # Find C' and D' boxes
    prime_cols = ['C_PRIME_MOTIF', 'C_PRIME_START', 'C_PRIME_END', 
                'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END']
    c_prime_motif_FN, c_prime_start_FN, c_prime_end_FN, d_prime_motif_FN, d_prime_start_FN, d_prime_end_FN = (
        ut.find_c_prime_d_prime_hamming(sequence, c_start_FN, d_start_FN)
    )
    score_c_prime_FN = ut.hamming(c_prime_motif_FN)
    score_d_prime_FN = ut.hamming(d_prime_motif_FN, consensus='CTGA')

    # Fix box coordinates so that if fits w/r to the fixed_length window
    c_start_FN = c_start_FN + left_ext
    c_end_FN = c_end_FN + left_ext
    d_start_FN = d_start_FN + left_ext
    d_end_FN = d_end_FN + left_ext
    c_prime_start_FN = c_prime_start_FN + left_ext
    c_prime_end_FN = c_prime_end_FN + left_ext
    d_prime_start_FN = d_prime_start_FN + left_ext
    d_prime_end_FN = d_prime_end_FN + left_ext

    
    # Compute box score
    box_score_FN = score_c_FN + score_d_FN + score_c_prime_FN + score_d_prime_FN
    fn_examples.append([c_motif_FN, c_start_FN, c_end_FN, d_motif_FN, d_start_FN, d_end_FN,
                c_prime_motif_FN, c_prime_start_FN, c_prime_end_FN, d_prime_motif_FN, d_prime_start_FN, d_prime_end_FN,
                gene_id, gene_bio, score_c_FN, score_d_FN, score_c_prime_FN, score_d_prime_FN, box_score_FN, 
                species_name, 'annotation_defined'])


# Create df out of all collected box locations
FN_df = pd.DataFrame(fn_examples, columns=cols)
cons_df = pd.concat([TP_df, FN_df])


# Get genomic location of snoRNA based on C and D positions
cons_df = cons_df.merge(cd_df[['gene_id', 'chr', 
                'start', 'end', 'strand',f'extended_{fixed_length}nt_sequence', 'sequence']], 
                how='left', on='gene_id')
cons_df = cons_df.rename(columns={'chr': 'chr_window', 'start':'start_old_seq', 
                    'end':'end_old_seq', 'strand':'strand_window'})
# get location of window around snoRNA using the snoRNA annotated sequence with position
cons_df['shift'] = cons_df.apply(find_shift, axis=1)
cons_df['start_window'] = cons_df['start_old_seq'] - cons_df['shift']  
cons_df['end_window'] = cons_df['start_window'] + fixed_length - 1
# get location of predicted snoRNA based on the window coordinates and box positions
cons_df[['chr', 'start', 'end', 'strand']] = cons_df.apply(
                            ut.get_sno_location, axis=1
)

# Get predicted sequence as well as predicted extended sequence
# These sequence will be used to determine the snoRNA structure stability and terminal stem stability
cons_df['predicted_sequence'] = cons_df.apply(ut.get_seq, axis=1)

# find extended seq (15 nt flanking the snoRNA) to compute the terminal stem stability
cons_df['predicted_extended_sequence'] = cons_df.apply(lambda row: ut.get_seq(row, extended=True), axis=1)
# Correct box position based on the predicted sequence, not the extended window
cons_df = cons_df.apply(ut.correct_box_pos, axis=1, 
                motifs=['C', 'D', 'C_PRIME', 'D_PRIME'])


# Find length difference between annotation and what snoBIRD predicts
cons_df['predicted_len'] = cons_df['predicted_sequence'].apply(lambda x: len(x))
cons_df['real_len'] = cons_df['sequence'].apply(lambda x: len(x))
cons_df['len_diff'] = cons_df.real_len - cons_df.predicted_len
print('Length difference between annotation and predicted snoRNA sequence (5th and 95th percentile)')
print(cons_df.len_diff.quantile(0.05, interpolation='nearest'))
print(cons_df.len_diff.quantile(0.95, interpolation='nearest'))

# Calculate overlap between the annotation of the snoRNA ends 
# and what snoBIRD predicts as the start/end
def calculate_overlap(start1, end1, start2, end2):
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap_length = max(0, overlap_end - overlap_start)
    min_len_element = min(end1 - start1, end2 - start2)
    overlap_proportion = overlap_length / min_len_element
    return overlap_proportion

cons_df['overlap_prop'] = cons_df.apply(lambda row: calculate_overlap(row['start_old_seq'], 
                                row['end_old_seq'], row['start'], row['end']), axis=1)


# Create density plot of the box score of all species or for human, mouse and 
# drosophila separated (only species with snoRNA_pseudogenes)
plt.rcParams['svg.fonttype'] = 'none'
fig, axes = plt.subplots(2, 2, figsize=(25, 25))
plt.subplots_adjust(hspace=0.5)
ax = axes.flatten()
ax_title = ['All species', 'H. sapiens', 'M. musculus', 'D. melanogaster']
dfs = [cons_df, cons_df[cons_df['species_name'] == 'H_sapiens'], 
        cons_df[cons_df['species_name'] == 'M_musculus'], 
        cons_df[cons_df['species_name'] == 'D_melanogaster']]
for i, df in enumerate(dfs):
    sns.kdeplot(x='box_score', hue='gene_biotype', data=df, ax=ax[i], palette=biotype_colors, fill=True)
    ax[i].set_title(ax_title[i], fontdict={'fontsize': 25}, x=0.5, y=1)
plt.savefig(snakemake.output.density_box_score)

# Create jointplot of length difference and overlap proportion between the annotation of snoRNA and the one predicted by SnoBIRD
# So remove the sno for which we found the boxes based on the annotated sequences (i.e. FN examples)
sns.set_style('white')
sns.jointplot(x='len_diff', y='overlap_prop', data=cons_df[cons_df['box_provenance'] != 'annotation_defined'], 
                cmap='Blues', fill=True, kind='kde')
plt.savefig(snakemake.output.jointplot)

sns.jointplot(x='len_diff', y='overlap_prop', data=cons_df[cons_df['box_provenance'] != 'annotation_defined'], 
                hue='species_name', fill=True, kind='kde', palette=species_colors)
plt.savefig(snakemake.output.jointplot_species)

sns.jointplot(x='len_diff', y='overlap_prop', data=cons_df[cons_df['box_provenance'] != 'annotation_defined'], 
                hue='gene_biotype', fill=True, kind='kde', palette=biotype_colors)
plt.savefig(snakemake.output.jointplot_biotype)

print('Overlap proportion between annotation and predicted snoRNA sequence (5th and 95th percentile)')
print(cons_df.overlap_prop.quantile(0.05, interpolation='nearest'))
print(cons_df.overlap_prop.quantile(0.95, interpolation='nearest'))

# Save df
cons_df.drop(columns=['len_diff', 'predicted_len', 'real_len', 'overlap_prop',
            'chr_window', 'strand_window', 'shift']).to_csv(snakemake.output.df, sep='\t', index=False)

