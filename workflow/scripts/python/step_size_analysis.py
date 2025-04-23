#!/usr/bin/python3
import pandas as pd
from glob import glob
import os
#import collections as coll
import utils as ut
from pybedtools import BedTool
import subprocess as sp
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# Define inputs, outputs and parameters
all_preds = glob(snakemake.input.preds_s1+'/positive_windows*')
input_fasta_dir = snakemake.input.input_fasta_dir
fixed_length = snakemake.params.fixed_length
prob_threshold = snakemake.params.prob_threshold
min_consecutive_windows = snakemake.params.min_consecutive_windows
step_size = [i for i in range(1,51)]    
species = snakemake.wildcards.species   
plot_output = snakemake.output.step_size_plot
output_df = snakemake.output.df

# Create bed of known annotated C/D that fit within the range of detection 
# of snoBIRD
cd_df = pd.read_csv(snakemake.input.known_cd, sep='\t', dtype={'chr': 'str'})
cd_df = cd_df[cd_df['sno_type'] == 'C/D']
if species == 'homo_sapiens':
    cd_df['start'] = cd_df['start'].astype(int)
    cd_df['end'] = cd_df['end'].astype(int)
    cd_df = cd_df.dropna(subset=['gene_id'])
cd_df = cd_df[(cd_df['end'] - cd_df['start']) <= fixed_length - 30]
cd_df['dot'] = '.'
cd_df[['chr', 'start', 'end', 'gene_id', 
            'dot', 'strand', 'sno_type']].to_csv(f'{species}_known_cd.bed', 
            sep='\t', index=False, header=False)

bed_cd = BedTool(f'{species}_known_cd.bed')



# Concat all predictions of first SnoBIRD model into 1 df (preds from all 
# chr and/or chunks)
""" Apply 1st filter on probability."""
dfs = []
for p in all_preds:
    df = pd.read_csv(p, sep='\t')
    df = df.rename(columns={'probs': 'probability_first_model'})
    df = df[df['probability_first_model'] > prob_threshold]
    dfs.append(df)
df_preds = pd.concat(dfs)

# For chunks of chr, rectify genomic location of positive windows based on 
# the number of previous chunks
# Preds on - strand have already been corrected during genome_prediction 
chunk_fa = [path for path in os.listdir(
            input_fasta_dir+'/fasta/') if '_chunk_' in path]

rectify_dict = {}
if len(chunk_fa) > 0:
    for chunk in chunk_fa:
        with open(input_fasta_dir+'/fasta/'+chunk, 'r') as f:
            first_line = f.readline().replace('>', '').replace('\n', '')
            chunk_name = first_line.split('  ')[0]
            chunk_number = int(chunk_name.split('_')[-1])
            prev_size = int(first_line.split(' ')[2].replace(
                        'prev_size=', ''))
            total_size = int(first_line.split(' ')[3].replace(
                            'total_size=', ''))
            rectify_dict[chunk_name] = (prev_size, total_size)

def rectify_pos(row, chr_dict):
    # Rectify position based on the number of nt before a given chunk
    if row['chr'] in chr_dict.keys():
        nt_number = chr_dict[row['chr']][0]
        row['start'] = row['start'] + nt_number 
        row['end'] = row['end'] + nt_number
        row['chr'] = row['chr'].split('_chunk_')[0]
    return row
df_preds = df_preds.apply(rectify_pos, axis=1, chr_dict=rectify_dict)
df_preds['chr'] = df_preds['chr'].astype(str)
df_preds = df_preds.sort_values(by=['chr', 'start', 'strand'])

# To fix memory issues (usually when s=1 and df is very large!)
df_preds.to_csv('df_preds_step1.tsv', sep='\t', index=False)
df_preds = pd.read_csv('df_preds_step1.tsv', sep='\t', dtype={'chr': 'str'})
df_preds_step1 = df_preds.copy()

# Depending on the synthetic step size s, remove windows that would not have 
# been predicted if s>1
total_, annotated_ = [], []
for step in step_size:
    if step > 1:
        df_preds = df_preds_step1[df_preds_step1['start'] % step == 1]
    
    # Calculate difference in start positions between consecutive rows
    # It creates a NA for each first row of a group in the groupby
    # This is to further filter on the nb of consecutive windows
    df_preds['diff'] = df_preds.groupby(['chr', 'strand'])['start'].diff()


    # Create boolean mask to highlight where diff >1 (
    # i.e. a new stretch of consecutive nt is created)
    # Use cumulative sum to assign a different number to each stretch 
    # of consecutive nts
    mask = (df_preds['start'].shift(-1) == df_preds['start'] + 1)

    # Manually change these NA to 2 because they should be counted in the 
    # block length (first element of the block)
    df_preds.loc[(df_preds['diff'].isna()) & (mask), 'diff'] = 2  

    # Manually change the remaining NA that are not part of stretches so 
    # that they are assigned a different stretch_id
    df_preds.loc[df_preds['diff'].isna(), 'diff'] = 2  
    df_preds['stretch_id'] = (df_preds['diff'] > 1).cumsum()
    df_preds['stretch_id'] =  'block_' + df_preds['stretch_id'].astype(str)
    df_preds = df_preds[['chr', 'start', 'end', 'stretch_id', 'strand', 
                        'probability_first_model']]
    df_preds.to_csv(f'temp_preds_{step}_{species}.bed', sep='\t', index=False, 
                    header=False)
    preds_bed = BedTool(f'temp_preds_{step}_{species}.bed')

    # Merge consecutive positive preds into one entry and create bed out of 
    # all blocks this is done via groupby over the stretch_id column & keeping 
    # the lowest/highest start/end of all windows in that block
    # and by averaging the probabilities over all windows in that block
    merged_blocks = preds_bed.groupby(g=[4], c=[2, 3, 4, 6, 5, 1], 
                    o=['min', 'max', 'first', 'mean', 'first', 'first'])
    # reorder columns of bed in right order
    merged_blocks = merged_blocks.cut([6, 1, 2, 3, 4, 5]).to_dataframe()
    merged_blocks['chrom'] = merged_blocks['chrom'].astype(str)  
    merged_blocks = merged_blocks.sort_values(by =
                                    ['chrom', 'strand', 'start', 'end'])

    # Get length of block
    merged_blocks['len'] = merged_blocks.end - merged_blocks.start + 1

    # Merge blocks that overlap reciprocally to at least 50%
    # Apply the merging function row-wise
    merged_rows = []
    current_row = None

    for i in range(len(merged_blocks)):
        if i < len(merged_blocks) - 1:
            if current_row is None:  # first row | after previous block is over
                current_row = merged_blocks.iloc[i].copy()
            next_row = merged_blocks.iloc[i + 1].copy()
            merged_row = ut.merge_rows(current_row, next_row)

            # Current_row (can be a block) is merged with next row/block
            if merged_row is not None:  
                current_row = pd.Series(merged_row)
            else:  # no merge
                merged_rows.append(current_row)
                current_row = None
        else:  # deal with last row if it's not included in last block
            last_row = merged_blocks.iloc[i]
            # Last block has ended and doesn't include last row
            if current_row is None:  
                merged_rows.append(last_row)
            # otherwise it has already been merged with before last row/block



    # Create a DataFrame from the merged rows
    result_df = pd.DataFrame(merged_rows)
    

    """ Apply 2nd filter: length of merged blocks"""
    result_df = result_df[result_df['len'] >= (
                                        fixed_length + min_consecutive_windows)]

    # Deal with large blocks that can contain more than 1 snoRNA (ex: clustered 
    # snoRNAs). They might not have a centered positively predicted window as there 
    # are likely two or more in the same block, so not necessarily centered
    # Include also step-size specific thresholds?



    # Identify the center window of fixed_length in the merged block
    center_window = result_df.reset_index(drop=True).copy()
    center_window[['centered_start', 'centered_end']] = center_window.apply(
                        lambda row: ut.centered_window(row, fixed_length), axis=1)
    # correct for 0-based bed
    center_window['centered_start'] = center_window['centered_start'] -1

    df1_windows = set(df_preds_step1.apply(lambda row: (
                    row['chr'], row['start'], row['end'], row['strand']), axis=1))

    
    # Check if rows in df2 are present in df1
    """ Apply 3rd filter: if center window is not predicted as CD in the merged 
        block, don't consider that prediction (if step_size=1, thus this centered 
        window was predicted by SnoBIRD, but not as a CD) or predict on that center 
        window (if step_size>1, thus this centered window was never predicted on by 
        SnoBIRD because the larger step size passed over it)."""
    center_window['is_present_in_df1'] = center_window.apply(lambda row: (
            row['chrom'], row['centered_start'], row['centered_end'], row['strand']
            ) in df1_windows, axis=1)
    
    center_window = center_window[center_window['is_present_in_df1'] == True]
    total_pred = len(center_window)



    # Add prediction id to each center window that are predicted as CD
    center_window = center_window.reset_index(drop=True)
    center_window['index_'] = center_window.index + 1
    center_window['prediction_id'] = 'CD_' + center_window.index_.astype(str)
    cols_bed = ['chrom', 'centered_start', 'centered_end', 'prediction_id', 
                'score', 'strand', 'name']
    center_window[cols_bed].reset_index(drop=True).to_csv(
                            f'center_window_{step}_{species}.bed', 
                            sep='\t', index=False, header=False)
    bed_center_window = BedTool(f'center_window_{step}_{species}.bed')


    # Find overlap between annotated C/D and predictions
    pos_overlap = bed_cd.intersect(bed_center_window, s=True, wb=True, 
                    f=0.5).to_dataframe(names=
                    ['chr', 'start', 'end', 'stretch_id', 'strand_cd', 
                        'probability_first_model']+cols_bed)
    num_pos = len(pos_overlap)

    print(step, total_pred, num_pos)
    total_.append(total_pred)
    annotated_.append(num_pos)



plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
ax.plot(step_size, total_, color='grey', label='Total predicted snoRNAs')
ax.axhline(y=len(cd_df), color='black', ls='--', linewidth=1, alpha=0.5, 
            label='Total annotated snoRNAs')
ax.plot(step_size, annotated_, color='lightblue', 
        label='Predicted annotated snoRNAs')
ax.axvspan(4.5, 5.5, color='lightgrey', alpha=0.5, label='Step size = 5')
fig.suptitle(f"Step size effect on SnoBIRD's\nprediction in {species}", 
            fontsize=25)
ax.set_xlabel("Step size (SnoBIRD's first step)", fontdict={'fontsize': 25})
ax.set_ylabel("Number of predictions", fontdict={'fontsize': 25})
ax.set_yscale('log')
#ax.set_yscale('symlog', linthresh=100)
ax.set_ylim(1, 100+max(total_))
ax.set_xlim(min(step_size), max(step_size))
ax.legend()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig(plot_output, bbox_inches='tight', dpi=500)

final_df = pd.DataFrame([step_size, total_, annotated_]).T
final_df.columns = ['step_size',  'total_nb_preds', 'total_predicted_annotated']
final_df.to_csv(output_df, sep='\t', index=False)

print(final_df)

sp.call(
    f'rm temp_preds_*_{species}.bed center_window_*_{species}.bed {species}_known_cd.bed df_preds_step1.tsv', 
        shell=True)
