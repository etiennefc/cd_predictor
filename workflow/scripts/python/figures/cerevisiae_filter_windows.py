#!/usr/bin/python3
import pandas as pd
import os
import collections as coll
import seaborn as sns
import matplotlib.pyplot as plt
import functions as ft
from pybedtools import BedTool
import subprocess as sp
from Bio import SeqIO

pred_dir = snakemake.input.pred_dir
#fasta = os.listdir(snakemake.input.fa)
#fasta_dir = snakemake.input.fa
all_preds = os.listdir(pred_dir)
output_df = snakemake.output.df
block_length_output = snakemake.output.density_block_length

# Create bed of C/D snoRNAs in yeast
all_cd = pd.read_csv(snakemake.input.all_cd, sep='\t')
cd_yeast = all_cd[all_cd['species_name'] == 'S_cerevisiae']
#cd_yeast = all_cd[all_cd['species_name'] == 'D_melanogaster']
cd_yeast['score'] = '.'
#cd_yeast = cd_yeast[cd_yeast['chr'] == '2R']
cd_yeast[['chr', 'start', 'end', 'gene_id', 'score', 'strand']].to_csv('temp_cd_yeast.bed', sep='\t', 
                                                                        index=False, header=False)
print(cd_yeast[['chr', 'start', 'end', 'strand']].sort_values(by=['start', 'strand']))     
print(len(cd_yeast))                                                    
cd_yeast_bed = BedTool('temp_cd_yeast.bed')


# Concat all predictions into 1 df for S. cerevisiae and other small genomes
dfs = []
for p in all_preds:
    df = pd.read_csv(f'{pred_dir}/{p}', sep='\t')
    dfs.append(df)
df_preds = pd.concat(dfs)
print(df_preds)
''''''


"""For Droso only (big chr)"""
'''
droso_sklias = BedTool(snakemake.input.cd_bed_sklias)
split_sizes = {}
for f in fasta:
    rec = list(SeqIO.parse(f'{fasta_dir}/{f}', 'fasta'))[0]
    chr_split_name = f.split('.fa')[0]
    split_sizes[chr_split_name] = len(rec.seq)

def get_real_size_chr(chr_name, num, size_dict):
    # Get the real start and end for a given chr that was split
    # by adding the size of previous splits to the one we're interested in
    split_num_before = [i for i in range(1, num)]
    chr_original = chr_name.rsplit('_')[0]
    sum_size = 0
    for n in split_num_before:
        sum_size += size_dict[f'{chr_original}_{n}']
    return sum_size

dfs = []
for p in all_preds:
    df = pd.read_csv(f'{pred_dir}/{p}', sep='\t')
    chr_split = p.split('windows_')[1].split('.tsv')[0]
    num_split = int(chr_split.split('_')[1])
    actual_size = get_real_size_chr(chr_split, num_split, split_sizes)
    df['start'] = df['start'] + actual_size
    df['end'] = df['end'] + actual_size
    dfs.append(df)
df_preds = pd.concat(dfs)
df_preds['chr'] = df_preds['chr'].replace(r'_[0-9]*[0-9*]', '', regex=True)
df_preds = df_preds.sort_values(by=['chr', 'strand', 'start'])
df_preds.to_csv('pos_dros', sep='\t', index=False, header=False)
print(df_preds)
'''


# Filter predictions based on probability
df_preds = df_preds[df_preds['probs'] > 0.999].sort_values(by=['chr', 'start', 'strand'])
print(df_preds)

# Calculate difference in start positions between consecutive rows
# It creates a NA for each first row of a group in the groupby
df_preds['diff'] = df_preds.groupby(['chr', 'strand'])['start'].diff()



# Create boolean mask to highlight where diff >1 (i.e. a new stretch of consecutive nt is created)
# Use cumulative sum to assign a different number to each stretch of consecutive nts
mask = (df_preds['start'].shift(-1) == df_preds['start'] + 1)
df_preds.loc[(df_preds['diff'].isna()) & (mask), 'diff'] = 2  # manually change these NA to 2 because they should be counted in the block length (first element of the block)
df_preds.loc[df_preds['diff'].isna(), 'diff'] = 2  # manually change the remaining NA that are not part of stretches so that they are assigned a different stretch_id
df_preds['stretch_id'] = (df_preds['diff'] > 1).cumsum()
df_preds['stretch_id'] =  'block_' + df_preds['stretch_id'].astype(str)
df_preds = df_preds[['chr', 'start', 'end', 'stretch_id', 'strand', 'probs']]
df_preds.to_csv('temp_preds.bed', sep='\t', index=False, header=False)
preds_bed = BedTool('temp_preds.bed')

# Create density plot of the block length of positive predictions
vals = dict(coll.Counter(df_preds.stretch_id)).values()  # get the length of block/stretch per stretch id
################################################################ft.density(vals, 'Length of positive predictions block (nt)', 'Density', 'Filtered prob > 0.999', block_length_output)

print(df_preds)
print(len(vals))




# Merge consecutive positive preds into one entry and create bed out of all blocks
# this is done via groupby over the stretch_id column and keeping the lowest/highest start/end of all windows in that block
# and by averaging the probs over all windows in that block
merged_blocks = preds_bed.groupby(g=[4], c=[2, 3, 4, 6, 5, 1], o=['min', 'max', 'first', 'mean', 'first', 'first'])
merged_blocks = merged_blocks.cut([6, 1, 2, 3, 4, 5]).to_dataframe()  # reorder columns of bed in right order

# Correct bed coordinates (since they are 1-based)
merged_blocks['start'] = merged_blocks['start'] + 1
print('\n\n')
print(merged_blocks.sort_values(by='start'))
print('\n\n')
merged_blocks = merged_blocks.sort_values(by=['chrom', 'strand', 'start', 'end'])

merged_blocks.to_csv('temp_test_merged.bed', sep='\t', index=False, header=False)



# Merge blocks that overlap reciprocally to at least 50%
merged_blocks['len'] = merged_blocks.end - merged_blocks.start + 1


'''Test filtering for len at that stage before merging overlapping blocks not consecutive??'''
merged_blocks = merged_blocks[merged_blocks['len'] > 217]
print(merged_blocks)
print('STARTTTTTTTTTTTT\n\n\n\n\n\n')

def merge_rows(row1, row2):
    # Function to merge two rows if they overlap reciprocally by at least 50%
    overlap_start = max(row1['start'], row2['start'])
    overlap_end = min(row1['end'], row2['end'])
    overlap_len = max(0, overlap_end - overlap_start)
    max_len = max(row1['len'], row2['len'])
    
    reciprocal_overlap = (overlap_len / max_len) >= 0.5  # test only for the longest block (the shortest will be also true 
                                                         # if the longest overlaps at least at 50%)
    
    if (row1['chrom'] == row2['chrom']) & (row1['strand'] == row2['strand']) & (reciprocal_overlap == True):
        merged_row = {
            'chrom': row1['chrom'],
            'start': row1['start'],
            'end': row2['end'],
            'name': str(row1['name']) + '_' + str(row2['name']),
            'score': (row1.score + row2.score) / 2,
            'strand': row1['strand'], 
            'len': row2['end'] - row1['start'] + 1
        }
        #print(row1)
        #print(row2)
        #print(merged_row)
        return merged_row
    else:
        return None

# Apply the merging function row-wise
merged_rows = []
current_row = None

for i in range(len(merged_blocks) - 1):
    if current_row is None:
        current_row = merged_blocks.iloc[i].copy()
    next_row = merged_blocks.iloc[i + 1].copy()
    
    merged_row = merge_rows(current_row, next_row)
    
    if merged_row is not None:
        current_row = pd.Series(merged_row)
        #merged_rows.append(current_row)
    else:
        merged_rows.append(current_row)
        current_row = None

## Add the last row if it was not merged
#if current_row is not None:
 #   merged_rows.append(current_row)

''' DEAL with last block if merge or not to return!!!!!!!'''

# Create a DataFrame from the merged rows
result_df = pd.DataFrame(merged_rows)

# Display the result
print(result_df.reset_index(drop=True)[['chrom', 'start', 'end', 'score', 'len', 'name']])
#result_df = result_df[result_df['score'] > 0.9]
#result_df = result_df[(result_df['len'] >= 218) & (result_df['len'] <= 269)]
result_df = result_df[(result_df['len'] >= 218) & (result_df['len'] <= 269) & (result_df['score'] > 0.9997)]
print(result_df.reset_index(drop=True)[['chrom', 'start', 'end', 'score', 'len', 'name']])
result_df.reset_index(drop=True)[['chrom', 'start', 'end', 'name', 'score', 'strand', 'len']].to_csv('test_test.bed', sep='\t', index=False, header=False)
#print(coll.Counter(result_df.len))
#print(result_df[(result_df['len'] >= 205) & (result_df['len'] <= 298)])
#print(result_df.reset_index(drop=True)[['chrom', 'start', 'end', 'score', 'len', 'name']].sort_values(by='len'))


merged_blocks_bed = BedTool('test_test.bed')






#a = BedTool('temp_test_merged.bed')

# If block overlaps with annotated C/D, return as positive
# enforce strandedness and make sure 100% (F=1) of the snoRNA is included in the block
pos = merged_blocks_bed.intersect(cd_yeast_bed, u=True, s=True)  
#pos = merged_blocks_bed.intersect(droso_sklias, u=True, s=True)
#pos = a.intersect(cd_yeast_bed, u=True, s=True)
print(pos.to_dataframe().sort_values(by='thickStart'))
pos = cd_yeast_bed.intersect(merged_blocks_bed, u=True, s=True)
#pos = droso_sklias.intersect(merged_blocks_bed, u=True, s=True, f=1)
pos = pos.to_dataframe()
#print(len(cd_yeast))
#print(cd_yeast[['chr', 'strand', 'start', 'end']].sort_values(by='start'))
print(pos.sort_values(by='start'))
print(len(pos))
#print(pos.sort_values(by='thickStart').reset_index(drop=True))
pos.to_csv('test_pos.tsv', sep='\t', index=False, header=False)

"""
'''

ft.density_x([result_df['len'], pos['thickStart']], 'Length of all merged blocks', 'Density', 'linear', 
            '', ['lightblue', 'red'], ['All TP/FP', 'Only TP'], 'length_blocks.svg')

ft.density_x([result_df['score'], pos['score']], 'Score of all merged blocks', 'Density', 'log', 
            '', ['lightblue', 'red'], ['All TP/FP', 'Only TP'], 'score_blocks.svg')


result_df.loc[result_df['name'].isin(list(pos['name'])), 'pos_type'] = 'TP'
result_df['pos_type'] = result_df['pos_type'].fillna('FP')
result_df = result_df.reset_index(drop=True)

# Filter result_df based on length
#result_df = result_df[(result_df['len'] >= 209) & (result_df['len'] <= 306)]
#result_df = result_df[result_df['score'] > 0.9994]

grid = sns.JointGrid(x='len', y='score', data=result_df)

g = grid.plot_joint(sns.scatterplot, data=result_df, hue='pos_type')
sns.kdeplot(result_df.loc[result_df['pos_type']=='TP', 'len'], ax=g.ax_marg_x, color='tab:orange', fill=True, legend=False)
sns.kdeplot(result_df.loc[result_df['pos_type']=='FP', 'len'], ax=g.ax_marg_x, color='tab:blue', fill=True, legend=False)
sns.kdeplot(result_df.loc[result_df['pos_type']=='TP', 'score'], ax=g.ax_marg_y, vertical=True, color='tab:orange', fill=True, legend=False)
sns.kdeplot(result_df.loc[result_df['pos_type']=='FP', 'score'], ax=g.ax_marg_y, vertical=True, color='tab:blue', fill=True, legend=False)
g.ax_marg_y.set(yscale='log')
#from matplotlib.ticker import ScalarFormatter
#formatter = ScalarFormatter(useMathText=False)
#g.ax_marg_y.yaxis.set_major_formatter(formatter)

#g.ax_marg_y.tick_params(axis='y', labelsize=2)
#penguins = sns.load_dataset("penguins")
#print(result_df)
#print(result_df[['pos_type', 'len', 'score']])
#print(penguins)
#print(penguins.columns)
#sns.jointplot(data=result_df, x="len", y="score", hue="pos_type")
plt.savefig('joint_plot_blocks_after_filter.svg', dpi=600, bbox_inches='tight')




#sp.call('rm temp_merged_blocks.bed temp_cd_yeast.bed temp_preds.bed', shell=True)
# Find if block overlaps with annotated C/D snoRNA


#for i in range (1, max(vals) + 1):
 #   print(i)



#filter consecutive blocks of 8 nt (based on snr69) and merge afterwards if overlap is >50%.

'''
"""


