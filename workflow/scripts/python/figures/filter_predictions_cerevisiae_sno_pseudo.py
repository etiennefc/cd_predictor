#!/usr/bin/python3
import pandas as pd
import subprocess as sp
from pybedtools import BedTool
import collections as coll

# Load C/D bed
all_cd = pd.read_csv(snakemake.input.all_cd, sep='\t')
cerevisiae_cd = all_cd[all_cd['species_name'] == 'S_cerevisiae']
cerevisiae_cd['score'] = '.'
cerevisiae_cd = cerevisiae_cd[['chr', 'start', 'end', 'gene_id', 'score', 'strand']]
cerevisiae_cd.to_csv('temp_cerevisiae_cd.bed', sep='\t', index=False, header=False)
print(cerevisiae_cd.sort_values(by=['chr', 'start', 'strand']))
cd_bed = BedTool('temp_cerevisiae_cd.bed')

# Load preds and their probability as bed
preds = pd.read_csv(snakemake.input.sno_pseudo_windows, sep='\t')
probs = pd.read_csv(snakemake.input.sno_pseudo_preds, sep='\t')

df = preds.merge(probs, how='left', on='gene_id')
df['block_id'] = df['gene_id'].str.rsplit('_', n=1).apply(lambda x: x[0])

# Select only positive predictions
df = df[df['predicted_label'] == 1]
df = df[['chrom', 'start', 'end', 'gene_id', 'probability', 'strand', 'block_id']]


# Add filtering step (score/consecutive windows)
#df = df[df['probability'] > 0.9]


print(len(pd.unique(df.block_id)))
#pred_bed = BedTool('preds_bed.bed')

## Merge consecutive stretches of positive predictions and filter on size
df['diff'] = df.groupby(['chrom', 'strand'])['start'].diff()
mask = (df['start'].shift(-1) == df['start'] + 1)
df.loc[(df['diff'].isna()) & (mask), 'diff'] = 2  # manually change these NA to 2 because they should be counted in the block length (first element of the block)
df.loc[df['diff'].isna(), 'diff'] = 2  # manually change the remaining NA that are not part of stretches so that they are assigned a different stretch_id
df['stretch_id'] = (df['diff'] > 1).cumsum()
df['stretch_id'] =  'new_block_' + df['stretch_id'].astype(str)

df['stretch_count'] = df.groupby('stretch_id')['stretch_id'].transform('count')
#df.to_csv('stretch_test.tsv', sep='\t', index=False)
##print(coll.Counter(df.stretch_id))


#Filter on stretch count (number of initial consecutive windows)
df = df[df['stretch_count'] > 11]
print(len(pd.unique(df.block_id)))
print(df.columns)


#df = df[['chrom', 'start', 'end', 'stretch_id', 'probability', 'strand', 'stretch_count']]
df.to_csv('preds_bed.bed', index=False, sep='\t', header=False)
#print(df)
pred_bed = BedTool('preds_bed.bed')



# Merge entries in same stretch
merged_blocks = pred_bed.groupby(g=[9], c=[2, 3, 4, 5, 6, 1, 9], o=['min', 'max', 'first', 'mean', 'first', 'first', 'first'])
#print(merged_blocks.to_dataframe())
merged_blocks = merged_blocks.cut([6, 1, 2, 0, 4, 5, 3]).to_dataframe() # reorder columns of bed in right order
merged_blocks['len'] = merged_blocks.end - merged_blocks.start + 1
print(merged_blocks)
print(merged_blocks['score'])

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

""" DEAL with last block if merge or not to return!!!!!!!"""

# Create a DataFrame from the merged rows
result_df = pd.DataFrame(merged_rows)
result_df = result_df[result_df['len'] > 201]
result_df = result_df[result_df['score'] > 0.92]
print(len(result_df))
result_df[['chrom', 'start', 'end', 'name', 'score', 'strand', 'len']].to_csv('oktest.bed', sep='\t', index=False, header=False)

merged_consec = BedTool('oktest.bed')



# Intersect predicted windows with 
# enforce strandedness and make sure at least 25% (F=1) of the snoRNA is included in the block
pos = merged_consec.intersect(cd_bed, u=True, s=True, F=0.25)  

pos = pos.to_dataframe()
pos.to_csv('pos_expcd.tsv', sep='\t', index=False)

#print(pos)
print(len(pd.unique(pos['name'])))

print(pos.sort_values(by=['thickStart'], ascending=False).reset_index(drop=True))

# See if a sno is predicted as pseudo/expressed on all the block or by parts
# Filter by that on the % of positive across all the windows in initial pos block






#sp.call('rm temp_cerevisiae_cd.bed preds_bed.bed', shell=True)