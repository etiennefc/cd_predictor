#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
import numpy as np


pred_dir = snakemake.input.pred_dir
cd_cerevisiae = str(snakemake.wildcards.cd_cerevisiae)
output = snakemake.output.lineplot
all_cd = pd.read_csv(snakemake.input.all_cd, sep='\t')

# Keep only the given C/D of interest from S. cerevisiae
filtered_cd = all_cd[all_cd['species_name'] == 'S_cerevisiae']
filtered_cd = filtered_cd[filtered_cd['gene_id'] == cd_cerevisiae]

cd_chr = filtered_cd['chr'].values[0]
cd_strand = filtered_cd['strand'].values[0]
cd_start = filtered_cd['start'].values[0]
cd_end = filtered_cd['end'].values[0]
print(cd_cerevisiae, cd_start, cd_end, cd_strand)

sno_len = cd_end - cd_start + 1
diff = 190 - sno_len
if diff % 2 != 0:  # odd number
    left_ext = int(diff // 2)
    right_ext = int(left_ext + 1)
else:
    left_ext = int(diff / 2)
    right_ext = int(diff / 2)

# Get predictions for the chr of that C/D
preds = pd.read_csv(f'{pred_dir}/pos_windows_{cd_chr}.tsv', sep='\t')

surrounding_preds = preds[(preds['start'] >= cd_start - left_ext - 100) & (preds['end'] <= cd_end + right_ext + 100)]

# Create total df with all possible ranges of window start/end between cd_start/cd_end +/- 100
all_combinations = [[cd_chr] * 201, [cd_strand] * 201, [start for start in range(cd_start - left_ext - 100, cd_start - left_ext + 100 + 1)], 
                    [start + 190 for start in range(cd_start - left_ext - 100, cd_start - left_ext + 100 + 1)], [0]*201]
full_range_df = pd.DataFrame(all_combinations).T
full_range_df.columns = ['chr', 'strand', 'start', 'end', 'probs']

# Concat both dfs and keep first occurence of duplicates (first being the real preds, not the missing ones)
concat_df = pd.concat([surrounding_preds, full_range_df]).drop_duplicates(subset=['start', 'end', 'strand'], keep='first').sort_values(by=['start', 'end', 'strand'])
concat_df = concat_df[concat_df['strand'] == cd_strand]
#concat_df.to_csv('win_test.tsv', index=False, sep='\t')

# Create lineplot
rc = {'ytick.labelsize': 25, 'xtick.labelsize': 25}
plt.rcParams.update(**rc)
plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(1, 1, figsize=(12, 10))
ax.plot([-i for i in reversed(range(1, 100+1))] + [i for i in range(0, 100+1)], 
        list(concat_df.probs), color='black', label=f'{cd_cerevisiae}\nstrand {cd_strand}\n{sno_len} nt')

ax.set_xlabel('Distance to original window (nt)', fontdict={'fontsize': 30})
ax.set_ylabel('Prediction probability of windows\nsurrounding the C/D snoRNA (%)', fontdict={'fontsize': 30})
plt.legend()
fig.suptitle(f'Window probability around the\nC/D snoRNA {cd_cerevisiae} in S. cerevisiae', fontsize=35, weight='bold', x=0.5, y=1)
plt.savefig(output, dpi=600, bbox_inches='tight')

