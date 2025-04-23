#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import matplotlib.pyplot as plt
from upsetplot import UpSet
from glob import glob
import subprocess as sbp

# Load beds of predictions
bed_files_human, bed_files_pombe = {}, {}
pombe_dfs, human_dfs = {}, {}
snoBIRD_human = pd.read_csv(glob([p for p in snakemake.input.snoBIRD_preds if 'homo_sapiens' in p][0] + '/*tsv')[0], sep='\t')
snoBIRD_pombe = pd.read_csv(glob([p for p in snakemake.input.snoBIRD_preds if 'schizosaccharomyces_pombe' in p][0] + '/*tsv')[0], sep='\t')
snoBIRD_human['chr'] = 'chr' + snoBIRD_human['chr'].astype(str)
snoBIRD_pombe['chr'] = 'chr' + snoBIRD_pombe['chr'].astype(str)
snoBIRD_human['score'] = '.'
snoBIRD_pombe['score'] = '.'
human_dfs['snoBIRD'] = snoBIRD_human
pombe_dfs['snoBIRD'] = snoBIRD_pombe


bed_files_human['snoBIRD'] = BedTool.from_dataframe(snoBIRD_human[['chr', 'start', 'end', 'gene_id', 'score', 'strand']])
bed_files_pombe['snoBIRD'] = BedTool.from_dataframe(snoBIRD_pombe[['chr', 'start', 'end', 'gene_id', 'score', 'strand']])
for p in snakemake.input.other_preds:
    tool_ = p.split('/')[2]
    sbp.call(f'sort -k1,1 -k2,2n {p} > temp', shell=True)
    bed = BedTool('temp')
    merged_bed = bed.merge(c=[4,5,6], o=['distinct', 'distinct', 'distinct'], delim='_', s=True)
    if 'homo_sapiens' in p:
        bed_files_human[tool_] = merged_bed
        human_dfs[tool_] = merged_bed.to_dataframe(names=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])
    elif 'schizosaccharomyces_pombe' in p:
        bed_files_pombe[tool_] = merged_bed
        pombe_dfs[tool_] = merged_bed.to_dataframe(names=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])
    sbp.call('rm temp', shell=True)

# Compute overlaps (50% reciprocal) for human and pombe
overlap_matrix_pombe, overlap_matrix_human = {}, {}
bed_files = [bed_files_pombe, bed_files_human]
for i, sp in enumerate(['schizosaccharomyces_pombe']):
    print(bed_files[i]['snoBIRD'].to_dataframe())
    intersect_snoBIRD = bed_files[i]['snoBIRD'].intersect([v for k,v in bed_files[i].items() if k != 'snoBIRD'], f=0.5, wa=True, wb=True, s=True, loj=True)
    print(intersect_snoBIRD.to_dataframe(names=['chr', 'start', 'end', 'gene_id', 'score', 'strand', 'chr1', 'start1', 'end1', 'gene_id1', 'score1', 'strand1']).drop_duplicates('gene_id'))
    intersect_snoBIRD_merge = intersect_snoBIRD.merge(c=[4, 10], o=['distinct', 'distinct'], delim=';', s=True) # delimit tools that overlap the same snoBIRD pred
    intersect_df_snoBIRD = intersect_snoBIRD_merge.to_dataframe(names=
        ['chr', 'start', 'end', 'gene_id', 'gene_id_other'])
    print(snoBIRD_pombe[~snoBIRD_pombe['gene_id'].isin(intersect_df_snoBIRD.gene_id)])
    print(intersect_df_snoBIRD)
    print(intersect_df_snoBIRD[intersect_df_snoBIRD.gene_id.str.contains(';')])
#seen_tool = []
#for tool1, bed1 in bed_files_pombe.items():
#    for tool2, bed2 in bed_files_pombe.items():
#        if (tool1 == tool2) | (tool2 in seen_tool):
#            continue
#        # Intersect the bed with at least 50% reciprocal overlap
#        intersect = bed1.intersect(bed2, f=0.5, wa=True, wb=True, s=True)
#        intersect_df = intersect.to_dataframe(names=
#            [f'chr_{tool1}', f'start_{tool1}', f'end_{tool1}', f'gene_id_{tool1}', f'score_{tool1}', f'strand_{tool1}',
#            f'chr_{tool2}', f'start_{tool2}', f'end_{tool2}', f'gene_id_{tool2}', f'score_{tool2}', f'strand_{tool2}'])
#        print(tool1, tool2)
#        print(intersect_df)
#        overlap_matrix_pombe[tool1] = intersect_df  # Store intersected predictions
#    seen_tool.append(tool1)

# Convert to binary matrix for UpSet plot
binary_matrix = []
index_labels = []

for tool, overlaps in overlap_matrix_pombe.items():
    unique_entries = set.union(*overlaps) if overlaps else set()
    for entry in unique_entries:
        row = [entry in overlap for overlap in overlaps]
        binary_matrix.append([tool] + row)
        index_labels.append(entry)

# Convert to DataFrame
columns =  list(bed_files_pombe.keys())
df = pd.DataFrame(binary_matrix, columns=columns)
print(df)

# Convert to numeric (1 for presence, 0 for absence)
df = df.astype(int)

# Create UpSet plot
upset = UpSet(df, subset_size="sum", show_percentages=True)
upset.plot()
plt.savefig('upset_pombe_tools.svg', dpi=600, bbox_inches='tight')



