#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import os
import subprocess as sp
from functools import reduce

# Load and format inputs and outputs
fixed_length = snakemake.params.fixed_length
species = snakemake.wildcards.species
bedgraph_path = snakemake.input.bedgraph
bedgraph_names = os.listdir(bedgraph_path)
bedgraphs = [bedgraph_path+'/'+p for p in bedgraph_names]
bedgraphs_sorted = [
    bedgraph_path+'/'+'final_'+p for p in bedgraph_names if p.endswith(
        '.bedgraph')]
simple_bedgraph_names = [
    b for b in bedgraph_names if 'final' not in b and 'sorted' not in b]
for i, b in enumerate(bedgraphs):  # remove unwanted last line
    sp.call(f"sed '/type=bedGraph/d' {b} | LC_COLLATE=C sort -k1,1 -k2,2n > {bedgraphs_sorted[i]}", 
    shell=True)
pred_dir = snakemake.input.snoBIRD_preds
preds = pd.read_csv(pred_dir+'/'+os.listdir(pred_dir)[0], sep='\t')
if species != 'tetrahymena_thermophila':
    preds['chr'] = 'chr' + preds['chr'].astype(str)
else:
    preds['chr'] = preds['chr'].astype(str)

preds['gene_length'] = preds['end'] - preds['start']
bed_cols = ['chr', 'start', 'end', 'gene_id', 'gene_length', 'strand']
preds[bed_cols].to_csv(f'preds_{species}.bed', index=False, sep='\t', 
        header=False)
preds_bed = BedTool(f'preds_{species}.bed')
output_df = snakemake.output.coverage

# Create bed of known annotated C/D that fit within the range of detection 
# of snoBIRD
if species == 'tetrahymena_thermophila':
    cd_df = pd.read_csv(snakemake.input.known_cd_tetrahymena, sep='\t', dtype={'chr': 'str'})
else:
    cd_df = pd.read_csv(snakemake.input.known_cd, sep='\t', dtype={'chr': 'str'})
cd_df = cd_df[cd_df['sno_type'] == 'C/D']    
cd_df['start'] = cd_df['start'].astype(int)
cd_df['end'] = cd_df['end'].astype(int)
if species == 'schizosaccharomyces_pombe':
    cd_df['chr'] = 'chr' + cd_df['chr'].astype(str)
if species == 'homo_sapiens':
    cd_df = cd_df.dropna(subset=['gene_id'])
cd_df = cd_df[(cd_df['end'] - cd_df['start']) <= fixed_length - 30]
cd_df['dot'] = '.'
cd_df[['chr', 'start', 'end', 'gene_id', 
            'dot', 'strand', 'sno_type']].to_csv(f'{species}_known_cd.bed', 
            sep='\t', index=False, header=False)

bed_cd = BedTool(f'{species}_known_cd.bed')

#if species == 'drosophila_melanogaster':
#    # Get regions upstream and downstream of predictions
#    sp.call(
#        """awk -v OFS="\t" '{print $1, $2-11, $2-1, $4, $5, $6}' preds_"""+
#        f"""{species}.bed | sort -k1,1 -k2,2n > preds_{species}_upstream.bed""", 
#        shell=True)
#    sp.call(
#        """awk -v OFS="\t" '{print $1, $3+1, $3+11, $4, $5, $6}' preds_"""+
#        f"""{species}.bed | sort -k1,1 -k2,2n > preds_{species}_downstream.bed""", 
#        shell=True)
#
#    up_bed = BedTool(f'preds_{species}_upstream.bed')
#    down_bed = BedTool(f'preds_{species}_downstream.bed')
#
#    # Find avg for these up/downstream regions
#    ups, downs = [], []
#    for bg in bedgraphs_sorted:
#        up = up_bed.map(b=bg, c=4, o="mean").to_dataframe(
#            names=bed_cols+['avg_upstream_cov'])
#        up['avg_upstream_cov'] = up['avg_upstream_cov'].replace('.', 0)
#        down = down_bed.map(b=bg, c=4, o="mean").to_dataframe(
#                            names=bed_cols+['avg_downstream_cov'])
#        down['avg_downstream_cov'] = down['avg_downstream_cov'].replace('.', 0)
#        up = up.rename(columns={
#            'avg_upstream_cov': 'avg_upstream_cov_'+bg.split('final_')[-1].split(
#                '_sorted')[0]})
#        down = down.rename(columns={
#            'avg_downstream_cov': 'avg_downstream_cov_'+bg.split(
#                'final_')[-1].split('_sorted')[0]})
#        ups.append(up.filter(regex='gene_id|^avg_upstream*'))
#        downs.append(down.filter(regex='gene_id|^avg_downstream*'))
#
#    # Find avg up/down across samples
#    all_ups = reduce(lambda left, right: pd.merge(left, right, on='gene_id', 
#                        how='left'), ups)
#    all_ups[[c for c in all_ups.columns if 'avg_up' in c]] = all_ups.filter(
#                                    regex='^avg_upstream').astype(float)
#    all_ups['avg_upstream_cov'] = all_ups.filter(regex='^avg_upstream').mean(axis=1)
#
#    all_downs = reduce(lambda left, right: pd.merge(left, right, on='gene_id', 
#                        how='left'), downs)
#    all_downs[[c for c in all_downs.columns if 'avg_down' in c]] = all_downs.filter(
#                                    regex='^avg_downstream').astype(float)
#    all_downs['avg_downstream_cov'] = all_downs.filter(
#                            regex='^avg_downstream').mean(axis=1)
#

# Intersect predictions with bedgraphs (bg) of expression 
val = preds_bed.intersect(b=bedgraphs, wa=True, wb=True, 
                        names=simple_bedgraph_names)
col_names = bed_cols + ['bg_name', 'chr_bg', 'start_bg', 'end_bg', 'coverage']
#print(val.to_dataframe())
df = val.to_dataframe(names=col_names).sort_values(['gene_id', 'bg_name'])
#print(pd.unique(df.chr))

# Find per prediction what is the average coverage over the whole sequence 
df['cov_length'] = df['end_bg'] - df['start_bg']
df['weighted_coverage'] = df['cov_length'] * df['coverage']

groupby_cols = ['gene_id', 'chr', 'start', 'end', 'strand', 'bg_name', 
                'gene_length']
# find the sum of bedgraph coverage that intersects with the prediction
grouped_df = df.groupby(groupby_cols).agg(total_weighted_coverage=(
    'weighted_coverage', 'sum'), total_cov_length=('cov_length', 'sum')
    ).reset_index()
## flag uneven coverage (coverage < 50 % of pred length)
#uneven_id = list(pd.unique(grouped_df[
#    grouped_df['total_cov_length'] < 0.9 * grouped_df['gene_length']][
#        'gene_id']))
## flag expanded coverage surrounding prediction (so not really expressed, but 
## the overlapping transcript is)
#expanded_cov_id = list(pd.unique(grouped_df[
#    grouped_df['total_cov_length'] > 20 + grouped_df['gene_length']][
#        'gene_id']))
# get avg coverage over the prediction total length
grouped_df['avg_coverage'] = grouped_df[
                            'total_weighted_coverage'] / grouped_df[
                                                                'gene_length']

# Pivot df to create 1 avg_coverage col per bedgraph
pivot_df = grouped_df.pivot_table(index=[
                        'gene_id', 'chr', 'start', 'end', 'strand'], 
                        columns='bg_name', values='avg_coverage', 
                        fill_value=0).reset_index()

# Avg coverage across all samples/bedgraphs
pivot_df['avg_coverage_samples'] = pivot_df.filter(regex='bedgraph$').mean(
                                                                        axis=1)
pivot_df = pivot_df.rename(columns=lambda name: name.split(
                    '.bedgraph')[0] if name.endswith('.bedgraph') else name)

# Concat missing predictions (not detected in any condition)
not_expressed_pred = preds[~preds['gene_id'].isin(pivot_df.gene_id)]
not_expressed_pred = not_expressed_pred[['gene_id', 'chr', 'start', 
                                        'end', 'strand']]
new_cols = [n.split(
    '.bedgraph')[0] for n in simple_bedgraph_names] + ['avg_coverage_samples']
not_expressed_pred[new_cols] = 0  # give 0 for coverage value for these 
                                  # non-detected predictions
    
final_df = pd.concat([pivot_df, not_expressed_pred]).sort_values(
                                            ['chr', 'strand', 'start', 'end'])



# Add filters to consider if prediction is expressed and a potential snoRNA
exp_preds = final_df[final_df['avg_coverage_samples'] >= 5]
#print(exp_preds)
exp_preds['dot'] = '.'
exp_preds[['chr', 'start', 'end', 'gene_id', 
            'dot', 'strand', 'avg_coverage_samples']].to_csv(f'exp_preds_{species}.bed', index=False, sep='\t', 
        header=False)

exp_preds_bed = BedTool(f'exp_preds_{species}.bed')

# Find overlap between annotated C/D and expressed predictions
pos_overlap = bed_cd.intersect(exp_preds_bed, s=True, wb=True, 
                f=0.5).to_dataframe()
num_pos = len(pos_overlap)
#print(pos_overlap)


print(species)
print(f'Total number of SnoBIRD predictions: {len(final_df)}')
print('Total number of expressed SnoBIRD predictions (avg coverage >=5): '+
        f'{len(exp_preds)},')
print(f'with {num_pos} being already annotated and {len(exp_preds) - num_pos} potential novel C/D snoRNAs')
print(f'Total number of known annotated C/D snoRNAs: {len(cd_df)}')

# Info df
cols = ['species', 'total_annotated_cd', 'total_snoBIRD_predictions', 'expressed_snoBIRD_predictions (>=5 reads)', 
        'expressed_annotated_and_predicted', 'potential_novel_expressed']
final_vals = [species, len(cd_df), len(final_df), len(exp_preds), num_pos, len(exp_preds) - num_pos]
temp_df = pd.DataFrame([final_vals], columns=cols)
temp_df.to_csv(snakemake.output.coverage, sep='\t', index=False)

"""
if species == 'drosophila_melanogaster':
    # length of avg_cov must cover at least 50% of the prediction
    exp_preds = exp_preds[~exp_preds['gene_id'].isin(uneven_id)]
    print('Total number of expressed SnoBIRD predictions (even coverage >=5): '+
            f'{len(exp_preds)}')

    # length of avg_cov must cover at least 50% of the prediction
    exp_preds = exp_preds[~exp_preds['gene_id'].isin(expanded_cov_id)]
    print('Total number of expressed SnoBIRD predictions (not expanded coverage):'+
            f' {len(exp_preds)}')

    # avg_cov must of 10 nt up/downstream regions must be <=25% of the predicted 
    # region, otherwise most likely not a snoRNA but a prediction within a longer 
    # expressed transcript
    # flag if up/downstream coverage is similar to avg coverage, meaning it's not 
    # a block so most likely not a snoRNA

    up_down_df = final_df.merge(all_ups[['gene_id', 'avg_upstream_cov']], 
                how='left', on='gene_id')
    up_down_df = up_down_df.merge(all_downs[['gene_id', 'avg_downstream_cov']],
                how='left', on='gene_id')
    up_down_id = list(pd.unique(up_down_df[(
        (up_down_df['avg_upstream_cov'] / up_down_df['avg_coverage_samples']
        ) <= 0.25) & (
        (up_down_df['avg_downstream_cov'] / up_down_df['avg_coverage_samples']
        ) <= 0.25)]['gene_id']))
    exp_preds = exp_preds[exp_preds['gene_id'].isin(up_down_id)]
    print('Total number of expressed SnoBIRD predictions (block of reads): '+
            f'{len(exp_preds)}')

final_df.loc[final_df['gene_id'].isin(exp_preds.gene_id), 'real_exp'] = 'real_exp'
final_df['real_exp'] = final_df['real_exp'].fillna('uneven_or_expanded')
final_df.to_csv(output_df, sep='\t', index=False)
"""
sp.call(f'rm preds_{species}* exp_preds_{species}*', shell=True)
sp.call(f'rm results/bedgraph_TGIRT/{species}/final_*', shell=True)
