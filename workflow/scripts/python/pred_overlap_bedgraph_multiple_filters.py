#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import os
import subprocess as sp
from functools import reduce
import warnings 
warnings.filterwarnings("ignore")

# Load and format inputs and params
low_rep_species = ['schizosaccharomyces_pombe', 'macaca_mulatta', 'gallus_gallus', 'tetrahymena_thermophila', 'arabidopsis_thaliana']  # species with low number of different samples (use .all())
high_rep_species = ['homo_sapiens', 'drosophila_melanogaster', 'plasmodium_falciparum', 'danio_rerio']  # large number of conditions, so not expressed necessarily in all
screenshot_dir = str(snakemake.params.screenshot_dir)
gtf = snakemake.input.gtf
chr_size = snakemake.input.chr_size
fixed_length = snakemake.params.fixed_length
species = snakemake.wildcards.species
bedgraph_path = snakemake.input.bedgraph
bedgraph_names = [b for b in os.listdir(bedgraph_path) if 'final_' not in b]
bedgraphs_pos = [bedgraph_path+'/'+p for p in bedgraph_names if '_fwd' in p]
bedgraphs_neg = [bedgraph_path+'/'+p for p in bedgraph_names if '_rev' in p]
bedgraphs_sorted_pos = [
    bedgraph_path+'/'+'final_'+p for p in bedgraph_names if p.endswith(
        'fwd.bedgraph')]
bedgraphs_sorted_neg = [
    bedgraph_path+'/'+'final_'+p for p in bedgraph_names if p.endswith(
        'rev.bedgraph')]
simple_bedgraph_names_pos = [
    b for b in bedgraph_names if 'final' not in b and 'sorted' not in b and '_fwd' in b]
simple_bedgraph_names_neg = [
    b for b in bedgraph_names if 'final' not in b and 'sorted' not in b and '_rev' in b]

print('POS')
if species == 'homo_sapiens':
    for i, b in enumerate(bedgraphs_pos):
        print(b)
        sp.call(f"LC_COLLATE=C sort -V {b} > {bedgraphs_sorted_pos[i]}", 
        shell=True)
    for i, b in enumerate(bedgraphs_neg):  
        print(b)
        sp.call(f"LC_COLLATE=C sort -V {b} > {bedgraphs_sorted_neg[i]}", 
        shell=True)
else:
    for i, b in enumerate(bedgraphs_pos):  # remove unwanted last line
        sp.call(f"sed '/type=bedGraph/d' {b} | LC_COLLATE=C sort -k1,1 -k2,2n > {bedgraphs_sorted_pos[i]}", 
        shell=True)
    for i, b in enumerate(bedgraphs_neg):  # remove unwanted last line
        sp.call(f"sed '/type=bedGraph/d' {b} | LC_COLLATE=C sort -k1,1 -k2,2n > {bedgraphs_sorted_neg[i]}", 
        shell=True)
pred_dir = snakemake.input.snoBIRD_preds
preds = pd.read_csv(pred_dir+'/'+os.listdir(pred_dir)[0], sep='\t', 
            dtype={'chr': 'str'})
if species != 'tetrahymena_thermophila':
    preds['chr'] = 'chr' + preds['chr'].astype(str)
else:
    preds['chr'] = preds['chr'].astype(str)

preds['gene_length'] = preds['end'] - preds['start']
bed_cols = ['chr', 'start', 'end', 'gene_id', 'gene_length', 'strand']
preds_bed = BedTool.from_dataframe(preds[bed_cols])
preds[preds['strand'] == '+'][bed_cols].to_csv('preds_pos.tsv', sep='\t', index=False, header=False)
preds[preds['strand'] == '-'][bed_cols].to_csv('preds_neg.tsv', sep='\t', index=False, header=False)
sp.call('sort -V preds_pos.tsv > preds_pos_sorted.tsv', shell=True)
sp.call('sort -V preds_neg.tsv > preds_neg_sorted.tsv', shell=True)


preds_bed_pos = BedTool('preds_pos_sorted.tsv')
preds_bed_neg = BedTool('preds_neg_sorted.tsv')

print('FINISHED SORTING')
def merge_dicts_with_max(d, other_dict):
    # Merge key:val pairs in an existing dict only if val> existing val
    for key, val in other_dict.items():
        if key not in d or val > d[key]:
            d[key] = val
    return d


# Create bed of known annotated C/D in RNAcentral that fit within the range of detection 
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

if species in ['arabidopsis_thaliana', 'danio_rerio']:
    # the TGIRT-Seq samples of these species were size-selected for 
    # RNAs <= 100 nt (designed for tRNAs initially)
    cd_df = cd_df[(cd_df['end'] - cd_df['start']) <= 100]


cd_df['dot'] = '.'
cd_df = cd_df[['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'sno_type']]


# Add also snoRNAs from the ensembl gtf that might complement what is missing in RNAcentral 
# (except for human (already complete from snoDB) and T. thermophila (no snoRNA annotated in ensembl gtf))
# Remove H/ACAs starting with "SNORA" and snoRNAs with length>fixed_length+30
if species not in ['homo_sapiens', 'tetrahymena_thermophila', 'schizosaccharomyces_pombe']:
    col = '$14'
    names = "SNORD|U3|U8"
    if species == 'drosophila_melanogaster':
        col = '$12'
        names = "CD|snoRNA:Me"
    awk_cmd = """awk -v OFS='\t' 'NR>6 && $3=="gene" && $0 ~ "gene_name" && $0 ~ "gene_biotype .snoRNA." && $0 !~ "SNORA" && $5-$4+1<"""+ \
            str(fixed_length-30) + f""" {{print "chr"$1,$4,$5,$10,$6,$7,"C/D_potentially",{col}}}' """+ gtf + \
            """ | sed 's/"//g; s/;//g' > """+ f"""{species}_ensembl.bed"""
    sp.call(awk_cmd, shell=True)
    if species == "gallus_gallus":
    # Add these 2 known snoRNAs as CD (they don't have gene_names in the gtf)
        other_cd = """awk -v OFS='\t' 'NR>6 && $3=="gene" && $0 ~ "ENSGALG00010021006|ENSGALG00010018970" """+ \
            """{print "chr"$1,$4,$5,$10,$6,$7,"C/D_potentially",$10}' """+ gtf + \
            """ | sed 's/"//g; s/;//g' >> """+ f"""{species}_ensembl.bed"""
        sp.call(other_cd, shell=True)
    ensembl_cd = pd.read_csv(f"{species}_ensembl.bed", 
        names=['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'sno_type', 'gene_name'], sep='\t')
    # Filter based on gene_name to keep only validated C/D
    ensembl_cd = ensembl_cd[(ensembl_cd['gene_name'].str.contains(names)) | (
                ensembl_cd['gene_id'].isin(["ENSGALG00010021006","ENSGALG00010018970"]))]
    ensembl_cd = ensembl_cd.drop(columns='gene_name')
    cd_df = pd.concat([cd_df, ensembl_cd])


# Remove duplicates and overlapping annotated snoRNAs with slightly different coordinates (overlap at >=75%)
cd_df = cd_df.drop_duplicates(['chr', 'start', 'end', 'strand']).sort_values(
        ['chr', 'strand', 'start']).reset_index(drop=True)
def filter_group(group):
    """ Removes overlapping snoRNA entries."""
    non_overlapping = []
    for i, row in group.iterrows():
        is_overlapping = False
        for j, prev_row in enumerate(non_overlapping):
            # Calculate overlap percentage
            overlap_start = max(row["start"], prev_row["start"])
            overlap_end = min(row["end"], prev_row["end"])
            overlap_len = max(0, overlap_end - overlap_start)
            row_len = row["end"] - row["start"]
            if overlap_len / row_len >= 0.75:
                is_overlapping = True
                break
        if not is_overlapping:
            non_overlapping.append(row)
    return pd.DataFrame(non_overlapping)


cd_df = cd_df.groupby(['chr', 'strand'], group_keys=False).apply(filter_group)

cd_df_pos = cd_df[cd_df['strand'] == '+'].sort_values(
                                    by=['chr', 'start']).reset_index(drop=True)
cd_df_neg = cd_df[cd_df['strand'] == '-'].sort_values(
                                    by=['chr', 'start']).reset_index(drop=True)

print(cd_df_pos)
print(cd_df_neg)

cd_df_pos.to_csv('cd_test_pos.tsv', sep='\t', index=False, header=False)
cd_df_neg.to_csv('cd_test_neg.tsv', sep='\t', index=False, header=False)
sp.call('sort -V cd_test_pos.tsv > cd_test_pos_sorted.tsv', shell=True)
sp.call('sort -V cd_test_neg.tsv > cd_test_neg_sorted.tsv', shell=True)

bed_cd = BedTool.from_dataframe(cd_df)
bed_cd_pos = BedTool('cd_test_pos_sorted.tsv')
bed_cd_neg = BedTool('cd_test_neg_sorted.tsv')
print('POSNEG')
print(len(bedgraphs_pos))

# Intersect annotated C/D and bedgraph to see how many annotated C/D are expressed
col_names = ['chr', 'start', 'end', 'known_cd_id', 'score', 'strand', 
                'sno_type', 'bg_name', 'chr_bg', 'start_bg', 'end_bg', 'coverage']
if species == 'homo_sapiens':
    annot_pos = bed_cd_pos.intersect(b=bedgraphs_sorted_pos, wa=True, wb=True, 
                            names=simple_bedgraph_names_pos, sorted=True, g=chr_size)

    annot_neg = bed_cd_neg.intersect(b=bedgraphs_sorted_neg, wa=True, wb=True, 
                            names=simple_bedgraph_names_neg, sorted=True, g=chr_size)    

else:
    annot_pos = bed_cd_pos.intersect(b=bedgraphs_sorted_pos, wa=True, wb=True, 
                            names=simple_bedgraph_names_pos)

    annot_neg = bed_cd_neg.intersect(b=bedgraphs_sorted_neg, wa=True, wb=True, 
                            names=simple_bedgraph_names_neg)    

    
annot_df = pd.concat([annot_pos.to_dataframe(names=col_names), 
            annot_neg.to_dataframe(names=col_names)]).reset_index(drop=True)
    
print(annot_pos.to_dataframe())
print(annot_neg.to_dataframe())
print(annot_df)
print('intersecttCD')
# Correct start/end of bedgraph so that we keep only fully overlapping regions 
# with the known C/D (not extended after the snoRNA)
annot_df.loc[annot_df['end_bg'] > annot_df['end'], 'end_bg'] = annot_df['end']
annot_df.loc[annot_df['start_bg'] < annot_df['start'], 'start_bg'] = annot_df['start']

# Find per prediction what is the average coverage over the whole sequence 
annot_df['cov_length'] = annot_df['end_bg'] - annot_df['start_bg']
annot_df['weighted_coverage'] = annot_df['cov_length'] * annot_df['coverage']


# Find also the longest stretch of consecutive coverage over the prediction 
# Merge overlapping bedgraph intervals nd keep the longest stretch across a given snoRNA
cov_dict_cd = {}
for n in pd.unique(annot_df.bg_name):
    temp_cd = annot_df[annot_df['bg_name'] == n]
    temp_bed_cd = BedTool.from_dataframe(temp_cd[['known_cd_id', 'start_bg', 'end_bg', 'chr', 
            'coverage', 'strand', 'bg_name']].sort_values(
            by=['known_cd_id', 'start_bg']).reset_index(drop=True))
    # Merge based on the known_cd_id (as "chr" for the 1st bed column) so that we 
    # merge per snoRNA and not across snoRNAs across a given chr
    merged_cd = temp_bed_cd.merge(s=True, c=[4], o=['distinct'])
    merged_cd = merged_cd.to_dataframe(names=['known_cd_id', 'start_bg', 'end_bg', 'chr'])
    merged_cd['stretch_len'] = merged_cd.end_bg - merged_cd.start_bg
    cov_dict_cd = merge_dicts_with_max(cov_dict_cd, merged_cd.groupby('known_cd_id')['stretch_len'].max().to_dict())


# Find the sum of bedgraph coverage that intersects with the prediction 
# (i.e. sum all parts of bedgraphs that overlap entirely with the snoRNA)
groupby_cols = ['known_cd_id', 'chr', 'start', 'end', 'strand', 'bg_name']
annot_grouped_df = annot_df.groupby(groupby_cols).agg(total_weighted_coverage=(
    'weighted_coverage', 'sum'), total_cov_length=('cov_length', 'sum')
    ).reset_index()
annot_grouped_df['gene_length'] = annot_grouped_df.end - annot_grouped_df.start + 1

""" Apply 1st filter for annotated C/D: read must cover at least 50% of the snoRNA 
    length (for it to be assigned to the C/D) and the longest stretch of consecutive 
    reads must also cover at least 50% of the snoRNA predicted length"""
annot_grouped_df = annot_grouped_df[
            annot_grouped_df['total_cov_length'] >= annot_grouped_df['gene_length'] * 0.5]
annot_grouped_df['max_stretch_cov'] = annot_grouped_df['known_cd_id'].map(cov_dict_cd)
annot_grouped_df = annot_grouped_df[annot_grouped_df['max_stretch_cov'] >= annot_grouped_df['gene_length'] * 0.5]
annot_grouped_df = annot_grouped_df.drop(columns='max_stretch_cov')

# Avg coverage per bedgraph
annot_grouped_df['avg_coverage'] = annot_grouped_df[
                            'total_weighted_coverage'] / annot_grouped_df[
                                                                'gene_length']
                                                          

# Pivot df to create 1 avg_coverage col per bedgraph
annot_pivot_df = annot_grouped_df.pivot_table(index=[
                        'known_cd_id', 'chr', 'start', 'end', 'strand'], 
                        columns='bg_name', values='avg_coverage', 
                        fill_value=0).reset_index()

""" Apply 2nd filter for annotated C/D: Avg per condition must be at least >0 read or at least 
    one avg > 10 in the case of multiple samples per species (that ensures 
    consistency across replicates).But only based either on fwd (+ strand) or rev (- strand) bedgraphs, 
    because, if expressed from strand +, we don't expect read coverage on the opposite strand."""
if species in low_rep_species:  # low nb of TGIRT-Seq samples(<=2 conditions), thus we expect the C/D to be expressed in all samples
    annot_pivot_df = annot_pivot_df[((annot_pivot_df.filter(regex='fwd.bedgraph$') > 0).all(axis=1)) | 
                                ((annot_pivot_df.filter(regex='rev.bedgraph$') > 0).all(axis=1))]
if species in high_rep_species:  
    # High nb of TGIRT-Seq samples (>=3 conditions), thus at least one avg 
    # condition should be > 10 (and later, the 3rd filter checks if overall avg is > 5)
    _names = pd.unique([b.split('.')[0].rsplit('_', maxsplit=2)[0] for b in bedgraph_names if '_fwd' in b])
    for n in _names:
        annot_pivot_df[n+'_fwd__AVG'] = annot_pivot_df.filter(
                            regex=n+'_[0-9]*(.bam)?_fwd.bedgraph$').mean(axis=1)
        annot_pivot_df[n+'_rev__AVG'] = annot_pivot_df.filter(
                            regex=n+'_[0-9]*(.bam)?_rev.bedgraph$').mean(axis=1)
    annot_pivot_df = annot_pivot_df[(annot_pivot_df.filter(regex='__AVG$') > 10).any(axis=1)]
    annot_pivot_df = annot_pivot_df.drop(columns=[c for c in annot_pivot_df.columns if '__AVG' in c])

# Avg coverage across all samples/bedgraphs in which the C/D is expressed > 0
bg_cols = annot_pivot_df.filter(regex='bedgraph$')
annot_pivot_df['avg_coverage_samples'] = bg_cols.where(bg_cols > 0).mean(axis=1)
annot_pivot_df = annot_pivot_df.rename(columns=lambda name: name.split(
                    '.bedgraph')[0] if name.endswith('.bedgraph') else name)
    
""" Apply 3rd filter for annotated C/D: Avg across expressed conditions must be at least >5 reads"""
# Add filters to consider if prediction is expressed and a potential snoRNA
#annot_pivot_df = annot_pivot_df[annot_pivot_df['avg_coverage_samples'] >= 5]
num_annot_expressed = len(annot_pivot_df)




# Intersect annotated C/D and predictions (regardless of expression here)
annotated_and_pred = bed_cd.intersect(preds_bed, s=True, wb=True, 
                f=0.5).to_dataframe(names=['chr', 'start', 'end', 
                'known_cd_id', 'score', 'strand', 'sno_type', 'chr_pred', 
                'start_pred', 'end_pred', 'gene_id', 'gene_length_pred', 
                'strand_pred'])
num_annotated_and_pred = len(annotated_and_pred)


print('FINISHED KNOWN CD')

# Intersect predictions with bedgraphs (bg) of expression
if species == 'homo_sapiens': 
    val_pos = preds_bed_pos.intersect(b=bedgraphs_sorted_pos, wa=True, wb=True, 
                            names=simple_bedgraph_names_pos, sorted=True, g=chr_size)
    val_neg = preds_bed_neg.intersect(b=bedgraphs_sorted_neg, wa=True, wb=True, 
                            names=simple_bedgraph_names_neg, sorted=True, g=chr_size)
else:
    val_pos = preds_bed_pos.intersect(b=bedgraphs_sorted_pos, wa=True, wb=True, 
                            names=simple_bedgraph_names_pos)
    val_neg = preds_bed_neg.intersect(b=bedgraphs_sorted_neg, wa=True, wb=True, 
                            names=simple_bedgraph_names_neg)
col_names = bed_cols + ['bg_name', 'chr_bg', 'start_bg', 'end_bg', 'coverage']
df_pos = val_pos.to_dataframe(names=col_names)
df_neg = val_neg.to_dataframe(names=col_names)
df = pd.concat([df_pos, df_neg]).sort_values(
                ['gene_id', 'bg_name']).reset_index(drop=True)

# Correct start/end of bedgraph so that we keep only fully overlapping regions 
# with the predictions (not extended after the pred)
df.loc[df['end_bg'] > df['end'], 'end_bg'] = df['end']
df.loc[df['start_bg'] < df['start'], 'start_bg'] = df['start']

# Find per prediction what is the average coverage over the whole sequence 
df['cov_length'] = df['end_bg'] - df['start_bg']
df['weighted_coverage'] = df['cov_length'] * df['coverage']


# Find also the longest stretch of consecutive coverage over the prediction 
# Merge overlapping bedgraph intervals nd keep the longest stretch across a given snoRNA
cov_dict = {}
for n in pd.unique(df.bg_name):
    temp_ = df[df['bg_name'] == n]
    temp_bed = BedTool.from_dataframe(temp_[['gene_id', 'start_bg', 'end_bg', 'chr', 
            'coverage', 'strand', 'bg_name']].sort_values(
            by=['gene_id', 'start_bg']).reset_index(drop=True))
    # Merge based on the gene_id (as "chr" for the 1st bed column) so that we 
    # merge per snoRNA and not across snoRNAs across a given chr
    merged = temp_bed.merge(s=True, c=[4], o=['distinct'])
    merged = merged.to_dataframe(names=['gene_id', 'start_bg', 'end_bg', 'chr'])
    merged['stretch_len'] = merged.end_bg - merged.start_bg
    cov_dict = merge_dicts_with_max(cov_dict, merged.groupby('gene_id')['stretch_len'].max().to_dict())

# find the sum of bedgraph coverage that intersects with the prediction
# (i.e. sum all parts of bedgraphs that overlap entirely with the prediction)
groupby_cols = ['gene_id', 'chr', 'start', 'end', 'strand', 'bg_name', 
                'gene_length']
grouped_df = df.groupby(groupby_cols).agg(total_weighted_coverage=(
    'weighted_coverage', 'sum'), total_cov_length=('cov_length', 'sum')
    ).reset_index()


""" Apply 1st filter: reads must cover at least 50% of the snoRNA 
    predicted length and the longest stretch of consecutive reads must also cover 
    at least 50% of the snoRNA predicted length"""
grouped_df = grouped_df[
            grouped_df['total_cov_length'] >= grouped_df['gene_length'] * 0.5]

grouped_df['max_stretch_cov'] = grouped_df['gene_id'].map(cov_dict)
grouped_df = grouped_df[grouped_df['max_stretch_cov'] >= grouped_df['gene_length'] * 0.5]
grouped_df = grouped_df.drop(columns='max_stretch_cov')
# get avg coverage (per bedgraph) over the prediction total length
grouped_df['avg_coverage'] = grouped_df[
                            'total_weighted_coverage'] / grouped_df[
                                                                'gene_length']
# Pivot df to create 1 avg_coverage col per bedgraph
pivot_df = grouped_df.pivot_table(index=[
                        'gene_id', 'chr', 'start', 'end', 'strand'], 
                        columns='bg_name', values='avg_coverage', 
                        fill_value=0).reset_index()
                    
""" Apply 2nd filter: Avg per condition must be at least >0 read or at least 
    one avg > 10 in the case of multiple samples per species (that ensures 
    consistency across replicates).But only based either on fwd (+ strand) or rev (- strand) bedgraphs, 
    because, if expressed from strand +, we don't expect read coverage on the opposite strand."""

if species in low_rep_species:  # low nb of TGIRT-Seq samples(<=2 conditions), thus we expect the C/D to be expressed in all samples
    pivot_df = pivot_df[((pivot_df.filter(regex='fwd.bedgraph$') > 0).all(axis=1)) | 
                                ((pivot_df.filter(regex='rev.bedgraph$') > 0).all(axis=1))]
if species in high_rep_species:  
    # High nb of TGIRT-Seq samples (>=3 conditions), thus at least one avg 
    # condition should be > 10 (and later, the 3rd filter is overall avg should be > 5)
    _names = pd.unique([b.split('.')[0].rsplit('_', maxsplit=2)[0] for b in bedgraph_names if '_fwd' in b])
    for n in _names:
        pivot_df[n+'_fwd__AVG'] = pivot_df.filter(
                            regex=n+'_[0-9]*(.bam)?_fwd.bedgraph$').mean(axis=1)
        pivot_df[n+'_rev__AVG'] = pivot_df.filter(
                            regex=n+'_[0-9]*(.bam)?_rev.bedgraph$').mean(axis=1)
    pivot_df = pivot_df[(pivot_df.filter(regex='__AVG$') > 10).any(axis=1)]
    pivot_df = pivot_df.drop(columns=[c for c in pivot_df.columns if '__AVG' in c])

# Avg coverage across all samples/bedgraphs in which the prediction is expressed > 0
bedg_cols = pivot_df.filter(regex='bedgraph$')
pivot_df['avg_coverage_samples'] = bedg_cols.where(bedg_cols > 0).mean(axis=1)
pivot_df = pivot_df.rename(columns=lambda name: name.split(
                    '.bedgraph')[0] if name.endswith('.bedgraph') else name)

# Concat missing predictions (not detected in any condition)
not_expressed_pred = preds[~preds['gene_id'].isin(pivot_df.gene_id)]
not_expressed_pred = not_expressed_pred[['gene_id', 'chr', 'start', 
                                        'end', 'strand']]
new_cols = [n.split(
    '.bedgraph')[0] for n in simple_bedgraph_names_pos + simple_bedgraph_names_neg] + ['avg_coverage_samples']
not_expressed_pred[new_cols] = 0  # give 0 for coverage value for these 
                                  # non-detected predictions
    
all_preds = pd.concat([pivot_df, not_expressed_pred]).sort_values(
                                            ['chr', 'strand', 'start', 'end'])



""" Apply 3rd filter: Avg across conditions must be at least >5 reads"""
# Add filter to consider if prediction is expressed and a potential snoRNA
exp_preds = all_preds[all_preds['avg_coverage_samples'] >= 5]
exp_preds['dot'] = '.'
exp_preds_bed = BedTool.from_dataframe(exp_preds[['chr', 'start', 'end', 'gene_id', 
            'dot', 'strand', 'avg_coverage_samples']].sort_values(
            by=['chr', 'start']).reset_index(drop=True))

exp_preds[exp_preds['strand'] == '+'][['chr', 'start', 'end', 'gene_id', 
        'dot', 'strand', 'avg_coverage_samples']].to_csv('exp_preds_plus.tsv', sep='\t', header=False, index=False)
exp_preds[exp_preds['strand'] == '-'][['chr', 'start', 'end', 'gene_id', 
        'dot', 'strand', 'avg_coverage_samples']].to_csv('exp_preds_minus.tsv', sep='\t', header=False, index=False)
sp.call('sort -V exp_preds_plus.tsv > exp_preds_plus_sorted.tsv', shell=True)
sp.call('sort -V exp_preds_minus.tsv > exp_preds_minus_sorted.tsv', shell=True)

exp_preds_bed_pos = BedTool('exp_preds_plus_sorted.tsv')
exp_preds_bed_neg = BedTool('exp_preds_minus_sorted.tsv')


# Find overlap between annotated C/D and expressed predictions
pos_overlap = bed_cd.intersect(exp_preds_bed, s=True, wb=True, 
                f=0.5)
final_pos = pos_overlap.to_dataframe(names=['chr', 'start_cd', 'end_cd', 
                'gene_id', 'score', 'strand', 'sno_type', 'chr2', 'start_pred', 
                'end_pred', 'prediction_id', 'score2', 'strand2', 'avg_coverage_samples'])
final_pos = final_pos[['chr', 'start_cd', 'end_cd', 'strand', 'gene_id', 
                'start_pred', 'end_pred', 'prediction_id', 'avg_coverage_samples']]
final_pos.to_csv(snakemake.output.known_cd_df, sep='\t', index=False)

# Create IGV script to automate screenshots at all predictions that overlap with known expressed CD
# extend by 10 nt the screenshot window on each side
pos_overlap_extended = pos_overlap.slop(b=10, g=chr_size)
pos_overlap_extended.igv(name=True, path=screenshot_dir+"known_cd", img="png")
sp.call(f'mv {pos_overlap_extended.igv_script} {snakemake.output.batch_script_known_cd}', shell=True)
pos_overlap_extended = pos_overlap_extended.to_dataframe(names=['chr', 'start', 'end', 
                'known_cd_id', 'score', 'strand', 'sno_type', 'chr_pred', 
                'start_pred', 'end_pred', 'gene_id', 'score+pred', 
                'strand_pred', 'avg_coverage'])
num_pos = len(pos_overlap_extended)


# Add additional filter for potential new C/D predictions
""" Filter out predictions which have read coverage that extend outside of the 
    prediction for 20 nt (flanking coverage > 1/3 pred coverage)"""
not_known_cd = exp_preds[~exp_preds['gene_id'].isin(pos_overlap_extended['gene_id'])]

flanking_pos = exp_preds_bed_pos.flank(b=20, g=chr_size)
flanking_neg = exp_preds_bed_neg.flank(b=20, g=chr_size)
if species == 'homo_sapiens':
    flanking_pos.to_dataframe().to_csv('flanking_pos.tsv', sep='\t', index=False, header=False)
    sp.call('sort -k1,1V -k2,2n flanking_pos.tsv > flanking_pos_sorted.tsv', shell=True)
    flanking_pos = BedTool('flanking_pos_sorted.tsv')
    flanking_neg.to_dataframe().to_csv('flanking_neg.tsv', sep='\t', index=False, header=False)
    sp.call('sort -k1,1V -k2,2n flanking_neg.tsv > flanking_neg_sorted.tsv', shell=True)
    flanking_neg = BedTool('flanking_neg_sorted.tsv')
    cov_flanking_pos = flanking_pos.intersect(b=bedgraphs_sorted_pos, wa=True, wb=True, sorted=True,
                        names=simple_bedgraph_names_pos, g=chr_size)
    cov_flanking_neg = flanking_neg.intersect(b=bedgraphs_sorted_neg, wa=True, wb=True, sorted=True,
                        names=simple_bedgraph_names_neg, g=chr_size)
else:
    cov_flanking_pos = flanking_pos.intersect(b=bedgraphs_sorted_pos, wa=True, wb=True,
                            names=simple_bedgraph_names_pos)
    cov_flanking_neg = flanking_neg.intersect(b=bedgraphs_sorted_neg, wa=True, wb=True,
                            names=simple_bedgraph_names_neg)
col_names = ['chr', 'start', 'end', 'gene_id', 'score', 'strand', 'fake_cov', 
            'bg_name', 'chr_bg', 'start_bg', 'end_bg', 'coverage']
df_cov_flanking_pos = cov_flanking_pos.to_dataframe(names=col_names)
df_cov_flanking_neg = cov_flanking_neg.to_dataframe(names=col_names)
df_cov_flanking = pd.concat([df_cov_flanking_pos, df_cov_flanking_neg]).sort_values(
                    ['gene_id', 'bg_name']).reset_index(drop=True)

# Predictions that are expressed and that do not have any accumlation flanking the prediction (so potential real C/D)
no_flanking_but_expressed = not_known_cd[~not_known_cd['gene_id'].isin(df_cov_flanking.gene_id)]

# Correct start/end of bedgraph so that we keep only fully overlapping regions 
# with the predictions (not extended after the pred)
df_cov_flanking.loc[df_cov_flanking['end_bg'] > df_cov_flanking['end'], 'end_bg'] = df_cov_flanking['end']
df_cov_flanking.loc[df_cov_flanking['start_bg'] < df_cov_flanking['start'], 'start_bg'] = df_cov_flanking['start']

# Find per prediction what is the average coverage over the flanking sequence 
df_cov_flanking['cov_length'] = df_cov_flanking['end_bg'] - df_cov_flanking['start_bg']
df_cov_flanking['weighted_coverage'] = df_cov_flanking['cov_length'] * df_cov_flanking['coverage']


# find the sum of bedgraph coverage that intersects with the prediction
groupby_cols = ['gene_id', 'chr', 'start', 'end', 'strand', 'bg_name']
grouped_df_cov_flank = df_cov_flanking.groupby(groupby_cols).agg(total_weighted_coverage=(
    'weighted_coverage', 'sum'), total_cov_length=('cov_length', 'sum')
    ).reset_index()


grouped_df_cov_flank['avg_coverage'] = grouped_df_cov_flank[
                            'total_weighted_coverage'] / 20  # 20 nt flanking region



# Pivot df to create 1 avg_coverage col per bedgraph
pivot_df_cov_flank = grouped_df_cov_flank.pivot_table(index=[
                        'gene_id', 'chr', 'start', 'end', 'strand'], 
                        columns='bg_name', values='avg_coverage', 
                        fill_value=0).reset_index()

# Avg coverage flanking across all samples/bedgraphs in which the flanking region is expressed > 0
bedgra_cols = pivot_df_cov_flank.filter(regex='bedgraph$')
pivot_df_cov_flank['avg_coverage_samples'] = bedgra_cols.where(bedgra_cols > 0).mean(axis=1)

# Filter predictions so that the average coverage of the flanking regions must 
# be max at a third of the avg coverage of the prediction
cov_dict = dict(zip(exp_preds.gene_id, exp_preds.avg_coverage_samples))
pivot_df_cov_flank['avg_cov_pred'] = pivot_df_cov_flank['gene_id'].map(cov_dict)
fake_extended_sno_id = list(pivot_df_cov_flank[pivot_df_cov_flank['avg_coverage_samples'] > pivot_df_cov_flank['avg_cov_pred'] /3]['gene_id'])
pivot_df_cov_flank = pivot_df_cov_flank[pivot_df_cov_flank['avg_coverage_samples'] <= pivot_df_cov_flank['avg_cov_pred'] / 3]
pivot_df_cov_flank = pivot_df_cov_flank[~pivot_df_cov_flank['gene_id'].isin(fake_extended_sno_id)] # to filter out preds if only one of left/right extension extends past the sno 


"""Filter to exclude the predictions overlapping known cd"""
potential_new_cd = pivot_df_cov_flank[~pivot_df_cov_flank['gene_id'].isin(pos_overlap_extended['gene_id'])].reset_index(drop=True)

# Concat predictions that are expressed w no flanking read accumulation to 
# those expressed with small read levels on their flanking regions (potential real candidates!)
cand_list = list(no_flanking_but_expressed.gene_id) + list(potential_new_cd.gene_id)
final_candidates = exp_preds[exp_preds['gene_id'].isin(cand_list)].sort_values('avg_coverage_samples').reset_index(drop=True)
final_candidates.drop(columns=['dot']).to_csv(snakemake.output.new_cd_df, sep='\t', index=False)
print(final_candidates)
bed_final = BedTool.from_dataframe(final_candidates[['chr', 'start', 'end', 'gene_id', 'dot', 'strand']])
bed_final_ext = bed_final.slop(b=10, g=chr_size)
bed_final_ext.igv(name=True, path=screenshot_dir+"new_cd", img="png")
sp.call(f'mv {bed_final_ext.igv_script} {snakemake.output.batch_script_new_cd}', shell=True)


print(species)
print(f'Total number of SnoBIRD predictions: {len(all_preds)}')
print('Total number of expressed SnoBIRD predictions (avg coverage >=5 and flanking reads filters): '+
        f'{num_pos + len(final_candidates)},')
print(f'with {num_pos} being already annotated and {len(final_candidates)} potential novel C/D snoRNAs')
print(f'Total number of known annotated C/D snoRNAs: {len(cd_df_pos) + len(cd_df_neg)}')

# Info df
cols = ['species', 'total_annot_cd', 'total_annot_exp', 'total_snoBIRD_preds', 
        'total_annot_preds', 'exp_snoBIRD_preds_5_reads_and_filters', 
        'exp_annot_preds', 'potential_novel_expressed']
final_vals = [species, len(cd_df_pos) + len(cd_df_neg), num_annot_expressed, 
                len(all_preds), num_annotated_and_pred, 
                num_pos+len(final_candidates), num_pos, len(final_candidates)]
temp_df = pd.DataFrame([final_vals], columns=cols)
print(temp_df)

temp_df.to_csv(snakemake.output.coverage, sep='\t', index=False)



sp.call(f'rm {species}_ensembl.bed cd_test_* preds_neg* preds_pos*', shell=True)
sp.call(f'rm results/bedgraph_TGIRT/{species}/final_* exp_preds_* flanking_pos* flanking_neg*', shell=True)
