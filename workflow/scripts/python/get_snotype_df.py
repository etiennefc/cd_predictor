#!/usr/bin/python3
import pandas as pd
from pybedtools import BedTool
import subprocess as sp

# Load dfs and inputs
gtf = snakemake.input.gtf
species = snakemake.wildcards.species
rnacentral = pd.read_csv(snakemake.input.rnacentral, sep='\t')
rnacentral['dot'] = '.'
rnacentral[['chr', 'start', 'end', 'gene_id', 'dot', 'strand']].to_csv(
    f'rnacentral_{species}.bed', sep='\t', index=False, header=None)

# Extract C/D snoRNAs from gtf
if species in ['arabidopsis_thaliana', 'drosophila_melanogaster', 'danio_rerio']:
    cd_symbols = 'SNOR108|SNOR105|SNOR30|SNOR37-2|U49.1|U460.[12]F|snoRNA:Me|scaRNA:Me|snoRNA:185|' + \
        'snoRNA:229|snoRNA:291|snoRNA:684|snoRNA:83E4-5|snoRNA:CD-24|snoRNA:kis-a|snoRNA:nop5-x16-a|' + \
        'snoRNA:Or-CD|snoRNA:snR38:54E[bc]|snoRNA:U[12347]'  # CD info taken manually from Ensembl and Flybase
    sp.call(
    f'''awk -v OFS="\t" 'NR>5 && $3== "gene" && /gene_biotype "snoRNA"/' {gtf} | ''' + \
    f'''grep -E "{cd_symbols}" | awk -v OFS="\t" '{{print $1,$4,$5,$10,$6,$7}}' | '''+ \
    f'''sed 's/;//g; s/"//g' > gtf_cd_{species}.bed''', shell=True)
else:  # P. falciparum does not have the biotype snoRNA, but has some snoRNAs annotated as snRNAs or ncRNAs
    # So we need to exlucde everything that is not related to C/D, plus exclude 
    # a couple of H/ACA based on manual inspection of their sequence
    non_cd_symbols = 'RUF|ACA|U[45612]|ase|signal|unknown|gene_name "ncRNA"|snoR26|snoR14"|snoR31|'+ \
                    'snoR27|snoR04|snoR03|snoR35|snoR36'
    sp.call(f'''awk -v OFS="\t" 'NR>5 && $3== "gene" && /gene_biotype "[nc]*[sn]*RNA"/' {gtf} | '''+ \
        f'''grep -vE '{non_cd_symbols}' | awk -v OFS="\t" '{{print $1,$4,$5,$10,$6,$7}}' | '''+ \
        f'''sed 's/;//g; s/"//g' > gtf_cd_{species}.bed''', shell=True)


# Concat and then Bedtools merge RNAcentral entries and those from the gtf
sp.call(f'cat rnacentral_{species}.bed gtf_cd_{species}.bed | sort -k1,1 -k2,2n > {species}_all_cd.bed', shell=True)
bed_cd = BedTool(f'{species}_all_cd.bed')
merged_cd = bed_cd.merge(s=True, c=[4, 6], o=['distinct', 'distinct'])
merged_cd = merged_cd.to_dataframe(names=['chr', 'start','end', 'gene_id', 'strand'])
merged_cd['sno_type'] = 'C/D'
merged_cd['gene_id'] = merged_cd['gene_id'].str.replace(',', '_')
merged_cd['chr'] = 'chr' + merged_cd['chr'].astype(str)
print(species)
print(merged_cd)
merged_cd.to_csv(snakemake.output.snotype_df, sep='\t', index=False)


sp.call(f'rm rnacentral_{species}.bed gtf_cd_{species}.bed {species}_all_cd.bed', shell=True)

