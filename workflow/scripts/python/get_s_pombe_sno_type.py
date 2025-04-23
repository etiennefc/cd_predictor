#!/usr/bin/python3
import pandas as pd
import subprocess as sp

url = snakemake.params.pombase_url
output = snakemake.output.snotype
gtf = pd.read_csv(snakemake.input.gtf, sep='\t', skiprows=5, names=[
        'chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 
        'attributes'])
gtf = gtf[gtf['feature'] == 'gene']

# Keep snoRNA entries and clean columns values
snoRNA_gtf = gtf[gtf['attributes'].str.contains('gene_biotype "snoRNA"')]
snoRNA_gtf[['gene_id', 'gene_name', 'gene_source', 'gene_biotype', 'empty']] = snoRNA_gtf[
    'attributes'].str.split(';', expand=True)
snoRNA_gtf['gene_id'] = snoRNA_gtf['gene_id'].str.replace(
                    'gene_id "', '').str.replace('"', '').str.strip()
snoRNA_gtf['gene_name'] = snoRNA_gtf['gene_name'].str.replace(
                    'gene_name "', '').str.replace('"', '').str.strip()


# Filter to keep all PomBase snoRNAs and two Ensembl snoRNAs that are more 
# plausible annotation than the two respective snoRANs annotated in PomBase
snoRNA_gtf = snoRNA_gtf[(~snoRNA_gtf['gene_name'].isin(['snoU24', 'snoR61'])) & (
    (snoRNA_gtf['source'] == 'PomBase') | (snoRNA_gtf['gene_name'].isin([
        "snosnR61", "SNORD24"])))].sort_values(by=['chr', 'strand', 'start'])
snoRNA_gtf = snoRNA_gtf[
    ['gene_id', 'gene_name', 'source', 'chr', 'start', 'end', 'strand']].reset_index(
        drop=True)

haca = ['H/ACA', 'pseudour', 'SNORA', 'Psi']
cd = ['C/D', 'methyl', 'SNORD', 'U3', 'U8', 'Nm']
# Based on Rfam family belonging or on the boxes in the sequence, one can 
# recognize C/D and H/ACA for the snoRNAs that miss this information on PomBase
known_cd = ['SPSNORNA.01', 'SPSNORNA.07', 'SPSNORNA.08', 'SPSNORNA.09', 
            'SPSNORNA.10', 'SPSNORNA.11', 'SPSNORNA.31', 'SPSNORNA.54', 
            'SPSNORNA.13', 'SPSNORNA.15', 'SPSNORNA.19', 'SPSNORNA.12', 
            'SPSNORNA.30', 'SPSNORNA.16', 'SPSNORNA.17', 'SPSNORNA.21', 
            'SPSNORNA.22', 'SPSNORNA.23', 'SPSNORNA.24', 'SPSNORNA.28', 
            'SPSNORNA.26', 'SPSNORNA.27', 'ENSRNA049676952', 'ENSRNA049678266']
type_dict = {}
for gene_id in snoRNA_gtf['gene_id']:
    sp.call(f'curl -s -o {gene_id} {url}{gene_id}', shell=True)
    cd_num, haca_num = 0, 0
    with open(gene_id, 'r') as f:
        for line in f:
            if any(s in line for s in cd):
                cd_num += 1
            if any(s in line for s in haca):
                haca_num += 1
        if cd_num == 0 and haca_num > 0:
            type_dict[gene_id] = 'H/ACA'
        elif cd_num > 0 and haca_num == 0:
            type_dict[gene_id] = 'C/D'
        else:
            if gene_id in known_cd:
                type_dict[gene_id] = 'C/D'
            else:
                type_dict[gene_id] = 'Unknown'
    sp.call(f'rm {gene_id}', shell=True)

snoRNA_gtf['sno_type'] = snoRNA_gtf['gene_id'].map(type_dict)
snoRNA_gtf.to_csv(output, sep='\t', index=False)


