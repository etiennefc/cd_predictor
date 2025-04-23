#!/usr/bin/python3
import pandas as pd
import subprocess as sp
from pybedtools import BedTool
import re

# Load inputs and outputs
snoscan_pred_dirs = snakemake.input.snoscan_preds
snoreport_pred_dirs = snakemake.input.snoreport_preds
infernal_pred_dirs = snakemake.input.infernal_preds
human_chr_df = pd.read_csv([p for p in snakemake.input.chr_size if 'homo_sapiens' in p][0], 
                    names=['chr', 'size'], sep='\t')
pombe_chr_df = pd.read_csv([p for p in snakemake.input.chr_size if 'schizosaccharomyces_pombe' in p][0], 
                    names=['chr', 'size'], sep='\t')
human_chr_dict = dict(zip(human_chr_df['chr'], human_chr_df['size']))
pombe_chr_dict = dict(zip(pombe_chr_df['chr'], pombe_chr_df['size']))


snoscan_output = [snakemake.output.snoscan_pombe, snakemake.output.snoscan_human]
snoreport_output = [snakemake.output.snoreport_pombe, snakemake.output.snoreport_human]
infernal_output = [snakemake.output.infernal_pombe, snakemake.output.infernal_human]

# Filter snoscan output
def filter_snoscan(fasta_dir, species):
    snoscan_preds, snoscan_bed = [], []
    sp.call(f'cat {fasta_dir}/*txt >> snoscan_all_{species}.fa', shell=True)
    with open(f'snoscan_all_{species}.fa', 'r') as f:
        for line in f:
            if line.startswith('>>'):
                chrom = line.split(' ')[1]
                if chrom.startswith('chr') == False:
                    chrom = 'chr' + chrom
                
                # get start, end and strand
                s, e = re.findall('\([0-9]+-[0-9]+\)', line)[0].strip('()').split('-')
                if int(s) > int(e):
                    strand = '-'
                    start, end = e, s
                else:
                    strand = '+'
                    start, end = s, e
                
                # if location not already present in snoscan_preds, add it 
                # (to filter out same pred but with different rRNA target)
                gene_id = f'snoscan_{chrom}_{start}_{end}_{strand}'
                if gene_id not in snoscan_preds:
                    snoscan_preds.append(gene_id)
                    snoscan_bed.append([chrom, start, end, gene_id, '.', strand])
    
    # Create df then save as bed
    df_snoscan = pd.DataFrame(snoscan_bed, columns=[
                    'chr', 'start', 'end', 'gene_id', 'score', 'strand']).sort_values(['chr', 'start'])
    df_snoscan[['start', 'end']] = df_snoscan[['start', 'end']].astype(int)
    print(df_snoscan)
    df_snoscan.to_csv([p for p in snoscan_output if species in p][0], sep='\t', 
                        index=False, header=False)

    sp.call(f'rm snoscan_all_{species}.fa', shell=True)




# Filter snoreport output
def filter_snoreport(fasta_dir, species, dictio):
    snoreport_bed = []
    if species == 'schizosaccharomyces_pombe':
        sp.call(f'cat {fasta_dir}/*fa >> snoreport_all_{species}.fa', shell=True)
    elif species == 'homo_sapiens':
        print('ok')
        sp.call(f'''for i in {fasta_dir}/*fa; do p=$(head -n1 $i | grep -Eow '>[^_]*_1'); '''+ \
            f'''awk -v pat=$p 'BEGIN {{count = 0; output = "snoreport_{species}_plus"}}$0 ~ pat" " {{count++}} count == 2 {{output = "snoreport_{species}_minus"}} {{print >> output}}' $i; done''', shell=True)
        sp.call(f'cat snoreport_{species}_plus snoreport_{species}_minus > snoreport_all_{species}.fa', shell=True)
    with open(f'snoreport_all_{species}.fa', 'r') as f:
        pos_gene_id, prev_strand = [], '+'
        for line in f:
            if line.startswith('>'):
                gene_id = line.split(' ')[0].strip('>')
                chrom = 'chr' + gene_id.split('_')[0]
                s = int(line.split(' ')[3])
                e = int(line.split(' ')[4])
                # To account for predictions on strand - appended after strand +
                if (gene_id in pos_gene_id) | (prev_strand == '-'): 
                    #print('NEG')
                    strand = '-'
                    # Adjust position based on chr size
                    start = dictio[chrom] - e
                    end = dictio[chrom] - s
                else:
                    strand = '+'
                    start = s
                    end = e
                real_gene_id = f'snoreport_{chrom}_{start}_{end}_{strand}'

                snoreport_bed.append([chrom, start, end, real_gene_id, '.', strand])
                pos_gene_id.append(gene_id)
                prev_strand = strand
    
    # Create df then save as bed
    df_snoreport = pd.DataFrame(snoreport_bed, columns=[
                    'chr', 'start', 'end', 'gene_id', 'score', 'strand']).sort_values(['chr', 'start'])
    df_snoreport[['start', 'end']] = df_snoreport[['start', 'end']].astype(int)
    print(df_snoreport)
    df_snoreport.to_csv([p for p in snoreport_output if species in p][0], sep='\t', 
                        index=False, header=False)

    sp.call(f'rm snoreport_all_{species}.fa snoreport_{species}_*', shell=True)



# Filter infernal output
def filter_infernal(fasta_dir, species):
    infernal_bed, infernal_id = [], []
    if species == 'schizosaccharomyces_pombe':
        # remove unwanted H/ACA snoRNAs predicted as snoRNAs by infernal (not to inflate TP)
        pombe_specific_haca = 'RF0143[456789]|RF0144[023456789]|RF0145[012]' 
        sp.call(f'cat {fasta_dir}/*tblout | grep -vE "{pombe_specific_haca}" >> infernal_all_{species}.tblout', shell=True)
    else:
        sp.call(f'cat {fasta_dir}/*tblout >> infernal_all_{species}.tblout', shell=True)
    with open(f'infernal_all_{species}.tblout', 'r') as f:
        for line in f:
            if ('#' and 'SNORA' not in line) & ('mall nucleolar RNA' in line):
                gene_id = re.split('RF[0-9]{5}   ', line)[1].split('-')[0].strip(' ')
                chrom = 'chr' + gene_id
                # get start, end and strand
                location = re.findall('([1-9][0-9]+[\t\s][\t\s]*[1-9][0-9]+)[\t\s]+[-\+]', line)[-1]
                s, e = location.split(' ', maxsplit=1)
                s, e = int(s.strip(' ')), int(e.strip(' '))

                if int(s) > int(e):
                    strand = '-'
                    start, end = e, s
                else:
                    strand = '+'
                    start, end = s, e
                real_gene_id = gene_id = f'infernal_{chrom}_{start}_{end}_{strand}'

                if real_gene_id not in infernal_id:
                    infernal_bed.append([chrom, start, end, real_gene_id, '.', strand])
                    infernal_id.append(real_gene_id)
    
    # Create df then save as bed
    df_infernal = pd.DataFrame(infernal_bed, columns=[
                    'chr', 'start', 'end', 'gene_id', 'score', 'strand']).sort_values(['chr', 'start'])
    df_infernal[['start', 'end']] = df_infernal[['start', 'end']].astype(int)
    print(df_infernal)
    df_infernal.to_csv([p for p in infernal_output if species in p][0], sep='\t', 
                        index=False, header=False)

    sp.call(f'rm infernal_all_{species}.tblout', shell=True)




# Filter results of the different tools
filter_snoscan([i for i in snoscan_pred_dirs if 'schizosaccharomyces_pombe' in i][0], 
                'schizosaccharomyces_pombe')
filter_snoscan([i for i in snoscan_pred_dirs if 'homo_sapiens' in i][0], 'homo_sapiens')


filter_snoreport([i for i in snoreport_pred_dirs if 'schizosaccharomyces_pombe' in i][0], 
                'schizosaccharomyces_pombe', pombe_chr_dict)
filter_snoreport([i for i in snoreport_pred_dirs if 'homo_sapiens' in i][0], 
                'homo_sapiens', human_chr_dict)


filter_infernal([i for i in infernal_pred_dirs if 'homo_sapiens' in i][0], 
                'homo_sapiens')
filter_infernal([i for i in infernal_pred_dirs if 'schizosaccharomyces_pombe' in i][0], 
                'schizosaccharomyces_pombe')

