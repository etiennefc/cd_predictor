#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft 
import seaborn as sns 
from glob import glob
import warnings
warnings.filterwarnings("ignore")
from Bio import SeqIO
import subprocess as sp
from pybedtools import BedTool
import matplotlib.patches as mpatches
from scipy.stats import pearsonr, mannwhitneyu
import numpy as np 

# Load gtf and chr size and params
colors = snakemake.params.color_dict
colors = dict(sorted(colors.items()))
gtfs = snakemake.input.gtf
chr_sizes = snakemake.input.chr_size

# Load total number of snoBIRD preds per species
preds = {}
species_all = []
for p in snakemake.input.snoBIRD_preds:
    spe = p.split('/')[4]
    preds_df = pd.read_csv(glob(p+'/*tsv')[0], sep='\t')
    preds[spe] = len(preds_df)
    species_all.append(spe)
species_all = sorted(species_all)
print('Total SnoBIRD preds per genome:\n', preds)

## Find total genome size and GC content per genome
genome_size, gc_content = {}, {}
for p in snakemake.input.genome:
    g_size, gc_count = 0, 0
    spe = p.split('/')[-1].replace('_genome.fa', '')
    # Add up chr size
    for record in SeqIO.parse(p, 'fasta'):
        g_size += len(record.seq)
        gc_count += record.seq.count('C') + record.seq.count('G')
    genome_size[spe] = g_size
    gc_content[spe] = (gc_count/g_size) * 100

print('GC content per genome:\n', gc_content)
print('Size of genomes:\n', genome_size)


# Find intron content in the genomes (% of intron, relative nb of introns and intron length distribution)
# To do so, find the complement of all exons_bed, and then substract intergenic_regions bed
def get_intergenic_regions(gtf, species_name, chr_size):
    """ Get intergenic regions in a genome by taking
        the bedtools complement of all genes (i.e. all the 
        regions not included in these genes)."""
    # Create bed of all genes
    chr_ = 'chr'
    if species_name == 'tetrahymena_thermophila':
        chr_ = ''
    sp.call(f'''awk -v OFS="\t" '$3=="gene" {{print "{chr_}"$1,$4,$5,$10,$6,$7}}' {gtf} | '''+ \
        f'''grep -vE 'chrKZ|chrKN' | sort -k1,1 -k2,2n > all_genes_{species_name}_sorted.bed''', shell=True)

    gene_bed = BedTool(f'all_genes_{species_name}_sorted.bed')

    ## Create temporary sorted chr size file
    sp.call(f'sort -k1,1 -k2,2n {chr_size} > temp_chr_size{species_name}.tsv', shell=True)

    # Get the complement of these genes
    complement = gene_bed.complement(g=f'temp_chr_size{species_name}.tsv')


    return complement


def get_intronic_regions(gtf, species_name, chr_size):
    """ Get the intronic regions by substracting intergenic regions to the 
        complement of exons_regions. Exclude retained_intron as exons to keep 
        more 'real' introns. """
    # Get intergenic regions
    intergenic_bed = get_intergenic_regions(gtf, species_name, chr_size)

    # Select intronic and intergenic regions using the exons of all genes as boundaries
    chr_ = 'chr'
    if species_name == 'tetrahymena_thermophila':
        chr_ = ''
    sp.call(f'''awk -v OFS="\t" '!/retained_intron/ && $3=="exon" {{print "{chr_}"$1,$4,$5,$10,$6,$7}}' {gtf} | '''+ \
    f'''grep -vE 'chrKZ|chrKN' | sort -k1,1 -k2,2n > all_exons_{species_name}_sorted.bed''', shell=True)
    gene_bed = BedTool(f'all_exons_{species_name}_sorted.bed')

    # Get the complement of these exons (i.e intron and intergenic regions)
    complement = gene_bed.complement(g=f'temp_chr_size{species_name}.tsv')

    # Subtract intergenic regions to the complement to get intronic regions only
    intronic_regions = complement.subtract(b=intergenic_bed)

    introns_df = intronic_regions.to_dataframe()
    introns_df['intron_length'] = introns_df.end - introns_df.start + 1
    introns_df['species_name'] = species_name

    sp.call(f'rm all_genes_{species_name}*.bed temp_chr_size{species_name}.tsv', shell=True)
    sp.call(f'rm all_exons_{species_name}*.bed', shell=True)
    
    return introns_df

intron_dfs = []
intron_nb = {}
rel_intron_nb = {}  # nb of intron per 100,000 nt
intron_percent = {}
danio, gallus = [], []
for gtf in gtfs:
    spe = gtf.split('/')[-1].replace('.gtf','')
    c_size = [i for i in chr_sizes if spe in i][0]
    introns = get_intronic_regions(gtf, spe, c_size)
    intron_dfs.append(introns)
    # total nb of introns
    intron_nb[spe] = len(introns)
    # total nb of intron per 100,000 nt
    rel_intron_nb[spe] = (len(introns) / genome_size[spe]) * 100000
    # % of intron in overall genome sequence
    total_intron_len = introns['intron_length'].sum()
    if spe == 'danio_rerio':
        danio.append(introns)
    elif spe == 'gallus_gallus':
        gallus.append(introns)
    intron_percent[spe] = (total_intron_len / genome_size[spe]) * 100
    
print('Relative intron nb / 100kb: \n', rel_intron_nb)
print('Total intron nb:\n', intron_nb)
print('Intron proportion in genome:\n', intron_percent)

# Create multiplot:
# density plot of intron length distribution and scatter plots showing the 
# relationship between the total nb of SnoBIRD predictions and other features 
# (genome size, GC content, intron nb, relative intron nb and intron percentage)
plt.rcParams['svg.fonttype'] = 'none'
dictio = [None, genome_size, gc_content, intron_nb, rel_intron_nb, intron_percent]
col_name = ['intron_length_log10', 'genome_size', 'GC_content', 'total_intron_nb', 'intron_nb_per_100kb', 'intron_percent']
fig, ax = plt.subplots(2, 3, figsize=(24, 12))
for j, ax in enumerate(ax.flat):
    if j == 0:
        for i, df in enumerate(intron_dfs):
            spe = list(pd.unique(df.species_name))[0]
            df['intron_length_log10'] = np.log10(df['intron_length'])
            sns.kdeplot(data=df[col_name[j]], fill=True, ax=ax, color=colors[spe])
        
        # Compute Mann-Whitney test between D. rerio and G. gallus distributions
        U1, p_mn = mannwhitneyu(list(danio[0]['intron_length']), list(gallus[0]['intron_length']))
        median_danio, median_gallus = danio[0]['intron_length'].median(), gallus[0]['intron_length'].median()
        ax.text(3, 4, 
            f"Mann-Whitney pval={p_mn:.2e}\nD. rerio median={median_danio:.2f}\nG.gallus median={median_gallus:.2f}", 
            fontsize=12)
        ax.set_xlabel('Intron length (log10(nt))', fontdict={'fontsize': 25})
        ax.set_ylabel('Density', fontdict={'fontsize': 25})
        #ax.set_xscale('log')
        ax.spines['right'].set_linewidth(0)
        ax.spines['top'].set_linewidth(0)
        legend_list = []
        for crit in species_all:
            legend_element = mpatches.Patch(color=colors[crit], label=crit)
            legend_list.append(legend_element)
        plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(-1.4,2),
                    fontsize=12)
    else:
        temp_df = pd.DataFrame(dictio[j], index=[0]).T.reset_index()
        temp_df.columns = ['species_name', col_name[j]]
        temp_df['SnoBIRD_total_preds'] = temp_df['species_name'].map(preds)
        temp_df = temp_df.sort_values('species_name').reset_index(drop=True)
        r_val, p_val = pearsonr(temp_df['SnoBIRD_total_preds'], temp_df[col_name[j]])
        sns.regplot(data=temp_df, x='SnoBIRD_total_preds', y=col_name[j], color='grey', ax=ax)
        sns.scatterplot(temp_df, x='SnoBIRD_total_preds', y=col_name[j], 
                hue='species_name', palette=colors, ax=ax, s=100, legend=False)
        ax.set_xlabel('Total number of SnoBIRD predictions', fontdict={'fontsize': 25})
        ax.set_ylabel(col_name[j], fontdict={'fontsize': 25})
        ax.text(25000, temp_df[col_name[j]].max(), f"Pearson_r={r_val:.2f}\npval={p_val:.2e}", fontsize=12)

plt.savefig(snakemake.output.multi_plot, bbox_inches='tight', dpi=500)

