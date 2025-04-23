#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import functions as ft 
from glob import glob
import os
import warnings
warnings.filterwarnings("ignore")
import subprocess as sp
from pybedtools import BedTool


# Load gtf and chr size and params
colors = snakemake.params.color_dict
gtfs = snakemake.input.gtf
chr_sizes = snakemake.input.chr_size
snoBIRD_pred_dirs = snakemake.input.snoBIRD_preds
snoscan_preds = snakemake.input.snoscan_preds
snoreport_preds = snakemake.input.snoreport_preds
infernal_preds = snakemake.input.infernal_preds
all_spe = ['schizosaccharomyces_pombe', 'plasmodium_falciparum', 'tetrahymena_thermophila', 
            'arabidopsis_thaliana', 'drosophila_melanogaster', 'gallus_gallus', 'danio_rerio', 
            'macaca_mulatta',  'homo_sapiens']

# Get SnoBIRD predictions as bed
total_snoBIRD_preds = {}
snoBIRD_beds = {}
for p in snoBIRD_pred_dirs:
    spe = p.split('/')[4]
    preds_df = pd.read_csv(glob(p+'/*tsv')[0], sep='\t')
    chr_ = list(pd.unique(preds_df['chr']))
    if str(chr_[0]).startswith('chr') == False:
        preds_df['chr'] = 'chr' + preds_df['chr'].astype(str)
    preds_df = preds_df[
        ['chr', 'start', 'end', 'gene_id', 'probability_CD', 'strand']
        ]
    total_snoBIRD_preds[spe] = len(preds_df)
    preds_df.to_csv(f'snoBIRD_{spe}.bed', index=False, header=False, sep='\t')
    sp.call(f'sort -V snoBIRD_{spe}.bed > snoBIRD_{spe}_sorted.bed', shell=True)
    snoBIRD_beds[spe] = BedTool(f'snoBIRD_{spe}_sorted.bed')

# Get snoscan predictions as bed
total_snoscan_preds = {}
snoscan_beds = {}
for p in snoscan_preds:
    spe = p.split('snoscan_')[1].replace('_filtered.tsv', '')
    sp.call(f'sort -V {p} > snoscan_{spe}_sorted.bed', shell=True)
    snoscan_beds[spe] = BedTool(f'snoscan_{spe}_sorted.bed')
    total_snoscan_preds[spe] = len(BedTool(f'snoscan_{spe}_sorted.bed').to_dataframe(
        names=['chr', 'start', 'end', 'gene_id', 'score', 'strand']))

# Get snoreport2 predictions as bed
snoreport_beds = {}
total_snoreport_preds = {}
for p in snoreport_preds:
    spe = p.split('snoreport2_')[1].replace('_filtered.tsv', '')
    sp.call(f'sort -V {p} > snoreport_{spe}_sorted.bed', shell=True)
    snoreport_beds[spe] = BedTool(f'snoreport_{spe}_sorted.bed')
    total_snoreport_preds[spe] = len(BedTool(f'snoreport_{spe}_sorted.bed').to_dataframe(
        names=['chr', 'start', 'end', 'gene_id', 'score', 'strand']))

# Get infernal predictions as bed
infernal_beds = {}
total_infernal_preds = {}
for p in infernal_preds:
    spe = p.split('infernal_rfam/infernal_rfam_')[1].replace('_filtered.tsv', '')
    sp.call(f'sort -V {p} > infernal_{spe}_sorted.bed', shell=True)
    infernal_beds[spe] = BedTool(f'infernal_{spe}_sorted.bed')
    total_infernal_preds[spe] = len(BedTool(f'infernal_{spe}_sorted.bed').to_dataframe(
        names=['chr', 'start', 'end', 'gene_id', 'score', 'strand']))



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
        f'''grep -vE 'chrKZ|chrKN|chrQNV' | sort -V > all_genes_{species_name}_sorted.bed''', shell=True)

    gene_bed = BedTool(f'all_genes_{species_name}_sorted.bed')

    ## Create temporary sorted chr size file
    sp.call(f'grep -v chrQNV {chr_size} | sort -V > temp_chr_size{species_name}.tsv', shell=True)

    # Get the complement of these genes
    complement = gene_bed.complement(g=f'temp_chr_size{species_name}.tsv')
    # Add fake gene_id, score and strand columns for downstream compatibility
    complement = complement.each(lambda x: f"{x.chrom}\t{x.start}\t{x.end}\t.\t.\t.").saveas()

    return complement


def get_genomic_regions_bed(gtf, species_name, chr_size):
    """ Get the intronic regions by substracting intergenic regions to the 
        complement of exons_regions. Exclude retained_intron as exons to keep 
        more 'real' introns. """
    # Get intergenic regions
    intergenic_bed = get_intergenic_regions(gtf, species_name, chr_size)

    # Select intronic and intergenic regions using the exons of all genes as boundaries
    chr_ = 'chr'
    if species_name == 'tetrahymena_thermophila':
        chr_ = ''
    # remove midsize ncRNA as "exons"
    mncRNA = 'Mt_tRNA|scaRNA|snoRNA|snRNA|sRNA|pre-tRNA|tRNA|miRNA|S_pombe_snR46'
    if species_name == 'plasmodium_faclciparum':
        mncRNA = mncRNA + '|small nucleolar RNA|spliceosomal RNA|snoR36'
    sp.call(f'''awk -v OFS="\t" '!/retained_intron|{mncRNA}/ && $3=="exon" {{print "{chr_}"$1,$4,$5,$10,$6,$7}}' {gtf} | '''+ \
    f'''grep -vE 'chrKZ|chrKN|chrQNV' | sort -V > all_exons_{species_name}_sorted.bed''', shell=True)
    exon_bed = BedTool(f'all_exons_{species_name}_sorted.bed')

    # Get the complement of these exons (i.e intron and intergenic regions)
    complement = exon_bed.complement(g=f'temp_chr_size{species_name}.tsv')

    # Subtract intergenic regions to the complement to get intronic regions only
    intron_bed = complement.subtract(b=intergenic_bed)
    # Add fake gene_id, score and strand columns for downstream compatibility
    intron_bed = intron_bed.each(lambda x: f"{x.chrom}\t{x.start}\t{x.end}\t.\t.\t.").saveas()
    
    return exon_bed, intron_bed, intergenic_bed


def genomic_overlap(pred_bed, gtf, species_name, chr_size):
    col_names = ['chr', 'start', 'end', 'pred_id', 'score', 'strand', 
                'genomic_region', 'chr_genomic', 'start_genomic', 'end_genomic', 
                'genomic_id','genomic_score', 'genomic_strand']
    # Get the overlap between predictions and exon, intron and intergenic regions
    if f'intergenic_{species_name}.bed' not in os.listdir('.'):
        exon_bed, intron_bed, intergenic_bed = get_genomic_regions_bed(
                                                    gtf, species_name, chr_size)
        exon_bed.to_dataframe().to_csv(f'exon_{species_name}.bed', sep='\t', index=False, header=False)
        intron_bed.to_dataframe().to_csv(f'intron_{species_name}.bed', sep='\t', index=False, header=False)
        intergenic_bed.to_dataframe().to_csv(f'intergenic_{species_name}.bed', sep='\t', index=False, header=False)
    
    # Get the proportion of the genome that is covered by each genomic element
    ex_bed = pd.read_csv(f'exon_{species_name}.bed', sep='\t', 
                        names=['chr', 'start', 'end', 'pred_id', 'score', 'strand'])
    intr_bed = pd.read_csv(f'intron_{species_name}.bed', sep='\t', 
                        names=['chr', 'start', 'end', 'pred_id', 'score', 'strand'])
    inter_bed = pd.read_csv(f'intergenic_{species_name}.bed', sep='\t', 
                        names=['chr', 'start', 'end', 'pred_id', 'score', 'strand'])
    
    total_genomic = []
    for df_ in [ex_bed, intr_bed, inter_bed]:
        df_['len'] = df_.end - df_.start + 1
        total_genomic.append(df_['len'].sum())
    
    percent_total_genomic = [i*100/sum(total_genomic) for i in total_genomic]


    
    genomic_beds = [f'exon_{species_name}.bed', f'intron_{species_name}.bed', f'intergenic_{species_name}.bed']
    genomic_names = ['Exons', 'Introns', 'intergenic'] # lowercase for intergenic so that it's the last one to be chosen
    
    # Intersect predictions with genomic regions
    overlap = pred_bed.intersect(b=genomic_beds, wa=True, wb=True, 
                            names=genomic_names, sorted=True, 
                            g=f'temp_chr_size{species_name}.tsv', f=0.5)
    overlap_df = overlap.to_dataframe(names=col_names).drop_duplicates()

    # Drop duplicate overlap (favoring exon over intron over intergenic)
    overlap_df = overlap_df.sort_values('genomic_region').drop_duplicates('pred_id').reset_index(drop=True)


    return overlap_df, percent_total_genomic, total_genomic
    



# Get per species the overlap between exon, intron and intergenic regions with the 
# predictions per tool
results_dict, percent_genomic_dict, total_genomic_dict = {}, {}, {}
tools = ['snoBIRD', 'snoscan', 'snoreport2', 'infernal_rfam']
for gtf in gtfs:
    spe = gtf.split('/')[-1].replace('.gtf','')
    c_size = [i for i in chr_sizes if spe in i][0]
    t_dict = {}
    for i, tool_dict in enumerate([snoBIRD_beds, snoscan_beds, 
                                snoreport_beds, infernal_beds]):
        if spe in tool_dict.keys():
            sp_bed = tool_dict[spe]
            df, percent_genomic, total_genomic = genomic_overlap(sp_bed, gtf, spe, c_size)
            t_dict[tools[i]] = df
            if spe not in percent_genomic_dict.keys():
                percent_genomic_dict[spe] = percent_genomic
                total_genomic_dict[spe] = total_genomic
    
    results_dict[spe] = t_dict

    sp.call(f'rm all_genes_{spe}*.bed temp_chr_size{spe}.tsv', shell=True)
    sp.call(f'rm all_exons_{spe}*.bed snoBIRD_{spe}.bed snoBIRD_{spe}_sorted.bed', shell=True)
    sp.call(f'rm exon_{spe}.bed intron_{spe}.bed intergenic_{spe}.bed', shell=True)
    sp.call(f'rm infernal_{spe}_sorted.bed snoscan_{spe}_sorted.bed snoreport_{spe}_sorted.bed', shell=True)

# Create donut chart to show the actual genome % that is exon/intron/intergenic per species
real_prop_list = [percent_genomic_dict[spe] for spe in all_spe]
#ft.pie_multiple(3, 3, donut_list, ['exonic', 'intronic', 'intergenic'], list(colors.values()), 
#                all_spe, 'Actual proportion of the genome that is exonic, intronic or intergenic across species', 
#                'Genomic element', snakemake.output.donut_actual_prop)
#ft.stacked_bar2(donut_list, all_spe, ['exonic', 'intronic', 'intergenic'], 
#    'Actual proportion of the genome that is exonic, intronic or intergenic across species', 
#    'Species',  'Proportion of genome (%)', colors, 0, 102, [''] * len(all_spe), 
#    snakemake.output.donut_actual_prop)

# For human/pombe and all tools:
# If resulting pred overlap < total pred, adjust to total by saying the rest is 
# intergenic (if they could not overlap with annotated genome locations)
models = {'snoBIRD': total_snoBIRD_preds, 'snoscan': total_snoscan_preds, 
        'snoreport2': total_snoreport_preds, 'infernal_rfam': total_infernal_preds}
final_vals = []
for spe in ['homo_sapiens', 'schizosaccharomyces_pombe']:
    vals = []
    print(spe)
    for mod in models.keys():
        print(mod)
        df = results_dict[spe][mod]
        diff_len = models[mod][spe] - len(df)
        df = pd.concat([df, pd.DataFrame([['']*6+['intergenic']+['']*6] * diff_len)])
        temp_vals = []
        for genomic in ['Exons', 'Introns', 'intergenic']:
            temp_vals.append(len(df[df['genomic_region'] == genomic]))
        print(temp_vals)
        temp_percent = [i/sum(temp_vals) * 100 for i in temp_vals]
        vals.append(temp_percent)
    final_vals.append(vals)


# For SnoBIRD on all species:
# If resulting pred overlap < total pred, adjust to total by saying the rest is 
# intergenic (if they could not overlap with annotated genome locations)
snoBIRD_species_val = []
for spe in all_spe:
    df = results_dict[spe]['snoBIRD']
    diff_len = models['snoBIRD'][spe] - len(df)
    df = pd.concat([df, pd.DataFrame([['']*6+['intergenic']+['']*6] * diff_len)])
    temp_vals = []
    for genomic in ['Exons', 'Introns', 'intergenic']:
        temp_vals.append(len(df[df['genomic_region'] == genomic]))
    temp_percent = [i/sum(temp_vals) * 100 for i in temp_vals]
    snoBIRD_species_val.append(temp_percent)

# Create stacked bar plots for human and pombe to compare the overlap of 
# predictions with exon/intron/intergenic per tool
rc = {'ytick.labelsize': 30, 'xtick.labelsize': 35}
plt.rcParams.update(**rc)
plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(2, 1, figsize=(20, 14))
plt.subplots_adjust(hspace=0.5)
genomic_elements = ['exonic', 'intronic', 'intergenic']
tools = ['snoBIRD', 'snoscan', 'snoreport2',  'infernal_rfam']
df_human = pd.DataFrame(final_vals[0], columns=genomic_elements, index=tools)
df_pombe = pd.DataFrame(final_vals[1], columns=genomic_elements, index=tools)
ax0 = df_human.plot.bar(rot=0, stacked=True, color=colors, ax=ax[0])
ax0.set_ylabel('Proportion of predictions (%)', fontsize=30)
ax0.set_title('Homo sapiens', fontdict={'fontsize': 25}, x=0.5, y=1.1)
ax1 = df_pombe.plot.bar(rot=0, stacked=True, color=colors, ax=ax[1])
ax1.set_xlabel('Tools', fontsize=30)
ax1.set_ylabel('Proportion of predictions (%)', fontsize=30)
ax1.set_title('Schizosaccharomyces pombe', fontdict={'fontsize': 25}, x=0.5, y=1.1)
fig.suptitle('Overlap of predictions with genomic elements per tool and species', 
            fontsize=32, x=0.5, y=1)
# Add total number of prediction per tool and per species
human_total_preds = [total_snoBIRD_preds['homo_sapiens'], 
                    total_snoscan_preds['homo_sapiens'], 
                    total_snoreport_preds['homo_sapiens'],  
                    total_infernal_preds['homo_sapiens']]
pombe_total_preds = [total_snoBIRD_preds['schizosaccharomyces_pombe'], 
                    total_snoscan_preds['schizosaccharomyces_pombe'], 
                    total_snoreport_preds['schizosaccharomyces_pombe'],  
                    total_infernal_preds['schizosaccharomyces_pombe']]                    
for rect, annotation in zip(ax0.patches, [f'({i})' for i in human_total_preds] * len(genomic_elements)):
    x = rect.get_x() + rect.get_width() / 2  # Center above the bar
    ax0.text(x, 106, annotation, ha='center', va='bottom', fontsize=15)
for rect, annotation in zip(ax1.patches, [f'({i})' for i in pombe_total_preds] * len(genomic_elements)):
    x = rect.get_x() + rect.get_width() / 2  # Center above the bar
    ax1.text(x, 106, annotation, ha='center', va='bottom', fontsize=15)

plt.savefig(snakemake.output.bar_human_pombe)


# Create stacked bar plot of SnoBIRD predictions overlap with 
# exon/intron/intergenic across all species
total_snoBIRD_preds_list = [f'{total_snoBIRD_preds[spe]}' for spe in all_spe]
pred_and_real = [snoBIRD_species_val, real_prop_list]
ft.grouped_stacked_bar2(pred_and_real, all_spe, genomic_elements, 
    'Overlap of SnoBIRD predictions with genomic elements (left), and actual genomic proportion (right) across species', 
    'Species',  'Proportion of predictions/genome (%)', colors, 'Genomic element', total_snoBIRD_preds_list, 
    snakemake.output.bar_other_species)
