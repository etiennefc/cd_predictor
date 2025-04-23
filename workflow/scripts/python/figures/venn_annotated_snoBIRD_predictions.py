#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
from matplotlib_venn import venn2
import matplotlib.patches as mpatches
from glob import glob

#inputs = snakemake.input.snoBIRD_preds
inputs = glob('results/predictions/snoBIRD/bedgraph_overlap/multiple_filters/*_TGIRT_coverage.tsv')

# Load dfs
dfs = []
for p in inputs:
    df = pd.read_csv(p, sep='\t')
    dfs.append(df)
sum_df = pd.concat(dfs).sort_values('species').reset_index(drop=True)

# Add manual inspection of the potential novel expressed, 
# i.e. real novel unannotated snoRNAs based on bedgraphs
d = {'drosophila_melanogaster': 0, 'gallus_gallus': 18, 'macaca_mulatta': 8, 
    'schizosaccharomyces_pombe': 8, 'tetrahymena_thermophila': 14, 'homo_sapiens': 22, 
    'plasmodium_falciparum': 0, 'danio_rerio': 32, 'arabidopsis_thaliana': 0}

# Filter total annotated C/D in macaca_mulatta < 164 nt after merging overalpping entries 
sum_df.loc[sum_df['species'] == 'macaca_mulatta', 'total_annot_cd'] = 422

sum_df['validated_new_cd'] = sum_df['species'].map(d)
print(sum_df)
c_list = sum_df[['total_annot_exp', 'validated_new_cd', 'exp_annot_preds']].values.tolist()


# get the difference to get which are annotated expressed but not predicted by SnoBIRD
c_list = [[c[0]-c[2], c[1], c[2]] for c in c_list]  
l = ['Annotated and expressed C/D', 'Expressed SnoBIRD prediction (new C/D)', 
    'Expressed annotated C/D and predicted by SnoBIRD']

# Adjust values after manual inspection
c_list[1][2] = 77
c_list[1][0] = 3
c_list[3][2] = 67
c_list[5][2] = 247

print(c_list)
species_l = sum_df['species'].values.tolist()


def venn_multiple(y, x, count_list, labels, colors, ax_title, title, path, **kwargs):
    """
    Creates x by y venn diagrams for which we have values for the left, right and middle (intersect) values
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, axes = plt.subplots(y, x, figsize=(25, 14))
    plt.subplots_adjust(hspace=0.5)
    ax = axes.flatten()
    for i, element in enumerate(count_list):
        counts = count_list[i]
        ax[i].set_title(ax_title[i], fontdict={'fontsize': 25}, x=0.5, y=1)
        v = venn2(subsets = counts, set_labels=(None, None), set_colors=colors, ax=ax[i], **kwargs)
        for text in v.subset_labels:
            if text:
                if text.get_text() == '0':
                    text.set_text('')
                text.set_fontsize(20)
        colors_ = [patch.get_facecolor() for patch in v.patches if patch is not None]
    fig.suptitle(title, x=0.5, y=1, fontsize=25)
    legend_list = []
    for i, crit in enumerate(labels):
        legend_element = mpatches.Patch(color=colors_[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, bbox_to_anchor=(0.4, 0.1), fontsize=20)
    plt.savefig(path, dpi=600)

# Create venn
venn_multiple(3, 3, c_list, l, ('magenta', 'lightgreen'), species_l, 
    'Overlap between expressed annotated C/D and expressed SnoBIRD predictions', 
    snakemake.output.venn)


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d} \n ({p:.1f}%)'.format(p=pct,v=val)
    return my_autopct

def donut_2(y, x, count_list, ax_title, labels, colors, title, legend_labels, legend_colors, path,
            **kwargs):
    """
    Creates a donut plot with two layers. Counts, labels and colors are nested
    lists for outer and inner layers attributes. The counts are represented as
    numerical and percentage values in parentheses.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    fig, axes = plt.subplots(y, x, figsize=(25, 14))
    plt.subplots_adjust(hspace=0.5)
    ax = axes.flatten()
    for i, counts in enumerate(count_list):
        ax[i].axis('equal')
        ax[i].set_title(ax_title[i], fontdict={'fontsize': 25}, x=0.5, y=1.06)
        outer_donut, _, _ = ax[i].pie(counts[0], radius=1.3, labels=labels[0],
                            colors=colors[0], autopct=make_autopct(counts[0]), pctdistance=0.85, **kwargs)
        plt.setp(outer_donut, width=0.4, edgecolor='white')

        inner_donut, _, _ = ax[i].pie(counts[1], radius=1.3-0.4, labels=labels[1],
                        colors=colors[1], autopct=make_autopct(counts[1]), pctdistance=0.8, **kwargs)
        plt.setp(inner_donut, width=0.4, edgecolor='white')
    legend_list = []
    for i, crit in enumerate(legend_labels):
        legend_element = mpatches.Patch(color=legend_colors[i], label=crit)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper left', bbox_to_anchor=(-0.1, 1),
                fontsize=25)
    fig.suptitle(title, y=1, fontsize=25)

    plt.savefig(path, dpi=600)


# Create donut of proportion of annotated expressed and not expressed that are 
# predicted (or not) by snoBIRD
sp_donut = []
for sp in species_l:
    sp_df = sum_df[sum_df['species'] == sp]
    exp_cd = sp_df['total_annot_exp'].values
    not_exp_cd = sp_df['total_annot_cd'].values - exp_cd
    exp_pred = sp_df['exp_annot_preds'].values
    exp_not_pred = exp_cd - exp_pred
    not_exp_pred = sp_df['total_annot_preds'].values - exp_pred
    not_exp_not_pred = not_exp_cd - not_exp_pred
    sp_donut.append([[exp_cd[0], not_exp_cd[0]], [exp_pred[0], exp_not_pred[0], not_exp_pred[0], not_exp_not_pred[0]]])

donut_2(3, 3, sp_donut, species_l, [['', ''], ['', '', '', '']], 
    [['#810f7c', '#8c6bb1'], ['beige', 'lightgrey', 'beige', 'lightgrey']], 
    'Proportion of annotated C/D that are expressed or not and predicted or not by SnoBIRD', 
    ['expressed_annotated_CD', 'not_expressed_annotated_CD', 'snoBIRD_predicted', 'not_predicted'], 
    ['#810f7c', '#8c6bb1', 'beige', 'lightgrey'], snakemake.output.donut)

