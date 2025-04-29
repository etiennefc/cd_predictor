#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from glob import glob

# Load prediction df and merge with relevant info
all_spe = ['schizosaccharomyces_pombe', 'plasmodium_falciparum', 
        'tetrahymena_thermophila', 'drosophila_melanogaster', 'gallus_gallus', 
        'danio_rerio', 'macaca_mulatta',  'homo_sapiens']
simplified_spe = [s.split('_')[0][0].upper() + '_' + s.split('_')[1] for s in all_spe]
colors = snakemake.params.target_colors
dfs = []
for p in snakemake.input.snoBIRD:
    path = glob(p)[0]
    species = path.split('/')[4]
    snoBIRD = pd.read_csv(path, sep='\t')
    snoBIRD['species'] = species
    dfs.append(snoBIRD)

snoBIRD_preds = pd.concat(dfs)
snoBIRD_preds['predicted_label'] = snoBIRD_preds['predicted_label'].replace(
                                'CD_snoRNA_pseudogene', 'snoRNA_pseudogene')

# Define score for specific boxes
def hamming(motif, consensus='RTGATGA'):
    """ 
    Compute Hamming distance with regards to consensus motif.
    """
    score = 0
    for i, s in enumerate(motif):
        if i == 0:
            if consensus == 'RTGATGA':
                if s not in ['A', 'G']:
                    score +=1
            else:  # for D box
                if s != consensus[i]:
                    score +=1
        else:
            if s != consensus[i]:
                score +=1
    return score

snoBIRD_preds['C_score'] = snoBIRD_preds['C_MOTIF'].apply(lambda x: hamming(x))
snoBIRD_preds['C_prime_score'] = snoBIRD_preds['C_PRIME_MOTIF'].apply(
                                        lambda x: hamming(x))
snoBIRD_preds['D_score'] = snoBIRD_preds['D_MOTIF'].apply(
                                        lambda x: hamming(x, consensus='CTGA'))
snoBIRD_preds['D_prime_score'] = snoBIRD_preds['D_PRIME_MOTIF'].apply(
                                        lambda x: hamming(x, consensus='CTGA'))




# Create density plots of box scores per species
cols = ['box_score', 'C_score', 'C_prime_score', 'D_score', 'D_prime_score', 
    'terminal_stem_score', 'normalized_sno_stability']
df_list = []
l = ['expressed_CD_snoRNA', 'snoRNA_pseudogene']
for col in cols:
    temp2 = []
    for sp in all_spe:
        temp_df = snoBIRD_preds[snoBIRD_preds['species'] == sp]
        temp = []
        for exp in l:
            temp.append(temp_df[temp_df['predicted_label'] == exp])
        temp2.append(temp)
    df_list.append(temp2)



plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(7, 8, figsize=(29, 33))
for i, col in enumerate(df_list):
    for j, spe in enumerate(col):
        ax[i, j].set_xlabel(cols[i], fontdict={'fontsize': 25})
        ax[i, j].tick_params(axis='x', labelsize=25)
        ax[i, j].tick_params(axis='y', labelsize=25)
        for k, expression_df in enumerate(spe):
            sns.kdeplot(data=expression_df[cols[i]], fill=True, ax=ax[i, j], color=colors[l[k]])
        if i == 0:
            ax[i, j].set_title(simplified_spe[j], fontdict={'fontsize': 28, 'fontweight':'bold'}, x=0.5, y=1.1)
        if j == 0:
            ax[i, j].set_ylabel('Density', fontdict={'fontsize': 25})
        else:
            ax[i, j].set_ylabel(None)
        
plt.tight_layout()
legend_list = []
for k, v in colors.items():
    legend_element = mpatches.Patch(color=v, label=k)
    legend_list.append(legend_element)
plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0,1.1),
            fontsize=20)

plt.savefig(snakemake.output.boxes_and_stability, bbox_inches='tight', dpi=500)


