#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches

# Load prediction df and merge with relevant info
colors = snakemake.params.target_colors
snoBIRD_test = pd.read_csv(snakemake.input.snoBIRD_test_set, sep='\t')
test_set = pd.read_csv(snakemake.input.test_set, sep='\t')
sp_dict = dict(zip(test_set.gene_id, test_set.species_name))
target_dict = dict(zip(test_set.gene_id, test_set.target))

# Select only examples that were accurately predicted as snoRNA pseudogenes 
# or expressed CD snoRNA
#snoBIRD_test['species'] = snoBIRD_test['chr'].map(sp_dict)
snoBIRD_test['predicted_label'] = snoBIRD_test['predicted_label'].replace(
                                'CD_snoRNA_pseudogene', 'snoRNA_pseudogene')
snoBIRD_test['target'] = snoBIRD_test['chr'].map(target_dict)
snoBIRD_test = snoBIRD_test[snoBIRD_test['target'] == snoBIRD_test['predicted_label']]

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

snoBIRD_test['C_score'] = snoBIRD_test['C_MOTIF'].apply(lambda x: hamming(x))
snoBIRD_test['C_prime_score'] = snoBIRD_test['C_PRIME_MOTIF'].apply(
                                        lambda x: hamming(x))
snoBIRD_test['D_score'] = snoBIRD_test['D_MOTIF'].apply(
                                        lambda x: hamming(x, consensus='CTGA'))
snoBIRD_test['D_prime_score'] = snoBIRD_test['D_PRIME_MOTIF'].apply(
                                        lambda x: hamming(x, consensus='CTGA'))


# Create density plots of box scores
cols = ['box_score', 'C_score', 'C_prime_score', 'D_score', 'D_prime_score']
df_list = []
l = ['expressed_CD_snoRNA', 'snoRNA_pseudogene']
for col in cols:
    temp = []
    for exp in l:
        temp.append(snoBIRD_test[snoBIRD_test['predicted_label'] == exp])
    df_list.append(temp)

print()


plt.rcParams['svg.fonttype'] = 'none'
#fig, ax = plt.subplots(1, 5, figsize=(25, 14))
output_scores = [snakemake.output.box_score, snakemake.output.c_score, snakemake.output.c_prime_score, snakemake.output.d_score, snakemake.output.d_prime_score]
for i, col in enumerate(df_list):
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    for j, df in enumerate(col):
        sns.kdeplot(data=df[cols[i]], fill=True, ax=ax, color=colors[l[j]])
    ax.set_xlabel(cols[i], fontdict={'fontsize': 40})
    ax.set_ylabel('Density', fontdict={'fontsize': 40})
    ax.tick_params(axis='x', labelsize=35)
    ax.tick_params(axis='y', labelsize=35)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    legend_list = []
    for k, v in colors.items():
        legend_element = mpatches.Patch(color=v, label=k)
        legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0,1.1),
            fontsize=20)

    plt.savefig(output_scores[i], bbox_inches='tight', dpi=500)


# Create density plots of normalized sno stability and terminal stem score
output_stability = [snakemake.output.terminal_stability, snakemake.output.structure_stability] 
cols = ['terminal_stem_score', 'normalized_sno_stability']
l = ['expressed_CD_snoRNA', 'snoRNA_pseudogene']
df_list = []
for col in cols:
    temp = []
    for exp in l:
        temp.append(snoBIRD_test[snoBIRD_test['predicted_label'] == exp])
    df_list.append(temp)


plt.rcParams['svg.fonttype'] = 'none'
for i, col in enumerate(df_list):
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    for j, df in enumerate(col):
        sns.kdeplot(data=df[cols[i]], fill=True, ax=ax, color=colors[l[j]])
    ax.set_xlabel(cols[i], fontdict={'fontsize': 40})
    ax.set_ylabel('Density', fontdict={'fontsize': 40})
    ax.tick_params(axis='x', labelsize=35)
    ax.tick_params(axis='y', labelsize=35)
    ax.spines['right'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    legend_list = []
    for k, v in colors.items():
            legend_element = mpatches.Patch(color=v, label=k)
            legend_list.append(legend_element)
    plt.legend(handles=legend_list, loc='upper right', bbox_to_anchor=(0,1.1), fontsize=20)
    plt.savefig(output_stability[i], bbox_inches='tight', dpi=500)




