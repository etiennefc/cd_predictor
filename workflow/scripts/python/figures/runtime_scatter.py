#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import os 
import seaborn as sns

genome_dir = snakemake.input.genome_dir
colors = snakemake.params.colors
output = snakemake.output.scatter
# Load dfs
snoBIRD = [p for p in snakemake.input.snoBIRD if 
            'schizosaccharomyces_pombe' not in p] 
other_tools = [p for p in snakemake.input.other_tool if 
            'schizosaccharomyces_pombe' not in p] 

dfs = []
for p in snoBIRD + other_tools:
    df = pd.read_csv(p, sep='\t')
    dfs.append(df)

df_tools = pd.concat(dfs).reset_index(drop=True)
df_tools['tool'] = df_tools['tool'].replace('snoreport', 'snoreport2')

# Add genome size as column
species = list(pd.unique(df_tools['species_genome']))
sp_dict = {}
for sp in species:
    if 'chr1' in sp:
        size = os.path.getsize(f'{genome_dir}/homo_sapiens/1.fa')
    else:
        size = os.path.getsize(f'{genome_dir}/{sp}_genome.fa')
    size = size / (1024 ** 2)  # convert bytes in MB
    sp_dict[sp] = size

df_tools['genome_size'] = df_tools['species_genome'].map(sp_dict)

# Sort df based on ascending genome size and tool(for visualization in scatter)
df_tools = df_tools.sort_values(by='genome_size')
print(df_tools)
temp_dfs = []
for i, df in enumerate(df_tools.groupby('genome_size')):
    print(df[1])
    temp_df = []
    df_ = df[1]
    for t in ['snoscan', 'infernal_rfam', 'snoreport2', 'snoBIRD']:
        temp_df.append(df_[df_['tool'] == t])
    temp_dfs.append(pd.concat(temp_df))
df_tools = pd.concat(temp_dfs)
print(df_tools)

# Add pseudocount to runtime so that we can see all tools at 3600 min or above
#df_tools.loc[(df_tools['tool'] == 'infernal_rfam') & (df_tools['species_genome'] == 'gallus_gallus'), 'runtime_min'] = df_tools['runtime_min'] + 10
#df_tools.loc[(df_tools['tool'] == 'snoscan') & (df_tools['species_genome'] == 'gallus_gallus'), 'runtime_min'] = df_tools['runtime_min'] + 20
#df_tools.loc[(df_tools['tool'] == 'infernal_rfam') & (df_tools['species_genome'] == 'macaca_mulatta'), 'runtime_min'] = df_tools['runtime_min'] + 10
#df_tools.loc[(df_tools['tool'] == 'snoscan') & (df_tools['species_genome'] == 'macaca_mulatta'), 'runtime_min'] = df_tools['runtime_min'] + 20
df_tools['runtime_hour'] = df_tools['runtime_min'] / 60
df_tools.loc[(df_tools['tool'] == 'infernal_rfam') & (df_tools['species_genome'] == 'gallus_gallus'), 'runtime_hour'] = df_tools['runtime_hour'] + 0.1
df_tools.loc[(df_tools['tool'] == 'snoscan') & (df_tools['species_genome'] == 'gallus_gallus'), 'runtime_hour'] = df_tools['runtime_hour'] + 0.2
df_tools.loc[(df_tools['tool'] == 'infernal_rfam') & (df_tools['species_genome'] == 'macaca_mulatta'), 'runtime_hour'] = df_tools['runtime_hour'] + 0.1
df_tools.loc[(df_tools['tool'] == 'snoscan') & (df_tools['species_genome'] == 'macaca_mulatta'), 'runtime_hour'] = df_tools['runtime_hour'] + 0.2



plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(1, 1, figsize=(12, 10))
#sns.scatterplot(data=df_tools, x='genome_size', y='runtime_min', hue='tool', 
#                palette=colors, ax=ax, s=100, zorder=2, edgecolors=None, linewidth=0)
#ax.axhline(y=3600, color='grey', ls='--', linewidth=1, 
#            label='Max time limit (3600 min)', zorder=1)
sns.scatterplot(data=df_tools, x='genome_size', y='runtime_hour', hue='tool', 
                palette=colors, ax=ax, s=100, zorder=2, edgecolors=None, linewidth=0)
ax.axhline(y=60, color='grey', ls='--', linewidth=1, 
            label='Max time limit (60 h)', zorder=1)
for v in sp_dict.values():
    ax.axvline(x=v, color='lightgrey', ls='--', linewidth=1, alpha=0.5, zorder=1)
ax.set_xlabel('Genome size (Mb)', fontdict={'fontsize': 35})
#ax.set_ylabel('Run time (min)', fontdict={'fontsize': 35})
ax.set_ylabel('Run time (hour)', fontdict={'fontsize': 35})

plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(fontsize=20)
plt.savefig(output, dpi=600, bbox_inches='tight')
