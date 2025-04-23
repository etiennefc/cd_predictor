#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
import numpy as np

df = pd.read_csv(snakemake.input.df, sep='\t')
output = snakemake.output.lineplot
output_strand = snakemake.output.lineplot_strand
output_density = snakemake.output.density

# Select only sno for which the central window is accurately predicted
df = df[df['window_0'] == 1]

# Create density to show the length distribution of + vs - strand snoRNA genes
df['SnoRNA length (nt)'] = df.end - df.start + 1
plt.rcParams['svg.fonttype'] = 'none'
sns.kdeplot(df, x='SnoRNA length (nt)', hue='strand', fill=True, 
            palette={'-': 'red', '+': 'blue'})
plt.savefig(output_density, dpi=600, bbox_inches='tight')

# Drop some cols
df = df.drop(columns=['start', 'end', 'SnoRNA length (nt)'])

# Subset according to different sno characteristics
sno = df[df['gene_biotype'] == 'expressed_CD_snoRNA']
pseudo = df[df['gene_biotype'] == 'snoRNA_pseudogene']
plus = df[df['strand'] == '+']
minus = df[df['strand'] == '-']

# Get proportion of predicted C/D (1) at each position 
# (from -50 to +50 around window)
prop = [i*100 for i in list(df.filter(regex='^window_').mean(
                                                        numeric_only=True))]
prop2 = [i*100 for i in list(sno.filter(regex='^window_').mean(
                                                        numeric_only=True))]
prop3 = [i*100 for i in list(pseudo.filter(regex='^window_').mean(
                                                        numeric_only=True))]
prop4 = [i*100 for i in list(plus.filter(regex='^window_').mean(
                                                        numeric_only=True))]
prop5 = [i*100 for i in list(minus.filter(regex='^window_').mean(
                                                        numeric_only=True))]

# Create lineplot for all CD, expressed vs pseudogenes
rc = {'ytick.labelsize': 25, 'xtick.labelsize': 25}
plt.rcParams.update(**rc)
plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(1, 1, figsize=(12, 10))
for i in range(-50, 50, 5):
    ax.vlines(x=i, ymin=0, ymax=100, 
            colors='lightgrey', linestyle='dashed')
for j in reversed(range(5, 100, 5)):
    ax.hlines(y=j, xmin=-50, xmax=50, 
            colors='lightgrey', linestyle='dashed')
ax.plot([-i for i in reversed(range(1, 50+1))] + [i for i in range(0, 50+1)], 
                            prop, color='black', label='All C/D')
ax.plot([-i for i in reversed(range(1, 50+1))] + [i for i in range(0, 50+1)], 
                            prop2, color='#abca77ff', label='Expressed C/D')
ax.plot([-i for i in reversed(range(1, 50+1))] + [i for i in range(0, 50+1)], 
                            prop3, color='#824c92ff', label='C/D pseudogene')

ax.set_xlabel('Distance to original window (nt)', fontdict={'fontsize': 30})
ax.set_ylabel('Proportion of windows\npredicted as C/D snoRNA (%)', 
                fontdict={'fontsize': 30})
plt.legend()
plt.margins(x=0, y=0)
fig.suptitle('Window step analysis', fontsize=35, weight='bold', x=0.5, y=1)
plt.savefig(output, dpi=600, bbox_inches='tight')

# Create lineplot for all CD and per strand in which the snoRNA is encoded
rc = {'ytick.labelsize': 25, 'xtick.labelsize': 25}
plt.rcParams.update(**rc)
plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(1, 1, figsize=(12, 10))
for i in range(-50, 50, 5):
    ax.vlines(x=i, ymin=0, ymax=100, 
            colors='lightgrey', linestyle='dashed')
for j in reversed(range(5, 100, 5)):
    ax.hlines(y=j, xmin=-50, xmax=50, 
            colors='lightgrey', linestyle='dashed')
ax.plot([-i for i in reversed(range(1, 50+1))] + [i for i in range(0, 50+1)], 
                            prop, color='black', label='All C/D')
ax.plot([-i for i in reversed(range(1, 50+1))] + [i for i in range(0, 50+1)], 
                            prop4, color='blue', label='+ strand C/D')
rev = []
for i in reversed(prop5):
    rev.append(i)
ax.plot([-i for i in reversed(range(1, 50+1))] + [i for i in range(0, 50+1)], 
                            rev, color='red', label='- strand C/D')

ax.set_xlabel('Distance to original window (nt)', fontdict={'fontsize': 30})
ax.set_ylabel('Proportion of windows\npredicted as C/D snoRNA (%)', 
                fontdict={'fontsize': 30})
plt.legend()
plt.margins(x=0, y=0)
fig.suptitle('Window step analysis per strand of snoRNA genes', fontsize=30, weight='bold', x=0.5, y=1)
plt.savefig(output_strand, dpi=600, bbox_inches='tight')



