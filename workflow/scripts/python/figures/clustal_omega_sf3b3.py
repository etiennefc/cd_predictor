#!/usr/bin/python3
import matplotlib.pyplot as plt 
import subprocess as sp
from Bio import Phylo



# Set input/output variables
input_nwk = snakemake.input.phylo_tree_nwk
tree_fig = snakemake.output.tree_fig

# Create phylo tree figure
tree_ = Phylo.read(input_nwk, 'newick')
tree_.root_at_midpoint()
plt.rcParams['svg.fonttype'] = 'none'
fig, ax = plt.subplots(1, 1, figsize=(8, 10))
Phylo.draw(tree_, axes=ax)
plt.xticks(fontsize=14)

# Increase font size of x-axis label
ax.set_xlabel("Evolutionary Distance", fontsize=16)

plt.savefig(tree_fig)
