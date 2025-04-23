#!/usr/bin/python3
import glob
import subprocess as sp
import pandas as pd

input_rDNA = snakemake.input.target_rDNA
genome = snakemake.input.genome
output = snakemake.output.predictions

# Split rDNA of all species in separate fasta files
sp.call("csplit -sz -f temp_rDNA_ "+input_rDNA+" '/>/' '{*}'", shell=True)

species = snakemake.wildcards.species
print(species)
species_rDNA = []
for fa in glob.glob('temp_rDNA_*'):
    with open(fa, 'r') as f:
        for i, line in enumerate(f):
            if species in line:
                species_rDNA.append(fa)
# For a given species, concat all of its rDNA fa into 1 fa (ex: 5.8S, 18S and 28S)
sp.call('cat '+' '.join(species_rDNA)+' > temp_rDNA_'+species+'.fa', shell=True)


# Run snoscan on each species 
sp.call("snoscan temp_rDNA_"+species+".fa "+genome+" -o "+output, shell=True)


sp.call('rm temp_rDNA_*', shell=True)