#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import subprocess as sp

df = pd.read_csv(snakemake.input.box_df, sep='\t')
fixed_length_col = [i for i in df.columns if 'nt_sequence' in i][0]
fasta_output = snakemake.output.fasta_mfe
density_biotype = snakemake.output.density_biotype
df_output = snakemake.output.mfe_df


seq_dict = dict(zip(df.gene_id, df.predicted_sequence))
df['len'] = df['predicted_sequence'].apply(lambda x: len(x))

# Create fasta required for RNAfold 
with open('temp_mfe.fa', 'w') as f:
    for id_, predicted_seq in seq_dict.items():
        f.write(f'>{id_}\n{predicted_seq}\n')

# Run RNAfold
sp.call("sed 's/>/>SNOBIRD_/g' temp_mfe.fa > SNOBIRD_rna_fold.fa", shell=True)
sp.call(f'RNAfold --infile=SNOBIRD_rna_fold.fa --outfile=temp_snoBIRD_structure.mfe', shell=True)
sp.call("sed -i 's/SNOBIRD_//g' temp_snoBIRD_structure.mfe", shell=True)
sp.call(f"mv temp_snoBIRD_structure.mfe {fasta_output}", shell=True)
sp.call('rm temp_mfe.fa SNOBIRD_rna_fold.fa *.ps', shell=True)


# Get structure stability (MFE)
mfe_dict = {}
with open(fasta_output, 'r') as f:
    for line in f:
        if line.startswith('>'):
            gene_id = line.strip('>\n')
        elif line.startswith('.') | line.startswith('('):
            mfe = float(line.split(' ', maxsplit=1)[1].strip(' \n()'))
            mfe_dict[gene_id] = mfe


# Create cols for structure stability and normalized MFE by the predicted length
df['sno_stability'] = df['gene_id'].map(mfe_dict)
df['normalized_sno_stability'] = df['sno_stability'] / df['len']
cols = ['gene_id', 'gene_biotype', 'species_name', 'sno_stability', 'normalized_sno_stability', 'predicted_sequence', fixed_length_col]
df[cols].to_csv(df_output, sep='\t', index=False)

sns.kdeplot(data=df, x='normalized_sno_stability', hue='gene_biotype', fill=True)
plt.savefig(density_biotype)



