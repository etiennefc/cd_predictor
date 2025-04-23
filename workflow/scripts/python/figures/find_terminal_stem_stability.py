#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import subprocess as sp
import re

df = pd.read_csv(snakemake.input.box_df, sep='\t')
fixed_length_col = [i for i in df.columns if 'predicted_extended' in i][0]
fasta_output = snakemake.output.fasta_terminal_stem
density_biotype = snakemake.output.density_biotype
df_output = snakemake.output.terminal_stem_df


# Create flanking seq columns
df['left_20nt_flanking'] = df['predicted_extended_sequence'].apply(
                                lambda x: x[0:20])
df['right_20nt_flanking'] = df['predicted_extended_sequence'].apply(
                                lambda x: x[-20:])

flanking_dict = dict(zip(df.gene_id, zip(df.left_20nt_flanking, 
                    df.right_20nt_flanking)))


# Create fasta required for RNAcofold 
# (reversed left seq + '&' + reversed right seq in fasta format)
with open('temp_terminal_stem.fa', 'w') as f:
    for id_, flanking_nts in flanking_dict.items():
        f.write(f'>{id_}\n')
        # Reverse both flanking seq so that it is 
        # correctly represented in the RNAcofold plot
        reverse_left = flanking_nts[0][::-1]
        reverse_right = flanking_nts[1][::-1]
        f.write(f'{reverse_left}&{reverse_right}\n')

# Run RNAcofold
sp.call(f'RNAcofold < temp_terminal_stem.fa > {fasta_output}', shell=True)
sp.call('rm temp_terminal_stem.fa *.ps', shell=True)


# Get terminal_stem stability and length score
# The terminal stem length score is defined as: 
# score = paired_nt - intramol_paired_nt - gap_number
terminal_mfe_dict = {}
terminal_stem_length_dict = {}
with open(fasta_output, 'r') as f:
    for line in f:
        if line.startswith('>'):
            gene_id = line.strip('>\n')
        elif line.startswith('.') | line.startswith('('):
            terminal_mfe = float(line.split(' ', maxsplit=1)[1].strip(' \n()'))
            terminal_mfe_dict[gene_id] = terminal_mfe

            # From the dot bracket, extract the number of paired nt
            dot_bracket = line.split(' ', maxsplit=1)[0]
            dot_bracket = dot_bracket[0:20].replace('\n', '')  # we select only the left flanking region (20 nt)
            paired_base = dot_bracket.count('(')
            intramolecular_paired = dot_bracket.count(')')  # these ')' are intramolecular paired nt
                                                            # (i.e. nt pair within left sequence only)
            # This finds all overlapping (and non-overlapping) gaps of 1 to 19 nt inside the left flanking region
            gaps = re.findall(r'(?=(\(\.{1,19}\())', dot_bracket)
            number_gaps = ''.join(gaps)  # join all gaps together in one string
            number_gaps = len(re.findall('\.', number_gaps))  # count the number of nt gaps in sequence
            
            stem_length_score = paired_base - intramolecular_paired - number_gaps
            if stem_length_score < 0:
                stem_length_score = 0
            terminal_stem_length_dict[gene_id] = stem_length_score

df['terminal_stem_stability'] = df['gene_id'].map(terminal_mfe_dict)
df['terminal_stem_length_score'] = df['gene_id'].map(terminal_stem_length_dict)
cols = ['gene_id', 'gene_biotype', 'species_name', 'left_20nt_flanking', 
        'right_20nt_flanking', 'terminal_stem_stability',
        'terminal_stem_length_score', fixed_length_col]
df[cols].to_csv(df_output, sep='\t', index=False)

sns.jointplot(data=df, x='terminal_stem_stability', 
                y='terminal_stem_length_score', hue='gene_biotype', 
                fill=True, kind='kde')
plt.savefig(density_biotype)



