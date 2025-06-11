#!/usr/bin/python3
import pandas as pd
from Bio import SeqIO 
import subprocess as sb 

# Load inputs / outputs
sp = str(snakemake.wildcards.species)
chr_fasta_dir = snakemake.output.chr_fasta_dir
output_gtf = snakemake.output.gtf
target_id_output = snakemake.output.target_ids
rDNA = snakemake.input.rDNA_species
rDNA_danio = snakemake.input.rDNA_danio
sno_df = pd.read_csv(snakemake.input.sno_table, sep='\t')
sno_df = sno_df[sno_df['species'] == sp]

# Create output fasta of snoRNA sequences
sno_dict = dict(zip(sno_df.prediction_id, sno_df.predicted_sequence))
with open(snakemake.output.sno_fasta, 'w') as f:
    for k,v in sno_dict.items():
        f.write(f'>{k}\n{v}\n')


def gtf_attribute(df, col_names):
    """
    Create an attribute column as last column of gtf file.
    """
    return df.apply(lambda row: '; '.join(
                    [f'{col} "{row[col]}"' for col in col_names]), axis=1)


def make_gtf_from_df(df):
    """
    Create a Gene transfer format (.gtf) file from a dataframe.
    """
    final_cols = ['gene_id', 'chr', 'source', 'feature', 'start', 'end', 
                'score', 'strand', 'frame', 'attributes']
    gene_cols = ['gene_id', 'gene_version', 'gene_name', 'gene_source', 
                'gene_biotype']
    transcript_cols = gene_cols + ['transcript_id', 'transcript_version', 
                'transcript_name', 'transcript_source', 'transcript_biotype', 
                'tag', 'transcript_support_level']
    exon_cols = transcript_cols + ['exon_id', 'exon_number', 'exon_version']
    given_cols = [gene_cols, transcript_cols, exon_cols]
    
    gtf = df.copy()
    gtf['gene_name'] = gtf['gene_id']
    gtf[['gene_version', 'transcript_version', 'exon_version']] = '1'
    gtf['transcript_id'] = gtf['gene_id'] + '.t1'
    gtf['transcript_name'] = gtf['gene_id'] + '-1'
    gtf[['score', 'frame']] = "."
    gtf[['source', 'gene_source', 'transcript_source']] = 'GenBank_NCBI'
    gtf['exon_id'] = gtf['gene_id'] + '.e1'
    gtf['exon_number'] = '1'
    gtf['tag'] = 'basic'
    gtf['transcript_support_level'] = 'NA'

    gtf[['gene_biotype', 'transcript_biotype']] = 'rRNA'

    # Build 3 gtf per snoRNA prediction (one gene line, one transcript line 
    # and one exon line)
    gtfs = []
    for i, feature in enumerate(['gene', 'transcript', 'exon']):
        temp_gtf = gtf.copy()
        temp_gtf['feature'] = feature
        if i == 0:  # add box info only for gene line
            temp_gtf['attributes'] = gtf_attribute(temp_gtf, given_cols[i])
        else:
            temp_gtf['attributes'] = gtf_attribute(temp_gtf, given_cols[i])
        temp_gtf['attributes'] = temp_gtf['attributes'] + ';'
        temp_gtf = temp_gtf[final_cols]
        gtfs.append(temp_gtf) 
    
    # Sort gtf by gene_id and drop gene_id column
    final_gtf = pd.concat(gtfs)
    final_gtf = final_gtf.sort_values('gene_id').drop(columns='gene_id')

    return final_gtf

# Create fasta dir for each rDNA as a chr and create a gtf of the rDNA 
# genes only
if sp == 'danio_rerio':
    file_rDNA = rDNA_danio 
else:
    file_rDNA = rDNA

sb.call(f'mkdir -p {chr_fasta_dir}', shell=True)


rDNA_dict = SeqIO.to_dict(SeqIO.parse(file_rDNA, "fasta"))
gtfs = []
for k, v in rDNA_dict.items():
    seq_ = v.seq
    desc = v.description
    if sp in desc:
        with open(f'{chr_fasta_dir}/{k}.fa', 'w') as f2:
            f2.write(f'>{k}\n{seq_}\n')
        
        # Create gtf of the rRNA entries per species
        initial_df = pd.DataFrame([[k, 1, len(seq_), '+', k]], 
                    columns=['chr', 'start', 'end', 'strand', 'gene_id'])
        gtf_df = make_gtf_from_df(initial_df)
        gtfs.append(gtf_df)

        # Create target id file so that snoglobe can predict only on the rRNA 
        # targets of interest
        sb.call(f'echo {k} >> {target_id_output}', shell=True)

final_gtf = pd.concat(gtfs)
final_gtf.to_csv(output_gtf, sep='\t', index=False, header=None)
sb.call(f'''sed -i 's/""/"/g; s/"gene_id/gene_id/g; s/;"/;/g' {output_gtf}''', 
    shell=True)
        




