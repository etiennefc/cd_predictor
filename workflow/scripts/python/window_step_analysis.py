#!/usr/bin/python3
import glob
import pandas as pd
import multiprocessing as mp
import subprocess as sp
from Bio import SeqIO
import numpy as np
from math import ceil
from pybedtools import BedTool
import time 
import torch
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification, logging
import utils as ut
#logging.set_verbosity_error()

genomes = glob.glob(f'{snakemake.input.genome}/*fa')

pretrained_model = str(snakemake.params.pretrained_model)
model_path = snakemake.input.snoBIRD
window_size = int(snakemake.wildcards.fixed_length)
step_size = snakemake.params.step_size
batch_size = 500
num_labels = 2
extension = 50

# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Load dfs
sno_df = pd.read_csv(snakemake.input.all_cd, sep='\t')
preds = pd.read_csv(snakemake.input.preds, sep='\t')

# Select only 1 example per snoRNA (no data augmentation)
preds = preds[(preds['y_true'] == 1) & (~preds['gene_id'].str.contains(
                                                            r'_[AB][0-5]$'))]
preds['gene_name'] = preds['gene_id']
preds = preds.merge(sno_df[['gene_id', 'species_name', 'chr', 'start', 'end', 
                'strand', 'sequence',f'extended_{window_size}nt_sequence', 
                'gene_biotype']])


# Load model and tokenizer
start_time = time.time()
tokenizer = AutoTokenizer.from_pretrained(pretrained_model)
model = BertForSequenceClassification.from_pretrained(pretrained_model, 
                            num_labels=num_labels)
model.load_state_dict(torch.load(model_path))
model.to(device)
model.classifier.to(device)
model.eval()
end_time = time.time()
print(f'elapsed time Loading model: {end_time -start_time}s')





# Change species_name for longer name
for k,v in snakemake.params.sp_name_dict.items():
    preds['species_name'] = preds['species_name'].replace(k, v)


# Define function to create sliding windows across a given interval
def sliding(start, end, left_shift, right_shift):
    slide_starts, slide_ends = [], []
    for i in reversed(range(1, left_shift+1)):
        slide_starts.append(start - i)
        slide_ends.append(end - i)
    slide_starts.append(start)
    slide_ends.append(end)
    for i in range(1, right_shift+1):
        slide_starts.append(start + i)
        slide_ends.append(end + i)
    return slide_starts, slide_ends



# Create bed for each species
l = len(preds)
bed_cols = ['chr', 'start', 'end', 'gene_id', 'score', 'strand', 'gene_name']
all_preds = []
for i, row in preds.iterrows():
    print('UPDATE:', (i+1)*100/l)
    sno_seq = row.sequence
    species_name = row.species_name
    genome = [g for g in genomes if species_name in g][0]
    ext_sno_seq = row[f'extended_{window_size}nt_sequence']
    ext_start = row.start - ext_sno_seq.find(sno_seq) 
    ext_end = ext_start + window_size - 1

    # Get sequence of the 50 sliding windows flanking the 
    # fixed_length nt original window
    all_starts, all_ends = sliding(ext_start, ext_end, 50, 50)
    temp_df = pd.DataFrame([[row.chr] * len(all_starts), all_starts, all_ends, 
                [row.gene_id]*len(all_starts), ['.']*len(all_starts), 
                [row.strand]*len(all_starts), 
                [row.gene_name]*len(all_starts)]).T
    temp_df.columns = bed_cols
    temp_df.to_csv('tempsliding.tsv', sep='\t', index=False, header=None)
    temp_bed = BedTool('tempsliding.tsv')
    fasta = temp_bed.sequence(fi=genome, nameOnly=True, s=True)
    temp_seqs = []
    with open(fasta.seqfn, 'r') as fasta_file:
        for line in fasta_file:
            if '>' not in line:
                temp_seqs.append(line.strip('\n'))

    # Kmerize and tokenize sequences
    kmer_seqs = [ut.seq2kmer(seq, 6) for seq in temp_seqs]
    eval_dataset = tokenizer(kmer_seqs, return_tensors='pt', 
                        padding=True).to(device)
    eval_dataset = TensorDataset(eval_dataset.input_ids, 
                        eval_dataset.attention_mask)
    eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)
    
    # Predict for each sequence if it contains sno or not
    ev_preds = [row.gene_id, row.chr, row.start, row.end, row.strand, 
                row.species_name, row.gene_biotype, row.y_true, row.y_pred]   
    for i, ev_batch in enumerate(eval_dataloader):
        #sp.call(f'echo EVAL BATCH {i+1}', shell=True)
        ev_input_ids, ev_attention_mask = ev_batch
        ev_input_ids = ev_input_ids.to(device)
        ev_attention_mask = ev_attention_mask.to(device)
        with torch.no_grad():  # nor gradient computation
            # Predict the labels of eval_dataset (returns logits here)
            # where logits are the model's prediction without applying any
            # activation function (pos: more probable; neg: less probable)
            output = model(ev_input_ids, attention_mask=ev_attention_mask)
            # Convert logits to probabilities via the softmax activation 
            # function (in 2nd dimension of output)
            probabilities = torch.softmax(output.logits, dim=1)
            # Get the predicted labels from these probabilities of each class
            # by getting the index (0 or 1) of the highest prob
            pred_labels = torch.argmax(probabilities, dim=1).tolist()  
            ev_preds.extend(pred_labels)
    print(ev_preds)
    all_preds.append(ev_preds)

sp.call('rm tempsliding.tsv', shell=True)

# Create final df
window_cols_left = ['window_-'+str(i) for i in reversed(range(1, 50+1))]
window_cols_right = ['window_'+str(i) for i in range(0, 50+1)]
other_cols = ['gene_id', 'chr', 'start', 'end', 'strand', 'species_name', 
            'gene_biotype', 'y_true', 'y_pred']
final_df = pd.DataFrame(all_preds, columns= (
                                                other_cols + 
                                                window_cols_left + 
                                                window_cols_right))

final_df.to_csv(snakemake.output.window_preds, sep='\t', index=False)
