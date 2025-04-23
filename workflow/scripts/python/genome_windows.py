#!/usr/bin/python3
import pandas as pd
import multiprocessing as mp
import subprocess as sp
from Bio import SeqIO
import numpy as np
from math import ceil
import time 
import torch
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification, logging
#logging.set_verbosity_error()

genome = snakemake.input.genome
pretrained_model = str(snakemake.params.pretrained_model)
model_path = snakemake.input.snoBIRD
window_size = snakemake.params.fixed_length
step_size = snakemake.params.step_size
batch_size = 16
num_labels = 2
#genome = 'data/references/genome_fa/homo_sapiens_genome.fa'

# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# Load model and tokenizer 
tokenizer = AutoTokenizer.from_pretrained(pretrained_model)
model = BertForSequenceClassification.from_pretrained(pretrained_model, num_labels=num_labels)
model.load_state_dict(torch.load(model_path)) 
model.to(device)
model.classifier.to(device)
model.eval()


def seq2kmer(seq, k):
    """
    Convert original str sequence to kmers (str of all kmers)
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers


start_time = time.time()
end_time = time.time()

# seqio human/candida: 8.6097s/0.0358s
print(f'elapsed time: {end_time -start_time}s')

########## Compare biopython, awk and julia split of fasta efficiency###########

###>2h (then stopped) for human with mem_bathc_size=10


chr_dict = {record.id: str(record.seq) 
            for record in SeqIO.parse(genome, "fasta")}

"""Increase len(kmer_seqs)??"""

# with yield mem_batch_size=10, batch_size=16 et window_size=190 --> 5,4 min; peak ram ~7GB (human)
def scan_fasta(chr_seq, window_size, mini_batch_size, step_size=1, mem_batch_size=100):
    id_chr = [chr_id for chr_id, s in chr_dict.items() if s == chr_seq][0]
    len_chr = len(chr_seq)

    if len_chr > window_size:
        for i in range(0, int(len_chr - window_size / step_size + 1), int(window_size+mini_batch_size*mem_batch_size)):
            print(f'{id_chr}: {i/int(len_chr - window_size / step_size + 1)*100}')
            big_batch, starts, ends = [], [], []
            for j in range(i, min(i+window_size+mini_batch_size*mem_batch_size, len_chr), step_size):
                window = chr_seq[j*step_size:j*step_size+window_size]
                if len(window) == window_size: # to account for the last batch which would contain smaller sequences
                    starts.append(j*step_size + 1)
                    ends.append(j*step_size + 1 + window_size)
                    big_batch.append(window)
            kmer_seqs = [seq2kmer(seq, 6) for seq in big_batch]
            eval_dataset = tokenizer(kmer_seqs, return_tensors='pt', padding=True).to(device)
            eval_dataset = TensorDataset(eval_dataset.input_ids, eval_dataset.attention_mask)
            eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)

            # where len(kmer_seqs) = window_size + mini_batch_size * mem_batch_size i.e. 350
            yield eval_dataloader, starts, ends   

        # Keep start/end and strand of sequences in batch
    elif len_chr == window_size: 
        n_windows = 1
        ###TO adapt to 190 nt sequences as input
    else:
        raise ValueError(f'Length of sequence {id_chr} ({len_chr}) must be >= than length of window size ({window_size})')
 





## 112s with multiprocessing (overhead of splitting data); 126 without
num_processes = mp.cpu_count()


# Instiate model and tokenizer



# How much time to create all 1 nt step windows in the human genome????
# In candida with eval_dataloader ~2h for the largest chr

def predict(chr_seq, window_size, mini_batch_size):
    for x in scan_fasta(chr_seq, window_size, mini_batch_size):
        # batch_dataloader = x[0]
        # kmer_starts = x[1]
        # kmer_ends = x[2]
        l = 2
        #print(len(x[0]), len(x[1]), len(x[2]), [chr_id for chr_id, s in chr_dict.items() if s == chr_seq][0])



# Multi-processing to split the chromosomes
start_time = time.time()
if len(chr_dict.keys()) > 1:  # or for strand - ???
    with mp.Pool() as pool:
        #result = pool.starmap(scan_fasta, [(seq_, window_size, batch_size, step_size) for seq_ in chr_dict.values()])
        allo = pool.starmap(predict, [(seq_, window_size, batch_size) for seq_ in chr_dict.values()])
    end_time = time.time()
    
    print(f'elapsed time: {end_time -start_time}s')

