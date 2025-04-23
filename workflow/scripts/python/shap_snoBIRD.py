#!/usr/bin/python3
import sys
import pandas as pd
import subprocess as sp
import os
import numpy as np
from math import ceil
import time 
import torch
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification, logging, TextClassificationPipeline
import matplotlib.pyplot as plt 
import seaborn as sns 
import shap
#import functions as ft

# Define model path and other variables
model_path = snakemake.input.model 
pretrained_model = snakemake.params.pretrained_model
num_labels = 2
fixed_length = snakemake.params.fixed_length

# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")




# Show packages versions
sp.call(f'echo PANDAS VERSION: {pd.__version__}', shell=True)
sp.call(f'echo TORCH VERSION: {torch.__version__}', shell=True)
sp.call(f'echo NUMPY VERSION: {np.__version__}', shell=True)
sp.call(f'echo TRANSFORMERS VERSION: {transformers.__version__}', shell=True)
sp.call(f'echo IS CUDA AVAILABLE?: {torch.cuda.is_available()}', shell=True)

# Load model and tokenizer 
start_time = time.time()
tokenizer = AutoTokenizer.from_pretrained(pretrained_model)
model = BertForSequenceClassification.from_pretrained(pretrained_model, num_labels=num_labels)
model.load_state_dict(torch.load(model_path))
model.to(device)
model.classifier.to(device)
model.eval()
end_time = time.time()
sp.call(f'echo LOAD INITIAL MODEL: {end_time-start_time}', shell=True)


def seq2kmer(seq, k):
    """
    Convert original str sequence to kmers (str of all kmers)
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers


# Predict on given sequence
pred_2545_2546 = "ATGTTCGCCGTATGCTGATTGATGTCTGAGTACAACGAAGCTGTTCTGATGATATTATTATGTCTTCTTTACCCGTCCTATCGACTATGACGATAACTTTGTCTCCACTATAGGTGGAGACCATACTTGCCAGTATCTTCGATGTCTGAATTGTTTTCGTTGTTCCAGTATCTTAGTAATCTCTGTATTGTGTA"


pipe = TextClassificationPipeline(model=model, tokenizer=tokenizer, device=device)



def kmer_score_per_nt(seq, scores, k):
    # Get the avg SHAP value per nt across 
    # the kmers in which that nt is found
    seq_scores = []
    L = len(seq)
    firsts = [i for i in range(k-1)]
    lasts = [i for i in range(L-k, L)]
    for i, nt in enumerate(seq):
        if i in firsts:
            sc = scores[0:i+1]
        elif i in lasts:
            sc = scores[-(L-i):]
        else:
            sc = scores[i+1-k:i+1]
        avg_score = np.mean(sc)
        seq_scores.append(avg_score)
    return seq_scores





    


def score_and_visualize(text, gene_id):
    st = time.time()
    prediction = pipe([seq2kmer(text, 6)])
    predicted_label = prediction[0]['label']
    print(predicted_label)
    pred_label_dict = {'LABEL_0': 'Other', 'LABEL_1': 'CD_snoRNA'}
    prob = prediction[0]['score']
    print(prediction[0])

    explainer = shap.Explainer(pipe, seed=42, output_names=['Other', 'CD_snoRNA'], 
                            algorithm='partition', max_evals=500)
    shap_values = explainer([seq2kmer(text, 6)])
    #shap_values = shap.Explainer(pipe([seq2kmer(text, 6)]), seed=42, output_names=['Other', 'CD_snoRNA'], 
    #                        algorithm='partition', max_evals=5)
    shap_vals = shap_values.values[0]
    #base_vals = shap_values.base_values
    #kmer_vals = shap_values.data[0]
    # Select SHAP values for the positive label (1: C/D snoRNA)
    pos_shap_vals = shap_vals[:, 1]
    CLS = pos_shap_vals[0]
    SEP = pos_shap_vals[-1]
    other_tokens = pos_shap_vals[1:-2]
    avg_shap_per_nt = kmer_score_per_nt(text, other_tokens, 6)
    avg_shap_per_nt = [gene_id, pred_label_dict[predicted_label], prob, CLS] + avg_shap_per_nt + [SEP]
    nd = time.time()
    sp.call(f'echo SHAP time {gene_id}: {nd-st}s', shell=True)
    
    #return avg_shap_per_nt

    # Plot the avg shap values across all tokens
    #xticklabels = ['CLS'] + [i for i in text] + ['SEP']
    #ft.shap_lineplot(range(len(avg_shap_per_nt)), avg_shap_per_nt, xticklabels, 'Sequence tokens', 
    #                'Average SHAP values', '', 'testpombe.svg')
    


    # Create a force plot stacked on a text plot
    #with open('shap_test2.html', 'w') as f2:
    #    f2.write(shap.plots.text(shap_values, display=False, num_starting_labels=10))




#cd = pd.read_csv('data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_194nt.tsv', sep='\t')
#all_seqs = list(cd[f'extended_{fixed_length}nt_sequence'])
#all_gene_ids = list(cd['gene_id'])
batch_size= 8
pipe2 = TextClassificationPipeline(model=model, tokenizer=tokenizer, device=device, batch_size=batch_size)




def batch_generator(sequences, batch_size):
    for i in range(0, len(sequences), batch_size):
        yield sequences[i:i+batch_size]

def shap_batch(text_batch, gene_ids):
    st = time.time()
    pred_label_dict = {'LABEL_0': 'Other', 'LABEL_1': 'CD_snoRNA'}
    
    # Perform prediction on the batch
    prediction_batch = pipe2([seq2kmer(text, 6) for text in text_batch])
    predicted_labels = [pred['label'] for pred in prediction_batch]
    probabilities = [pred['score'] for pred in prediction_batch]
    
    # Compute SHAP values for the batch
    explainer = shap.Explainer(pipe2, seed=42, output_names=['Other', 'CD_snoRNA'], 
                            algorithm='partition', max_evals=500)
    shap_values_batch = explainer([seq2kmer(text, 6) for text in text_batch])
    shap_values_batch = shap_values_batch.values
    
    # Process SHAP values for each sequence in the batch
    avg_shap_per_nt_batch = []
    for i, (shap_values, gene_id, predicted_label, prob) in enumerate(zip(shap_values_batch, gene_ids, predicted_labels, probabilities)):
        pos_shap_vals = shap_values[:, 1]
        CLS = pos_shap_vals[0]
        SEP = pos_shap_vals[-1]
        other_tokens = pos_shap_vals[1:-2]
        avg_shap_per_nt = kmer_score_per_nt(text_batch[i], other_tokens, 6)
        avg_shap_per_nt = [gene_id, pred_label_dict[predicted_label], prob, CLS] + avg_shap_per_nt + [SEP]
        avg_shap_per_nt_batch.append(avg_shap_per_nt)
    
    nd = time.time()
    sp.call(f'echo SHAP time per batch: {nd-st}s', shell=True)
    
    return avg_shap_per_nt_batch


#pred_2545_2546 = list(cd[f'extended_{fixed_length}nt_sequence'])[0]

all_scores = []
all_seqs = [pred_2545_2546, pred_2545_2546]
all_gene_ids = ['id1', 'id2']
for batch in batch_generator(all_seqs, batch_size):
    batch_scores = shap_batch(batch, all_gene_ids)
    all_scores.extend(batch_scores)

df = pd.DataFrame(all_scores, columns=['gene_id', 'predicted_label', 'probability', 'CLS'] + \
                    [f'SHAP_pos_{i}' for i in range(fixed_length)] + ['SEP'])
df.to_csv(snakemake.output.shap_df, sep='\t', index=False)

xticklabels = ['CLS'] + [i for i in pred_2545_2546] + ['SEP']
avg_shap_per_nt = list(df.iloc[0, 3:])
print(avg_shap_per_nt)

def shap_lineplot(x_val, y_val, xtick_labels, xlabel, ylabel, title, path, **kwargs):
    """ 
    Create a vertical connected dot plot or lineplot.
    """
    plt.rcParams['svg.fonttype'] = 'none'
    rc = {'ytick.labelsize': 25, 'xtick.labelsize': 25}
    plt.rcParams.update(**rc)    
    plt.subplots(1, 1, figsize=(30, 8))
    plt.plot(x_val, y_val)
    plt.xlabel(xlabel, fontdict={'fontsize': 30})
    plt.ylabel(ylabel, fontdict={'fontsize': 30})
    plt.xticks(range(len(xtick_labels)), xtick_labels, fontsize=12)
    plt.margins(x=0)
    #ax.set_ylim(0, 1.05)
    plt.title(title, fontsize=30, x=0.5, y=1)
    plt.savefig(path, dpi=600, bbox_inches='tight')

shap_lineplot(range(len(avg_shap_per_nt)), avg_shap_per_nt, xticklabels, 'Sequence tokens', 
                'Average SHAP values', '', 'testpombe.svg')


#score_and_visualize(pred_2545_2546, 'testid')
#d = {'id1': pred_2545_2546, 'id2': pred_2545_2546}
#
#avg_shaps = []
#for gene_id, seq in d.items():
#    avg_shaps.append(score_and_visualize(seq, gene_id))
#
#df = pd.DataFrame(avg_shaps, columns=['gene_id', 'predicted_label', 'probability', 'CLS'] + \
#                    [f'SHAP_pos_{i}' for i in range(fixed_length)] + ['SEP'])
#df.to_csv(snakemake.output.shap_df, sep='\t', index=False)