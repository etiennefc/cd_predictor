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
import functions as ft

# Define model path and other variables
model_path = snakemake.input.model 
pretrained_model = snakemake.params.pretrained_model
num_labels = 2

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
kmer_pred = seq2kmer(pred_2545_2546, 6)
kmer_list = [i for i in kmer_pred.split(' ')]
kmer_w_special_tok = ['[CLS]'] + kmer_list + ['[SEP]']
tok = tokenizer(kmer_pred, return_tensors="pt", padding=True)

outputs = model(tok.input_ids, output_attentions=True)
probabilities = torch.softmax(outputs.logits, dim=1).to(device)
pred_labels = torch.argmax(probabilities, dim=1).to(device)
#print(tok)
##print(outputs.logits)
##print(outputs.attentions)
#print(outputs.attentions[0].shape)
#print(tok.input_ids)
#print(type(tok.input_ids))
#print(kmer_w_special_tok)



# build an explainer using a token masker
#explainer = shap.Explainer(f, tokenizer, labels=['expressed'])
#shap_values = explainer(kmer_pred)


pipe = TextClassificationPipeline(model=model, tokenizer=tokenizer)
#def score_and_visualize(text):
#    prediction = pipe([seq2kmer(text, 6)])
#    print(prediction[0])
#
#    explainer = shap.Explainer(pipe, seed=42, output_names=['Other', 'CD_snoRNA'], max_evals=100)
#    shap_values = explainer([text])
#    print(shap_values)
#    with open('shap_test2.html', 'w') as f2:
#        f2.write(shap.plots.text(shap_values, display=False, num_starting_labels=10))
#score_and_visualize(pred_2545_2546)

droso_pred = "GACTTCAGAAGAACTATAAAGGTAAATTACCATTATTTTTTTGTAGTTATTTATTTAGAATGATGATATTTTCCGTCCGTGGAACTTATACAAAGTAAAATCATGCTGACATTTTATGTCTTCTTTCCCCAACTGATTATTAACTTCAAATTTACGCTTGATGGCACTTAAATATCAGGATGCTAATTGTTATT"
ENSG00000221611 = "ACAAAAAAAAAAAAAAAAAAAAAAAAGAGGAAATGAATAAGTTAAAAAATTACATTCCATGATGTCCAGCACTGGGCTCTGATCACTCCGGAGGACACAGTTTTCCCCAAGACCATGGCTACCTGGGGATCTGAGAGAAAAAGGAAAATAATTATACAGTGGTTTGTATTATGAAAAAAGCAAAAATATAATGG"
Ce211_C_ele_Zemann_2006_B2 = "CCGCGGAACTCGCTAAGTGTCGTCTGCTTCCCTGCAAGACCTCATGAAATAGTGCATTTGGCACTGGCAGTGATGATCACAAATCCGTGTTTCTGACAAGCGATTGACGATAGAAAACCGGCTGAGCCAAAAAATTTTGGAAATTAAACCTGAATTGGTTTGCTCTTTGAGAAATAAACTTCAATCCGAATACT"


a = 'ATCGATGTGTGTGAAC'

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




def score_and_visualize(text):
    prediction = pipe([seq2kmer(text, 6)])
    print(prediction[0])

    explainer = shap.Explainer(pipe, seed=42, output_names=['Other', 'CD_snoRNA'], 
                            algorithm='partition', max_evals=500)
    shap_values = explainer([seq2kmer(text, 6)])
    shap_vals = shap_values.values[0]
    base_vals = shap_values.base_values
    kmer_vals = shap_values.data[0]
    # Select SHAP values for the positive label (1: C/D snoRNA)
    pos_shap_vals = shap_vals[:, 1]
    CLS = pos_shap_vals[0]
    SEP = pos_shap_vals[-1]
    other_tokens = pos_shap_vals[1:-2]
    avg_shap_per_nt = kmer_score_per_nt(text, other_tokens, 6)
    avg_shap_per_nt = [CLS] + avg_shap_per_nt + [SEP]

    # Plot the avg shap values across all tokens
    xticklabels = ['CLS'] + [i for i in text] + ['SEP']
    ft.shap_lineplot(range(len(avg_shap_per_nt)), avg_shap_per_nt, xticklabels, 'Sequence tokens', 
                    'Average SHAP values', '', 'testallo.svg')
    
    neg_shap_vals = shap_vals[:, 0]
    CLS = neg_shap_vals[0]
    SEP = neg_shap_vals[-1]
    other_tokens = neg_shap_vals[1:-2]
    avg_shap_per_nt = kmer_score_per_nt(text, other_tokens, 6)
    avg_shap_per_nt = [CLS] + avg_shap_per_nt + [SEP]

    # Plot the avg shap values across all tokens
    xticklabels = ['CLS'] + [i for i in text] + ['SEP']
    ft.shap_lineplot(range(len(avg_shap_per_nt)), avg_shap_per_nt, xticklabels, 'Sequence tokens', 
                    'Average SHAP values', '', 'testallo2.svg')

    # Create a force plot stacked on a text plot
    #with open('shap_test2.html', 'w') as f2:
    #    f2.write(shap.plots.text(shap_values, display=False, num_starting_labels=10))

score_and_visualize(pred_2545_2546)
#score_and_visualize(pred_2545_2546)
'''
def predict_proba(input_text):
    s = seq2kmer(input_text, 6)
    inputs = tokenizer(s, return_tensors='pt', padding=True)
    print(inputs.input_ids)
    outputs = model(inputs.input_ids)
    return outputs.logits.detach()

print(tokenizer)
tokenized_input = tokenizer(kmer_pred, return_tensors='pt', padding=True)
print(tokenized_input.input_ids)
print(tokenized_input.input_ids.shape)
explainer = shap.Explainer(predict_proba, tokenizer)
shap_values = explainer(tokenized_input.input_ids)


# Test with first layer
attention_matrix = outputs.attentions
all_layers_matrix = [attention_matrix[i].detach() for i in range(0, 12)]
#print(attention_matrix)

# An attention score (0-1, 0 being not important; 1 being highly important)
# is given for each pair of tokens so 191X191
# We're interested in the attention score of all tokens with regards to the 
# [CLS] token which is used at the end by the classification head
# [CLS] is always the first token
cls_attention = [matr[:, :, :, 0].squeeze(0) for matr in all_layers_matrix]
#print(cls_attention.shape)
cols = 2
rows = int(12/cols)

fig, axes = plt.subplots(cols, rows, figsize = (30, 14))
axes = axes.flat
print (f'Attention weights for token [CLS]')
for i,att in enumerate(cls_attention):

    #im = axes[i].imshow(att, cmap='gray')
    ax1 = sns.heatmap(att, vmin=0, vmax=1, ax = axes[i], xticklabels = kmer_w_special_tok)
    ax1.tick_params(axis='x', labelsize=1)
    axes[i].set_title(f'head #{i}' )
    axes[i].set_ylabel('layers')
plt.savefig('att.svg')


print('Second plot')
avg_attention = torch.mean(torch.stack(cls_attention, dim=0), dim=0)
f, a = plt.subplots(1, 1, figsize = (12, 12))
ax2 = sns.heatmap(avg_attention, xticklabels= kmer_w_special_tok)
ax2.tick_params(axis='x', labelsize=1)
a.set_ylabel('layers')
a.set_title(f'Avg head')
plt.savefig('avg_att.svg')

'''

