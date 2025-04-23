#!/usr/bin/python3
import sys
import pandas as pd
import torch
import torch.nn as nn
import random
import numpy as np
import subprocess as sp
import sklearn
from sklearn.metrics import f1_score, accuracy_score
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification, logging
#logging.set_verbosity_error()

# Load inputs and params
pretrained_model = str(snakemake.params.pretrained_model)  # pretrained DNABert6 model
#best_hyperparams = pd.read_csv(sys.argv[4], sep='\t')
batch_size = int(16)  # nb of example per batch
num_labels = 2
model_path = str(snakemake.input.model)
fixed_length = str(snakemake.input.X_test).split('nt.ts')[0].split('_')[-1]

# Load test data and create tensor matrix (format used by pytorch)
X_test = pd.read_csv(snakemake.input.X_test, sep='\t')

X_tun = pd.read_csv(f'data/references/positives_and_negatives/initial/initial_tuning_set_fixed_length_{fixed_length}nt.tsv', sep='\t')
X_train = pd.read_csv(f'data/references/positives_and_negatives/initial/initial_training_set_fixed_length_{fixed_length}nt.tsv', sep='\t')

X_test = pd.concat([X_tun, X_test, X_train])
print(X_test)
X_test = X_test[X_test['target'] != 'other'].reset_index(drop=True)

# Select only mouse
X_test = X_test[(X_test['species_name'] == 'mus_musculus') & (X_test['target'] == 'snoRNA_pseudogene')]





y_test = pd.read_csv(snakemake.input.y_test, sep='\t')
y_tun = pd.read_csv(f'data/references/positives_and_negatives/initial/initial_tuning_target_fixed_length_{fixed_length}nt.tsv', sep='\t')
y_train = pd.read_csv(f'data/references/positives_and_negatives/initial/initial_training_target_fixed_length_{fixed_length}nt.tsv', sep='\t')
y_test = pd.concat([y_test, y_tun, y_train])

y_simple = y_test[y_test['gene_id'].isin(X_test.gene_id)]
y_simple = y_simple.drop(columns=['gene_id']).reset_index(drop=True)
y_simple['target'] = y_simple['target'].replace(1, 0)
y_simple['target'] = y_simple['target'].replace(2, 1)

print(X_test)
print(y_simple)

# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Load outputs
df_metrics = snakemake.output.df_metrics_on_test
df_preds = snakemake.output.test_predictions


# Transform sequence of examples in test set into kmers (6-mers)
def seq2kmer(seq, k):
    """
    Convert original str sequence to kmers (str of all kmers)
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers

test_seqs = list(X_test[f'extended_{fixed_length}nt_sequence'])
#test_seqs = list(X_test[f'test'])
kmer_seqs = [seq2kmer(s, 6) for s in test_seqs]

# Tokenize test data in right format and create dataloader
tokenizer = AutoTokenizer.from_pretrained(pretrained_model)  # BertTokenizerFast
eval_dataset = tokenizer(kmer_seqs, return_tensors='pt', padding=True).to(device)
eval_labels = torch.tensor(list(y_simple.target)).to(device)
eval_dataset = TensorDataset(eval_dataset.input_ids, eval_dataset.attention_mask, eval_labels)
eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)

# Load model
model = BertForSequenceClassification.from_pretrained(pretrained_model, num_labels=num_labels)
model.load_state_dict(torch.load(model_path)) 
model.to(device)
model.classifier.to(device)

# Test the model
model.eval()  # no more dropout
ev_preds, ev_labels = [], []
for i, ev_batch in enumerate(eval_dataloader):
    sp.call(f'echo EVAL BATCH {i+1}', shell=True)
    ev_input_ids, ev_attention_mask, ev_batch_labels = ev_batch
    ev_input_ids = ev_input_ids.to(device)
    ev_attention_mask = ev_attention_mask.to(device)
    ev_batch_labels = ev_batch_labels.to(device)
    with torch.no_grad():  # nor gradient computation
        # Predict the labels of eval_dataset (returns logits here)
        # where logits are the model's prediction without applying any
        # activation function (positive: more probable; negative: less probable)
        output = model(ev_input_ids, attention_mask=ev_attention_mask, labels=ev_batch_labels)
        
        # Convert logits to probabilities via the softmax activation function (in 2nd dimension of output 
        probabilities = torch.softmax(output.logits, dim=1)
        
        # Get the predicted labels from these probabilities of each class
        pred_labels = torch.argmax(probabilities, dim=1)  # get the index (0 or 1) of the highest prob
        ev_preds += pred_labels.tolist()
        ev_labels += ev_batch_labels.tolist()


# Save predictions
df = pd.DataFrame({'y_true': ev_labels, 'y_pred': ev_preds})
df = pd.concat([X_test[['gene_id']].reset_index(drop=True), df], axis=1)
df.to_csv(df_preds, sep='\t', index=False)

# Compute test prediction metrics
fscore = f1_score(ev_labels, ev_preds, average='macro')  # macro avg across 2 classes
sp.call(f'echo FSCORE {fscore}', shell=True)
accuracy = accuracy_score(df.y_true, df.y_pred)
sp.call(f'echo Accuracy {accuracy}', shell=True)
TP_sno = len(df[(df.y_true == 1) & (df.y_pred == 1)])
FP_sno = len(df[(df.y_true != 1) & (df.y_pred == 1)])
FN_sno = len(df[(df.y_true == 1) & (df.y_pred != 1)])
if TP_sno + FP_sno == 0:
    precision_sno = 0
else:
    precision_sno = TP_sno/(TP_sno + FP_sno)
if TP_sno + FN_sno == 0:
    recall_sno = 0
else:
    recall_sno = TP_sno/(TP_sno + FN_sno)
sp.call(f'echo PRECISION sno: {precision_sno}', shell=True)
sp.call(f'echo RECALL sno: {recall_sno}', shell=True)


TP_pseudosno = len(df[(df.y_true == 0) & (df.y_pred == 0)])
FP_pseudosno = len(df[(df.y_true != 0) & (df.y_pred == 0)])
FN_pseudosno = len(df[(df.y_true == 0) & (df.y_pred != 0)])
if TP_pseudosno + FP_pseudosno == 0:
    precision_pseudosno = 0
else:
    precision_pseudosno = TP_pseudosno/(TP_pseudosno + FP_pseudosno)
if TP_pseudosno + FN_pseudosno == 0:
    recall_pseudosno = 0
else:
    recall_pseudosno = TP_pseudosno/(TP_pseudosno + FN_pseudosno)
sp.call(f'echo PRECISION pseudosno: {precision_pseudosno}', shell=True)
sp.call(f'echo RECALL pseudosno: {recall_pseudosno}', shell=True)





metrics_df = pd.DataFrame([[accuracy, fscore, precision_sno,
                            recall_sno]],
                            columns=['accuracy_2_classes', 'f1_score_2_classes',
                            'precision_sno', 'recall_sno'])
metrics_df.to_csv(df_metrics, sep='\t', index=False)

