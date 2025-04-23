#!/usr/bin/python3
import sys
import subprocess as sp
import pandas as pd
import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from math import ceil
import time 
import torch
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import (
                    AutoTokenizer, BertForSequenceClassification, 
                    logging, TextClassificationPipeline
)
import shap
from scipy.signal import find_peaks, savgol_filter
import utils as ut

# Define model path and other variables
model_path = sys.argv[1] 
FP_df = pd.read_csv(sys.argv[2], sep='\t')
FP_seq = pd.read_csv(sys.argv[3], sep='\t')
second_model_path = sys.argv[4]
pretrained_model = sys.argv[5]
fixed_length = int(sys.argv[6])
half_window = int((fixed_length + 2) / 2)
batch_size = int(sys.argv[7])
output_df = sys.argv[8]
num_labels = 2
len_c_box, len_d_box = 7, 4

# The smallest reliably annotated C/D snoRNA is 50 nt long (NR_145814, a
# C/D pseudogene in human. To not miss any snoRNA, we define the minimal length
# to find the C or D box to be +-15 nt from the center of the predicted window
min_box_dist = 15
# Flanking nt extending after the snoRNA start/end are minimally of 15 nt
# (+ 5 nt to get to the C_start or D_end), so
# no C or D box should be found in the first and last 20 nt of the window
flanking_nt = 20

# Thus, the ranges in which a C and D boxes can be searched are defined below
## + 1 to account for the [CLS] ##
C_range = range(flanking_nt + 1, half_window-min_box_dist)
D_range = range(half_window+min_box_dist, fixed_length - flanking_nt + 1)

# Select only FP examples from the first model and add extended sequence
FP_df = FP_df[(FP_df['y_true'] == 0) & (FP_df['y_pred'] == 1)]
FP_df = FP_df.merge(FP_seq[['gene_id', f'extended_{194}nt_sequence']],
        how='left', on='gene_id')

# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Limit the number of threads that torch can spawn with (to avoid core
# oversubscription) i.e. set the number of threads to the number of CPUs
# requested (not all CPUs physically installed)
#N_CPUS = os.environ.get("SLURM_CPUS_PER_TASK")
#torch.set_num_threads(int(N_CPUS))

# Force to not parallelize tokenizing before dataloader
# (causes forking errors otherwise)
#os.environ["TOKENIZERS_PARALLELISM"] = "false"

# Allow TF32 on matrix multiplication to speed up computations
torch.backends.cuda.matmul.allow_tf32 = True

# Allow TF32 when using cuDNN library (GPU-related library usually
# automatically installed on the cluster)
torch.backends.cudnn.allow_tf32 = True


# Show packages versions
sp.call(f'echo PANDAS VERSION: {pd.__version__}', shell=True)
sp.call(f'echo TORCH VERSION: {torch.__version__}', shell=True)
sp.call(f'echo NUMPY VERSION: {np.__version__}', shell=True)
sp.call(f'echo TRANSFORMERS VERSION: {transformers.__version__}', shell=True)
sp.call(f'echo IS CUDA AVAILABLE?: {torch.cuda.is_available()}', shell=True)

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
sp.call(f'echo LOAD INITIAL MODEL: {end_time-start_time}', shell=True)


# Load all sequences, define batch size and 
# pipeline used for SHAP computation
all_seqs = list(FP_df[f'extended_{fixed_length}nt_sequence'])
all_gene_ids = list(FP_df['gene_id'])
pipe2 = TextClassificationPipeline(model=model, tokenizer=tokenizer, 
        device=device, batch_size=batch_size)


# Make prediction again with first model and compute SHAP values
all_scores = []
for batch in ut.batch_generator(all_seqs, batch_size):
    batch_scores = ut.shap_batch(batch, all_gene_ids, pipe2)
    all_scores.extend(batch_scores)

shap_df = pd.DataFrame(all_scores, columns =
                    ['gene_id', 'predicted_label', 'probability_CD', 'CLS'] + 
                    [f'SHAP_pos_{i}' for i in range(int(fixed_length))] +  
                    ['SEP'])


shap_df = shap_df.merge(FP_df[[
            'gene_id', f'extended_{fixed_length}nt_sequence']])
shap_cols = [i for i in shap_df.columns if i.startswith('SHAP_')]


# Based on SHAP values, find the C and D box
# Then find C' and D' boxes and the overall box score
FP_sno_df = ut.find_all_boxes(shap_df, fixed_length, shap_cols, C_range,
                            D_range, flanking_nt)

# Get snoRNA sequence, length and normalized_mfe
FP_sno_df['sno_seq'] =  FP_sno_df.apply(ut.get_seq, axis=1)
FP_sno_df['sno_length'] = (FP_sno_df.D_END + 5) - (FP_sno_df.C_START - 5)

FP_sno_df = ut.sno_mfe(FP_sno_df, f'extended_{fixed_length}nt_sequence',
                    'sno_length')

# Get terminal stem stability/score and combined
FP_sno_df = ut.terminal_stem(FP_sno_df, f'extended_{fixed_length}nt_sequence')


# Load second model 
start_time = time.time()
second_model = BertForSequenceClassification.from_pretrained(pretrained_model,
            num_labels=num_labels)
second_model.load_state_dict(torch.load(second_model_path))
second_model.to(device)
second_model.classifier.to(device)
second_model.eval()
end_time = time.time()
sp.call(f'echo LOAD SECOND MODEL: {end_time-start_time}', shell=True)

# Load sequences for the second model prediction 
test_seqs = list(FP_sno_df[f'extended_{fixed_length}nt_sequence'])
kmer_seqs = [ut.seq2kmer(s, 6) for s in test_seqs]

# Tokenize test data in right format and create dataloader
eval_dataset = tokenizer(kmer_seqs, return_tensors='pt', 
                padding=True).to(device)
# Fake labels because they are all FP of the first model
eval_labels = torch.tensor([0] * len(test_seqs)).to(device)  
eval_dataset = TensorDataset(eval_dataset.input_ids, 
                        eval_dataset.attention_mask, eval_labels)
eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)

# Predict with second snoBIRD model
ev_preds, probs = [], []
for i, ev_batch in enumerate(eval_dataloader):
    sp.call(f'echo 2nd MODEL EVAL BATCH {i+1}', shell=True)
    ev_input_ids, ev_attention_mask, ev_batch_labels = ev_batch
    ev_input_ids = ev_input_ids.to(device)
    ev_attention_mask = ev_attention_mask.to(device)
    ev_batch_labels = ev_batch_labels.to(device)
    with torch.no_grad():  # no gradient computation
        # Predict the labels of eval_dataset (returns logits here)
        # where logits are the model's prediction without applying any
        # activation function (pos: more probable; neg: less probable)
        output = second_model(ev_input_ids, attention_mask=ev_attention_mask,
                                labels=ev_batch_labels)
        
        # Convert logits to probabilities via the softmax activation 
        #  function (in 2nd dimension of output) 
        probabilities = torch.softmax(output.logits, dim=1)
        
        # Get the probability of prediction of the winning class
        highest_prob, _ = torch.max(probabilities, dim=1)
        probs += highest_prob.tolist()

        # Get the predicted labels from these probabilities of each class
        pred_labels = torch.argmax(
            probabilities, dim=1)  # get the index (0 or 1) of the highest prob
        ev_preds += pred_labels.tolist()
        


# Filter predictions using the feature values
df = pd.DataFrame({'second_model_prediction': ev_preds, 
                'second_model_probability': probs})
df = pd.concat([FP_sno_df[['gene_id']].reset_index(drop=True), df], axis=1)
FP_sno_df = FP_sno_df.merge(df, on='gene_id', how='left')

FP_sno_df = ut.feature_filters(FP_sno_df, 'second_model_probability', 
            'second_model_prediction', 'second_model_filtered_prediction')


# Save final df
final_cols = ['gene_id', 'probability_CD', 'C_MOTIF', 'C_START', 'C_END',
            'D_MOTIF', 'D_START', 'D_END', 'C_PRIME_MOTIF', 'C_PRIME_START', 
            'C_PRIME_END', 'D_PRIME_MOTIF', 'D_PRIME_START', 'D_PRIME_END', 
            'box_score', 'terminal_combined', 'normalized_sno_stability', 
            'second_model_probability', 'second_model_prediction', 
            'second_model_filtered_prediction']

FP_sno_df[final_cols].to_csv(output_df, sep='\t', index=False)


