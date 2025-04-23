#!/usr/bin/python3
import pandas as pd
import subprocess as sp
from pybedtools import BedTool
import numpy as np
import time
import torch
from torch.utils.data import TensorDataset, DataLoader
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification, logging
#logging.set_verbosity_error()
import onnx
import onnxruntime

# Load variables
pretrained_model = str(snakemake.params.pretrained_model)
model_path = snakemake.input.snoBIRD_pseudo
window_size = int(snakemake.params.fixed_length)
step_size = snakemake.params.step_size
batch_size = 16
num_labels = 2
output = snakemake.output.df

# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Create function to kmerize sequences
def seq2kmer(seq, k):
    """
    Convert original str sequence to kmers (str of all kmers)
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers


# Load model and tokenizer 
start_time = time.time()
tokenizer = AutoTokenizer.from_pretrained(pretrained_model)
model = BertForSequenceClassification.from_pretrained(pretrained_model, num_labels=num_labels)
model.load_state_dict(torch.load(model_path)) 
model.to(device)
model.classifier.to(device)
model.eval()
end_time = time.time()
print(f'LOAD INITIAL MODEL: {end_time-start_time}')

# Load sequences to predict in list and kmerize it into dataloader
df_window = pd.read_csv(snakemake.input.sno_pseudo_windows, sep='\t')
windows = list(df_window.seq)

kmer_seqs = [seq2kmer(seq, 6) for seq in windows]
eval_dataset = tokenizer(kmer_seqs, return_tensors='pt', padding=True).to(device)
eval_dataset = TensorDataset(eval_dataset.input_ids, eval_dataset.attention_mask)
eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)



# Load model in onnx for first batch only to speed up prediction after
start_time = time.time()
for i, ev_batch in enumerate(eval_dataloader):
    ev_input_ids, ev_attention_mask = ev_batch
    ev_input_ids = ev_input_ids.to(device)
    ev_attention_mask = ev_attention_mask.to(device)
    with torch.no_grad():  # nor gradient computation
        torch.onnx.export(model, (ev_input_ids, ev_attention_mask), 'model_temp_sno_pseudo.onnx', 
                        input_names=["input_ids", "attention_mask"], output_names=["output"])
        so = onnxruntime.SessionOptions()
        ####To setup running predictions on one core at a time per chr (but each predictions take longer)
        #so.execution_mode = onnxruntime.ExecutionMode.ORT_SEQUENTIAL
        #so.intra_op_num_threads = 1
        sp.call('export OMP_NUM_THREADS=1', shell=True)
        ort_session = onnxruntime.InferenceSession('model_temp_sno_pseudo.onnx', so)
        if i == 0:
            break
end_time = time.time()
print(f'Load model in onnx : {end_time -start_time}s')


# Predict with onnx on these windows
total_batch = len(eval_dataloader)
print(total_batch)
preds_probs = []
for i, ev_batch in enumerate(eval_dataloader):
    s_time = time.time()
    ev_input_ids, ev_attention_mask = ev_batch
    ev_input_ids = ev_input_ids.to(device)
    ev_attention_mask = ev_attention_mask.to(device)
    with torch.no_grad():  # nor gradient computation
        ev_input_ids_np = ev_input_ids.cpu().numpy()
        ev_attention_mask_np = ev_attention_mask.cpu().numpy()

        if len(ev_input_ids_np) != batch_size:  # to account for the last batch at the end of chromosome
            len_diff = batch_size - len(ev_input_ids_np)
            print('\n\nBATCH MISMATCH\n\n')
            # Pad it on the right side to be same size as mini_batch_size
            original_len = len(ev_input_ids_np)
            ev_input_ids_np = np.pad(ev_input_ids_np, ((0, len_diff), (0, 0)), 'constant', constant_values=0)
            ev_attention_mask_np = np.pad(ev_attention_mask_np, ((0, len_diff), (0, 0)), 'constant', constant_values=0)
            outputs = ort_session.run(None, {"input_ids": ev_input_ids_np, "attention_mask": ev_attention_mask_np})
            # Get only the prediction for the relevant examples, not the padding
            outputs = outputs[0][0:original_len]
            probabilities = torch.softmax(torch.from_numpy(outputs), dim=1).to(device)
            pred_labels = torch.argmax(probabilities, dim=1).to(device)
            
            # Save probabilities and pred_labels as one tensor converted to dataframe
            # get probabilities of 1 (expressed C/D)
            probs_expressed = probabilities[:,1].to(device)
            concat_tensor = torch.cat((pred_labels.unsqueeze(1), probs_expressed.unsqueeze(1)), dim=1).numpy()
            pred_prob_df = pd.DataFrame(concat_tensor, columns=['predicted_label', 'probability'])
            print(pred_prob_df)
            preds_probs.append(pred_prob_df)

        else:
            outputs = ort_session.run(None, {"input_ids": ev_input_ids_np, "attention_mask": ev_attention_mask_np})
            #print(outputs)
            probabilities = torch.softmax(torch.from_numpy(outputs[0]), dim=1).to(device)
            pred_labels = torch.argmax(probabilities, dim=1).to(device)

            # Save probabilities and pred_labels as one tensor converted to dataframe
            # get probabilities of 1 (expressed C/D)
            probs_expressed = probabilities[:,1].to(device)
            concat_tensor = torch.cat((pred_labels.unsqueeze(1), probs_expressed.unsqueeze(1)), dim=1).numpy()
            pred_prob_df = pd.DataFrame(concat_tensor, columns=['predicted_label', 'probability'])
            print(pred_prob_df)
            preds_probs.append(pred_prob_df)

    n_time = time.time()
    print(f'onnx pred time: {n_time -s_time}s BATCH {i}/{total_batch} ({i/total_batch*100}%)')

prediction_df = pd.concat(preds_probs)
print(prediction_df)
prediction_df = pd.concat([df_window['gene_id'].reset_index(drop=True), 
                prediction_df.reset_index(drop=True)], axis=1)
print(prediction_df)

prediction_df.to_csv(snakemake.output.df, sep='\t', index=False)
# Compute precision and recall based on the expressed C/D in cerevisiae
# We want in theory a perfect recall, and some FP are acceptable, but 
# it should be less than snoreport by a good margin!

# Find if expressed C/D overlap with sno predicted as sno by second predictor

sp.call('rm model_temp_sno_pseudo.onnx', shell=True)
