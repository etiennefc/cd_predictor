#!/usr/bin/python3
import sys
import pandas as pd
import torch
import torch.nn as nn
import os
import random
import numpy as np
import subprocess as sp
import sklearn
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, accuracy_score, precision_recall_fscore_support
from torch.utils.data import TensorDataset, DataLoader
import torch.optim.lr_scheduler as lr_scheduler
import transformers
from transformers import AutoTokenizer, BertForSequenceClassification, logging
#logging.set_verbosity_error()

# Load params
pretrained_model = snakemake.params.pretrained_model  # pretrained DNABert6 model
fold_num = str(10)  # training fold number
rs = 42  # random_state

# Load inputs
X_train = pd.read_csv(snakemake.input.X_train, sep='\t')
print(X_train)
#X_train = X_train[~X_train['gene_id'].isin(['URS0000531800_6239', 'URS00004F99EB_3702', 'URS000037BF83_6239', 'URS00006C129A_9606'])].drop_duplicates(subset=['extended_194nt_sequence'])
X_train = X_train.drop_duplicates(subset=['extended_194nt_sequence'])
y_train = pd.read_csv(snakemake.input.y_train, sep='\t')
y_train = y_train.iloc[X_train.index].reset_index(drop=True)
X_train = X_train.reset_index(drop=True)
y_simple = y_train.drop(columns=['gene_id']).reset_index(drop=True)
print(X_train)
print(y_simple)
best_hyperparams = pd.read_csv(snakemake.input.best_hyperparams, sep='\t')
fixed_length = 194

# Convert sno labels so that expressed and pseudogene 
# snoRNAs are considered the same label (i.e. 1)
y_simple = y_simple.replace(2, 1)
# Get path of outputs
output_model = snakemake.output.model
output_loss = snakemake.output.fold_loss
output_f1 = snakemake.output.fold_f1_score

# Show packages versions
sp.call(f'echo PANDAS VERSION: {pd.__version__}', shell=True)
sp.call(f'echo TORCH VERSION: {torch.__version__}', shell=True)
sp.call(f'echo NUMPY VERSION: {np.__version__}', shell=True)
sp.call(f'echo SKLEARN VERSION: {sklearn.__version__}', shell=True)
sp.call(f'echo TRANSFORMERS VERSION: {transformers.__version__}', shell=True)
sp.call(f'echo IS CUDA AVAILABLE?: {torch.cuda.is_available()}', shell=True)

# Set reproducible randomness
torch.backends.cudnn.deterministic = True
random.seed(rs)
torch.manual_seed(rs)
torch.cuda.manual_seed(rs)
np.random.seed(rs)

# Define if we use GPU (CUDA) or CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# Limit the number of threads that torch can spawn with (to avoid core oversubscription)
# i.e. set the number of threads to the number of CPUs requested (not all CPUs physically installed)
#N_CPUS = os.environ.get("SLURM_CPUS_PER_TASK")
N_CPUS = 20
torch.set_num_threads(int(N_CPUS))

# Force to not parallelize tokenizing before dataloader (causes forking errors otherwise)
os.environ["TOKENIZERS_PARALLELISM"] = "false"

# Allow TF32 on matrix multiplication to speed up computations
#torch.backends.cuda.matmul.allow_tf32 = True

# Allow TF32 when using cuDNN library (GPU-related library usually automatically installed on the cluster)
#torch.backends.cudnn.allow_tf32 = True


# Transform sequence of examples in training set into kmers (6-mers)
def seq2kmer(seq, k):
    """
    Convert original str sequence to kmers (str of all kmers)
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers

train_seqs = list(X_train[f'extended_{fixed_length}nt_sequence'])
kmer_seqs = [seq2kmer(s, 6) for s in train_seqs]


# Load pre-trained DNABERT model
tokenizer = AutoTokenizer.from_pretrained(pretrained_model)  # BertTokenizerFast
model = BertForSequenceClassification.from_pretrained(pretrained_model, num_labels=2)  # BertModel
print(model)
model.to(device)
model.classifier.to(device)

# Set number of batches (per epoch) and epochs
num_epochs = 20 
batch_size = int(best_hyperparams.batch_size.values[0])  # nb of example per batch
peak_lr = best_hyperparams.learning_rate.values[0]
#peak_lr = 0.00004
fraction = 0.1  # end factor of LR (e.g. 4e-5*0.25=1e-5)

# Define optimizer and loss function
optimizer = torch.optim.AdamW(model.parameters(), lr=peak_lr)
loss_fn = torch.nn.CrossEntropyLoss(weight=torch.tensor([1/5, 1]).to(device))


# Train over given fold (fold_num) in stratified 10-fold CV
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=rs)
fold_dict = {str(fold_index+1): [train_index, test_index]
            for fold_index, (train_index, test_index) in
            enumerate(skf.split(kmer_seqs, y_simple))}

train_index = fold_dict[fold_num][0]
test_index = fold_dict[fold_num][1]
print(len(train_index), len(test_index))
for k,v in fold_dict.items():
    print(k, len(v[0]), len(v[1]))
print(len(kmer_seqs))

print()

# Load train and eval (test) datasets
x_train = [k for k in kmer_seqs if kmer_seqs.index(k) in train_index]
x_test = [k for k in kmer_seqs if kmer_seqs.index(k) in test_index]
Y_train = [y_simple.loc[i, 'target'] for i in train_index]
Y_test = [y_simple.loc[i, 'target'] for i in test_index]
print(len(x_train), len(x_test), len(Y_train), len(Y_test))
print(min(train_index), max(train_index))
print(kmer_seqs.index(kmer_seqs[min(train_index)]), kmer_seqs.index(kmer_seqs[max(train_index)]))



# Load input sequences in right format (tokenize it for BERT)
inputs = tokenizer(x_train, return_tensors='pt', padding=True).to(device)
labels = torch.tensor(Y_train).to(device)
print(f'echo input_ids {inputs.input_ids.size()}')
print(f'echo att_mask {inputs.attention_mask.size()}')
print(f'echo labels {labels.size()}')
dataset = TensorDataset(inputs.input_ids, inputs.attention_mask, labels)
train_dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)


# Decrease linearly from peak value 4e-5 to 1e-5 over num_epochs
total_steps = len(train_dataloader) * num_epochs
scheduler = lr_scheduler.LinearLR(optimizer, start_factor=1.0, end_factor=fraction, total_iters=num_epochs)





# Iterate over epochs and batches per epoch
epoch_f_scores, epoch_loss = [], []
for epoch in range(num_epochs):
    p = f'EPOCH {epoch}'
    sp.call("echo " + p, shell=True)
    total_loss = 0.0
    
    # Start training
    for i, batch in enumerate(train_dataloader):
        sp.call(f'echo BATCH Train {i}', shell=True)
        input_ids, attention_mask, batch_labels = batch
        model.train()
        input_ids, batch_labels = input_ids.to(device), batch_labels.to(device)
        attention_mask = attention_mask.to(device)
        outputs = model(input_ids, attention_mask=attention_mask, labels=batch_labels)
        loss = loss_fn(outputs.logits.to(device), batch_labels)
        total_loss += loss.item()

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    before_lr = optimizer.param_groups[0]["lr"]
    scheduler.step()  # update learning rate
    after_lr = optimizer.param_groups[0]["lr"]
    sp.call(f"echo Before LR update {before_lr} After LR update {after_lr}", shell=True)

    avg_loss = total_loss / len(train_dataloader)  # across batches
    epoch_loss.append(avg_loss)
    sp.call(f"echo Epoch {epoch + 1}/{num_epochs}, Loss: {avg_loss}", shell=True)


    ## Compute f1_score on held-out fold
    model.eval()  # no more dropout

    # First, load input sequences in right format (tokenize it for BERT)
    eval_dataset = tokenizer(x_test, return_tensors='pt', padding=True).to(device)
    eval_labels = torch.tensor(Y_test).to(device)
    eval_dataset = TensorDataset(eval_dataset.input_ids, eval_dataset.attention_mask, eval_labels)
    eval_dataloader = DataLoader(eval_dataset, batch_size=batch_size)

    # Evaluate model on eval dataset per batch for 1 epoch
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

            # Convert logits to probabilities via the softmax activation function (in 2nd dimension of output)
            probabilities = torch.softmax(output.logits, dim=1)

            # Get the predicted labels from these probabilities of each class
            pred_labels = torch.argmax(probabilities, dim=1)  # get the index (0 or 1) of the highest probability
            ev_preds += pred_labels.tolist()
            ev_labels += ev_batch_labels.tolist()

    # Compute F1 score
    # The 'macro' makes it that each class (0 or 1) are taken into account equally (which is what we want)
    fscore = f1_score(ev_labels, ev_preds, average='macro')
    sp.call(f'echo Fscore: {fscore}', shell=True)

    # Save that f1_score for each epoch
    epoch_f_scores.append(fscore)
    if epoch == num_epochs - 1:
        sp.call('nvidia-smi', shell=True)


# Save those epoch f1_scores and loss for each fold
f_df = pd.DataFrame(epoch_f_scores, columns=['f1_score_per_epoch'])
f_df.to_csv(output_f1, sep='\t', index=False)
l_df = pd.DataFrame(epoch_loss, columns=['loss_per_epoch'])
l_df.to_csv(output_loss, sep='\t', index=False)




# Save model for that given fold (only save weights and parameters as it is lighter than saving the whole model)
model.to('cpu')
torch.save(model.state_dict(), output_model)

