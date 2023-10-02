#!/usr/bin/python3
print('allo')
import pandas as pd
from sklearn.preprocessing import StandardScaler
print('test')
tune = pd.read_csv(snakemake.input.X_tuning, sep='\t', index_col='gene_id')
train = pd.read_csv(snakemake.input.X_train, sep='\t', index_col='gene_id')
test = pd.read_csv(snakemake.input.X_test, sep='\t', index_col='gene_id')
cols = ['box_score', 'structure_mfe', 'terminal_stem_mfe', 'length']

# Keep only the expressed C/D and C/D pseudogenes
tune = tune[tune['target'] != 'other']
train = train[train['target'] != 'other']
test = test[test['target'] != 'other']

# Create and fit the scaler on the training data only
scaler = StandardScaler()
scaler.fit(train[cols])

# Apply scaler to all the sets separately
tune_scaled = scaler.transform(tune[cols]).reset_index()
train_scaled = scaler.transform(train[cols]).reset_index()
test_scaled = scaler.transform(test[cols]).reset_index()
print(tune_scaled)
print(train_scaled)
print(test_scaled)

# Save dfs
#tune_scaled.to_csv(snakemake.output.X_tuning_scaled, index=False, sep='\t')
#train_scaled.to_csv(snakemake.output.X_train_scaled, index=False, sep='\t')
#test_scaled.to_csv(snakemake.output.X_test_scaled, index=False, sep='\t')
