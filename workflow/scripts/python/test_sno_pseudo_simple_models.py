#!/usr/bin/python3
import pandas as pd
from sklearn.metrics import f1_score, accuracy_score, precision_score, recall_score
import pickle

""" Test each model performance on unseen test data and report their accuracy."""

X_test = pd.read_csv(snakemake.input.X_test_scaled, sep='\t', index_col='gene_id')
y_test = pd.read_csv(snakemake.input.y_test, sep='\t')
y_test = y_test[y_test.gene_id.isin(X_test.index)]
y_test_copy = y_test.copy()
y_test = y_test.drop(columns=['gene_id'])
transformer_preds = pd.read_csv(snakemake.input.transformer_CD_preds, sep='\t')
cd_preds = transformer_preds[transformer_preds['y_pred'] == 1]

# Unpickle and thus instantiate the trained model 
model = pickle.load(open(snakemake.input.model, 'rb'))

# Predict target label (1 (CD_pseudogene) or 2 (expressed_CD)) on whole test set
y_pred = model.predict(X_test)
acc = {}
acc['model'] = snakemake.wildcards.simple_models
acc['test_accuracy'] = accuracy_score(y_test.target, y_pred)
acc['precision'] = precision_score(y_test.target, y_pred, average='macro')
acc['recall'] = recall_score(y_test.target, y_pred, average='macro')
acc['test_f1_score'] = f1_score(y_test.target, y_pred, average='macro')
acc_df = pd.DataFrame(acc, index=[0])
acc_df.to_csv(snakemake.output.df_metrics_on_test_all, sep='\t', index=False)

# Save predictions on whole test set to file
target_pred_col = str(snakemake.wildcards.simple_models)+'_prediction'
y_pred_df = pd.DataFrame(y_pred, columns=[target_pred_col])
preds = pd.concat([y_test_copy, y_pred_df], axis=1)
preds.to_csv(snakemake.output.all_test_predictions, sep='\t', index=False)

# Filter to keep only the examples that this predictor would receive, 
# i.e. those predicted as CD (in the general sense) by the transformer
#X_test_filtered = X_test[X_test.index.isin()]


