#!/usr/bin/python3
import pandas as pd
import random
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, accuracy_score, precision_recall_fscore_support, confusion_matrix, classification_report
from sklearn.svm import OneClassSVM
import collections as coll
from sklearn.ensemble import IsolationForest

# Outputs
df_metrics = snakemake.output.df_metrics_on_test
df_preds = snakemake.output.test_predictions

# Load training data and create tensor matrix 
# Load SHAP values of first snoBIRD model and keep only TP
X = pd.read_csv(snakemake.input.shap_first_snobird, sep='\t')
shap_cols = [i for i in X.columns if ('SHAP' in i) | (i in ['CLS', 'SEP'])]
print(type(shap_cols))
#X = X.merge(cd_df[['gene_id', 'dataset', 'gene_biotype']], on='gene_id', how='left')
#X = X[(X['predicted_label'] == 'CD_snoRNA') & (X['dataset'] == 'Test')]
#y = X.copy()
#y['target'] = y['gene_biotype'].replace({'snoRNA_pseudogene': 0, 'expressed_CD_snoRNA': 1})
#
## Create tensor matrix (format used by pytorch)
#col_removal = ['gene_id', 'dataset', 'gene_biotype', 'predicted_label', 'probability']

all_features = ['score_c', 'score_d', 'box_score', 'score_c_prime', 'score_d_prime', 
            'len', 'sno_stability', 'normalized_sno_stability', 
            'terminal_stem_stability', 'terminal_stem_length_score']
features = ['box_score', 'normalized_sno_stability', 'terminal_stem_stability', 'terminal_stem_length_score']
all_cd = pd.read_csv(snakemake.input.all_examples, sep='\t')
all_cd['target'] = all_cd['gene_biotype'].replace({'snoRNA_pseudogene': 0, 'expressed_CD_snoRNA': 1})
all_cd = all_cd.merge(X[shap_cols + ['gene_id', 'probability']], on='gene_id', how='left')
print(all_cd)
print(list(all_cd.columns))

train = all_cd[(all_cd['gene_biotype'] == 'expressed_CD_snoRNA') & (all_cd['dataset'] == 'Training')][features + ['probability']]
test = all_cd[all_cd['dataset'] == 'Test'][features + ['probability']]
test_labels = all_cd[all_cd['dataset'] == 'Test']['target']
print(test)
print(test_labels)

clf = OneClassSVM(gamma='auto', kernel='poly', degree=4, nu=0.015)
clf.fit(train)
test_predictions = clf.predict(test)
print(test_predictions)
test_predictions[test_predictions == -1] = 0
print(test_predictions)
tn, fp, fn, tp = confusion_matrix(test_labels, test_predictions).ravel()
print(tn, fp, fn, tp)
print(classification_report(test_labels, test_predictions))

#clf = IsolationForest(contamination=0.2, random_state=42)
#clf.fit(train)
#test_predictions = clf.predict(test)
#print(test_predictions)
#test_predictions[test_predictions == -1] = 0
#print(test_predictions)
#tn, fp, fn, tp = confusion_matrix(test_labels, test_predictions).ravel()
#print(tn, fp, fn, tp)
#print(classification_report(test_labels, test_predictions))