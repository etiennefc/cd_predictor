#!/usr/bin/python3
import pandas as pd
import collections as coll
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
from figures import functions as ft
warnings.simplefilter('ignore')

feat_cols = ['probability', 'box_score', 'terminal_combined', 'normalized_sno_stability', 
            'gene_biotype', 'old_predicted_label', 'new_predicted_label', 'new_conf_val']

# Load and merge dfs
shap_df = pd.read_csv(snakemake.input.shap_df, sep='\t')
shap_df['predicted_label'] = shap_df['second_model_prediction'].replace(
                        {0: 'snoRNA_pseudogene', 1: 'expressed_CD_snoRNA'})
shap_df['probability'] = shap_df['probability_expressed_pseudogene']
feat_df = pd.read_csv(snakemake.input.all_examples, sep='\t')
feat_df = feat_df.merge(shap_df[['gene_id', 'probability', 'predicted_label']], 
            on='gene_id', how='left')
cols = ['gene_biotype', 'new_predicted_label', 'terminal_combined', 'box_score', 
        'normalized_sno_stability', 'gene_id', 'score_c', 'score_d', 'probability', 
        'new_conf_val']

# Create terminal_combined column
feat_df['terminal_combined'] = feat_df.terminal_stem_stability * feat_df.terminal_stem_length_score
feat_df['normalized_sno_stability'] = feat_df['normalized_sno_stability'].round(2)
feat_df['terminal_combined'] = feat_df['terminal_combined'].round(2)
# Create initial conf value column (based on pseudogene class and the 2nd snoBIRD model output only) 
# and select examples in test set
feat_df.loc[(feat_df['gene_biotype'] == 'snoRNA_pseudogene') & (feat_df['predicted_label'] == 'snoRNA_pseudogene'), 
            'conf_val'] = 'TP'
feat_df.loc[(feat_df['gene_biotype'] != 'snoRNA_pseudogene') & (feat_df['predicted_label'] == 'snoRNA_pseudogene'), 
            'conf_val'] = 'FP'
feat_df.loc[(feat_df['gene_biotype'] == 'expressed_CD_snoRNA') & (feat_df['predicted_label'] == 'expressed_CD_snoRNA'), 
            'conf_val'] = 'TN'
feat_df.loc[(feat_df['gene_biotype'] == 'snoRNA_pseudogene') & (feat_df['predicted_label'] != 'snoRNA_pseudogene'), 
            'conf_val'] = 'FN'
test_df = feat_df[feat_df['dataset'] == 'Test']


# Show prediction metrics before applying any additional features
print('\n\n\nPrediction metrics BEFORE applying feature-based filters:')
tp = len(test_df[test_df['conf_val'] == 'TP'])
fp = len(test_df[test_df['conf_val'] == 'FP'])
tn = len(test_df[test_df['conf_val'] == 'TN'])
fn = len(test_df[test_df['conf_val'] == 'FN'])
prec_p, reca_p = tp/(tp+fp), tp/(tp+fn)
print(f'precision pseudo: {prec_p}', f'recall pseudo: {reca_p}', f'f1 pseudo: {2*prec_p*reca_p/(prec_p+reca_p)}')
prec_e, reca_e = tn/(tn+fn), tn/(tn+fp)
print(f'precision expressed CD: {prec_e}', f'recall expressed CD: {reca_e}', 
        f'f1 expressed CD: {2*prec_e*reca_e/(prec_e+reca_e)}\n\n\n')


# Apply 1st filter post-snoBIRD's 2nd model
# Filter based on prediction probability (no change to prediction if pred probability > 0.999)
test_df['new_conf_val'] = test_df['conf_val']
test_df['new_predicted_label'] = test_df['predicted_label']
change = test_df[test_df['probability'] <= 0.999]
nochange = test_df[test_df['probability'] > 0.999]
changed_id = []

# FIRST FILTER: on prediction with low probability (<=0.999)
# If C box has 2 mutations or more and D box has >=1 mutation --> pseudo
change.loc[(~change['gene_id'].isin(changed_id)) & (change['score_c'] >= 2) & (change['score_d'] >= 1), 'new_predicted_label'] = 'snoRNA_pseudogene'
changed_id.extend(list(change[(~change['gene_id'].isin(changed_id)) & (change['score_c'] >= 2) & (change['score_d'] >= 1)]['gene_id']))
change.loc[(change['new_predicted_label'] == 'expressed_CD_snoRNA') & (change['gene_biotype'] == 'expressed_CD_snoRNA'), 'new_conf_val'] = 'TN'
change.loc[(change['new_predicted_label'] == 'snoRNA_pseudogene') & (change['gene_biotype'] == 'snoRNA_pseudogene'), 'new_conf_val'] = 'TP'
change.loc[(change['new_predicted_label'] == 'expressed_CD_snoRNA') & (change['gene_biotype'] == 'snoRNA_pseudogene'), 'new_conf_val'] = 'FN'
change.loc[(change['new_predicted_label'] == 'snoRNA_pseudogene') & (change['gene_biotype'] == 'expressed_CD_snoRNA'), 'new_conf_val'] = 'FP'
print('After 1st filter (score_c >= 2 & score_d >= 1):')
print(coll.Counter(pd.concat([change['new_conf_val'], nochange['new_conf_val']])), '\n')


# SECOND FILTER: on prediction with low probability (<=0.999)
# If terminal combined < -25 --> expressed
change.loc[(~change['gene_id'].isin(changed_id)) & (change['terminal_combined'] <= -25), 'new_predicted_label'] = 'expressed_CD_snoRNA'
changed_id.extend(list(change[(~change['gene_id'].isin(changed_id)) & (change['terminal_combined'] <= -25)]['gene_id']))
change.loc[(change['new_predicted_label'] == 'expressed_CD_snoRNA') & (change['gene_biotype'] == 'expressed_CD_snoRNA'), 'new_conf_val'] = 'TN'
change.loc[(change['new_predicted_label'] == 'snoRNA_pseudogene') & (change['gene_biotype'] == 'snoRNA_pseudogene'), 'new_conf_val'] = 'TP'
change.loc[(change['new_predicted_label'] == 'expressed_CD_snoRNA') & (change['gene_biotype'] == 'snoRNA_pseudogene'), 'new_conf_val'] = 'FN'
change.loc[(change['new_predicted_label'] == 'snoRNA_pseudogene') & (change['gene_biotype'] == 'expressed_CD_snoRNA'), 'new_conf_val'] = 'FP'
print('After 2nd filter (terminal_combined <= -25):')
print(coll.Counter(pd.concat([change['new_conf_val'], nochange['new_conf_val']])), '\n')


# THIRD FILTER: on prediction with low probability (<=0.999)
# If less than 2 features are "good" --> pseudo
change.loc[(~change['gene_id'].isin(changed_id)) & (((change['terminal_combined'] <= -25).astype(int) + 
            (change['normalized_sno_stability'] <= -0.2).astype(int) + 
            (change['box_score'] <= 5).astype(int) + (change['score_c'] <= 2).astype(int) + (change['score_d'] < 1).astype(int)) <= 2), 'new_predicted_label'] = 'snoRNA_pseudogene'
changed_id.extend(list(change.loc[(~change['gene_id'].isin(changed_id)) & (((change['terminal_combined'] <= -25).astype(int) + 
            (change['normalized_sno_stability'] <= -0.2).astype(int) + 
            (change['box_score'] <= 5).astype(int) + (change['score_c'] <= 2).astype(int) + (change['score_d'] < 1).astype(int)) <= 2)]['gene_id']))
change.loc[(change['new_predicted_label'] == 'expressed_CD_snoRNA') & (change['gene_biotype'] == 'expressed_CD_snoRNA'), 'new_conf_val'] = 'TN'
change.loc[(change['new_predicted_label'] == 'snoRNA_pseudogene') & (change['gene_biotype'] == 'snoRNA_pseudogene'), 'new_conf_val'] = 'TP'
change.loc[(change['new_predicted_label'] == 'expressed_CD_snoRNA') & (change['gene_biotype'] == 'snoRNA_pseudogene'), 'new_conf_val'] = 'FN'
change.loc[(change['new_predicted_label'] == 'snoRNA_pseudogene') & (change['gene_biotype'] == 'expressed_CD_snoRNA'), 'new_conf_val'] = 'FP'
print('After 3rd filter: <= 2 "good" features --> snoRNA_pseudogene')
print(coll.Counter(pd.concat([change['new_conf_val'], nochange['new_conf_val']])), '\n')


# FOURTH FILTER: on prediction with high probability (>0.999) (only on snoRNA_pseudogene predictions)
# If 4 out of 5 features are "good" --> expressed_CD_snoRNA
nochange['last_predicted_label'] = nochange['new_predicted_label']
nochange.loc[(nochange['last_predicted_label'] == 'snoRNA_pseudogene') & (((nochange['terminal_combined'] <= -25).astype(int) + 
            (nochange['normalized_sno_stability'] <= -0.2).astype(int) + 
            (nochange['box_score'] <= 5).astype(int) + (nochange['score_c'] <= 2).astype(int) + (nochange['score_d'] < 1).astype(int)) >= 4), 'new_predicted_label'] = 'expressed_CD_snoRNA'

nochange.loc[(nochange['new_predicted_label'] == 'expressed_CD_snoRNA') & (nochange['gene_biotype'] == 'expressed_CD_snoRNA'), 'new_conf_val'] = 'TN'
nochange.loc[(nochange['new_predicted_label'] == 'snoRNA_pseudogene') & (nochange['gene_biotype'] == 'snoRNA_pseudogene'), 'new_conf_val'] = 'TP'
nochange.loc[(nochange['new_predicted_label'] == 'expressed_CD_snoRNA') & (nochange['gene_biotype'] == 'snoRNA_pseudogene'), 'new_conf_val'] = 'FN'
nochange.loc[(nochange['new_predicted_label'] == 'snoRNA_pseudogene') & (nochange['gene_biotype'] == 'expressed_CD_snoRNA'), 'new_conf_val'] = 'FP'
print('After 4th filter: >= 4 "good" features --> expressed_CD_snoRNA')
print(coll.Counter(pd.concat([change['new_conf_val'], nochange['new_conf_val']])), '\n')

# Concat nochange and change dfs
concat_df = pd.concat([nochange, change])


# LAST filter: if C and/or D motif are NNNNNNN or NNNN --> pseudogene
concat_df.loc[(concat_df['C_MOTIF'] == 'NNNNNNN') | (concat_df['D_MOTIF'] == 'NNNN'), 'new_predicted_label'] = 'snoRNA_pseudogene'

# Get new confusion value
concat_df.loc[(concat_df['new_predicted_label'] == 'expressed_CD_snoRNA') & (concat_df['gene_biotype'] == 'expressed_CD_snoRNA'), 'new_conf_val'] = 'TN'
concat_df.loc[(concat_df['new_predicted_label'] == 'snoRNA_pseudogene') & (concat_df['gene_biotype'] == 'snoRNA_pseudogene'), 'new_conf_val'] = 'TP'
concat_df.loc[(concat_df['new_predicted_label'] == 'expressed_CD_snoRNA') & (concat_df['gene_biotype'] == 'snoRNA_pseudogene'), 'new_conf_val'] = 'FN'
concat_df.loc[(concat_df['new_predicted_label'] == 'snoRNA_pseudogene') & (concat_df['gene_biotype'] == 'expressed_CD_snoRNA'), 'new_conf_val'] = 'FP'

print('After 5th filter: if C = NNNNNNN | D = NNNN --> snoRNA_pseudogene')
print(coll.Counter(concat_df['new_conf_val']), '\n')

# Compute test set metrics after applying filters
tp = len(concat_df[concat_df['new_conf_val'] == 'TP'])
fp = len(concat_df[concat_df['new_conf_val'] == 'FP'])
tn = len(concat_df[concat_df['new_conf_val'] == 'TN'])
fn = len(concat_df[concat_df['new_conf_val'] == 'FN'])


prec_p, reca_p = tp/(tp+fp), tp/(tp+fn)  # metrics for pseudo
f1_p = 2*prec_p*reca_p/(prec_p+reca_p)
print('\n\n\nPrediction metrics AFTER applying feature-based filters:')
print(f'precision pseudo: {prec_p}', f'recall pseudo: {reca_p}', f'f1 pseudo: {f1_p}')
prec_e, reca_e = tn/(tn+fn), tn/(tn+fp)  # metrics for expressed C/D
f1_e = 2*prec_e*reca_e/(prec_e+reca_e)
print(f'precision expressed CD: {prec_e}', f'recall expressed CD: {reca_e}', 
        f'f1 expressed CD: {f1_e}\n\n\n')
accuracy = (tp + tn) / (tp + tn + fp + fn)
metrics_df = pd.DataFrame([[accuracy, prec_e, reca_e, f1_e, prec_p, reca_p, f1_p]], 
            columns=['accuracy_test', 'precision_sno_test', 'recall_sno_test', 'f1_sno_test', 
                    'precision_pseudo_test', 'recall_pseudo_test', 'f1_pseudo_test'])


# Save dfs
concat_df['old_predicted_label'] = concat_df['predicted_label']
concat_df[['gene_id', 'species_name'] + feat_cols].to_csv(snakemake.output.preds, sep='\t', index=False)

metrics_df.to_csv(snakemake.output.metrics, sep='\t', index=False)


