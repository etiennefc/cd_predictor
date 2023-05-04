#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, recall_score


snoreport = pd.read_csv(snakemake.input.snoreport, sep='\t')
snoscan = pd.read_csv(snakemake.input.snoscan, sep='\t')
infernal_rfam = pd.read_csv(snakemake.input.infernal_rfam, sep='\t')
pseudosno_preds = pd.read_csv(snakemake.input.pseudosno_preds, sep='\t')
output = snakemake.output.dotplot
color_dict = snakemake.params.colors

# Merge dfs
df = snoreport.merge(snoscan[['gene_id', 'snoscan_prediction']], 
                        how='left', on='gene_id')
df = df.merge(infernal_rfam[['gene_id', 'infernal_rfam_prediction']], 
                        how='left', on='gene_id')

# Compute the different metrics for each model
dfs = []
for model in color_dict.keys():
    accuracy = accuracy_score(df['target'], df[f'{model}_prediction'])
    precision, recall, f1_score, _ = precision_recall_fscore_support(
                                    df['target'], df[f'{model}_prediction'], 
                                    pos_label='expressed_CD_snoRNA', average='binary')
    # Specificity or true negative rate is the same as the recall but computed based on the negative label
    specificity =  recall_score(pseudosno_preds.target, pseudosno_preds[f'{model}_prediction'], pos_label="other")
    temp_df = pd.DataFrame([['accuracy', accuracy, model], 
                            ['precision', precision, model], 
                            ['recall', recall, model], 
                            ['f1_score', f1_score, model],
                            ['specificity on\nsnoRNA_pseudogenes', specificity, model]], 
                            columns=['score_name', 'score_value', 'predictor'])
    print(model, ' on test set')
    print('accuracy: ', accuracy)
    print('precision: ', precision)
    print('recall: ', recall)
    print('f1_score: ', f1_score)
    print('specificity on\nsnoRNA pseudogenes: ', specificity, '\n')
    dfs.append(temp_df)

final_df = pd.concat(dfs)
# Create graph
ft.lineplot(final_df, 'score_name', 'score_value', 'predictor', 
            'Metrics', 'Metrics value', '', color_dict, output)