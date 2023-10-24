#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft 
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, recall_score, matthews_corrcoef

dictio = {0: "other", 1: "snoRNA_pseudogene", 2: "expressed_CD_snoRNA"}
snoreport = pd.read_csv(snakemake.input.snoreport, sep='\t')
snoscan = pd.read_csv(snakemake.input.snoscan, sep='\t')
infernal_rfam = pd.read_csv(snakemake.input.infernal_rfam, sep='\t')
#pseudosno_preds = pd.read_csv(snakemake.input.pseudosno_preds, sep='\t')
#simple_models_preds = pd.concat([pd.read_csv(snakemake.input.simple_models_preds[0], sep='\t')[['gene_id', 'target']]] + 
#                                [pd.read_csv(path, sep='\t').drop(columns=['target', 'gene_id']) for path in 
#                                snakemake.input.simple_models_preds], axis=1).replace(dictio)

# Best transformer yet 
metrics_transformer = [['accuracy', 0.9859, 'CD_predictor'], 
                ['precision', 0.9327, 'CD_predictor'], 
                ['recall', 0.8818, 'CD_predictor'],
                ['f1_score', 0.9065, 'CD_predictor']] 
                #['specificity on\nsnoRNA_pseudogenes', None, 'CD_predictor']]
transformer = pd.DataFrame(metrics_transformer, columns=['score_name', 'score_value', 'predictor'])


# These pseudosno predictions df have different lengths because they have a different number of FP
#logreg = pd.read_csv([path for path in snakemake.input.simple_models_pseudosno_preds 
#                    if 'logreg' in path][0], sep='\t').replace(dictio)
#svc = pd.read_csv([path for path in snakemake.input.simple_models_pseudosno_preds 
#                    if 'svc' in path][0], sep='\t').replace(dictio)   
#rf = pd.read_csv([path for path in snakemake.input.simple_models_pseudosno_preds 
#                    if '/rf' in path][0], sep='\t').replace(dictio)                                     
#gbm = pd.read_csv([path for path in snakemake.input.simple_models_pseudosno_preds 
#                    if 'gbm' in path][0], sep='\t').replace(dictio)
#knn = pd.read_csv([path for path in snakemake.input.simple_models_pseudosno_preds 
#                    if 'knn' in path][0], sep='\t').replace(dictio)                    
#pseudo_df_dict = {'logreg': logreg, 'svc': svc, 'rf': rf, 'gbm': gbm, 'knn': knn}                

output = snakemake.output.dotplot
#output_simple_models = snakemake.output.dotplot_simple_models
color_dict = snakemake.params.predictors_colors
#simple_models_color_dict = snakemake.params.simple_models_colors

# Merge dfs
df = snoreport.merge(snoscan[['gene_id', 'snoscan_prediction']], 
                        how='left', on='gene_id')
df = df.merge(infernal_rfam[['gene_id', 'infernal_rfam_prediction']], 
                        how='left', on='gene_id')


# Keep only the prediction made on positives and "other" negatives, not snoRNA pseudogenes (i.e. a separate class to predict)
df = df[df['gene_biotype'] != 'snoRNA_pseudogene']
#simple_models_preds = simple_models_preds[simple_models_preds['target'] != 'snoRNA_pseudogene']

# Compute the different metrics for each existing CD predictor
dfs = []
for model in color_dict.keys():
    # Based here on the existing_CD_snoRNA and other classes
    accuracy = accuracy_score(df['target'], df[f'{model}_prediction'])
    precision, recall, f1_score, _ = precision_recall_fscore_support(
                                    df['target'], df[f'{model}_prediction'], 
                                    pos_label='expressed_CD_snoRNA', average='binary')
    mcc = matthews_corrcoef(df['target'], df[f'{model}_prediction'])
    # Specificity or true negative rate is the same as the recall but computed based on the negative label
    # Since the existing predictors were not designed to predict snoRNA pseudogenes, we compute only their specificity 
    # (i.e. whether they're able to not predict it as an expressed snoRNA)
    #specificity =  recall_score(pseudosno_preds.target, pseudosno_preds[f'{model}_prediction'], pos_label="other")
    temp_df = pd.DataFrame([['accuracy', accuracy, model], 
                            ['precision', precision, model], 
                            ['recall', recall, model], 
                            ['f1_score', f1_score, model]],
                            #['specificity on\nsnoRNA_pseudogenes', specificity, model]], 
                            columns=['score_name', 'score_value', 'predictor'])
    print(model, ' on expressed_CD_snoRNA/other in test set')
    print('accuracy: ', accuracy)
    print('precision: ', precision)
    print('recall: ', recall)
    print('f1_score: ', f1_score)
    #print('specificity on\nsnoRNA pseudogenes: ', specificity, '\n')
    dfs.append(temp_df)

final_df = pd.concat(dfs + [transformer]).reset_index(drop=True)
print(final_df)
color_dict = {'snoreport2': '#fda47a', 'snoscan': '#80b1d3', 'infernal_rfam': '#d73027', 'CD_predictor': '#66a61e'}

# Create graph
ft.lineplot(final_df, 'score_name', 'score_value', 'predictor', 
            'Metrics', 'Metrics value', '', color_dict, output, markersize=15, linewidth=8)

## For simple models (which were designed to predict the 3 classes)
# Compute the accuracy, precision, recall and f1_score based on the expressed_CD_snoRNA and other
#dfs = []
#for model in simple_models_color_dict.keys():
#    # Based here on the existing_CD_snoRNA as positives and other as negatives 
#    accuracy = accuracy_score(simple_models_preds['target'], simple_models_preds[f'{model}_prediction'])
#    print(pd.unique(simple_models_preds['target']))
#    print(pd.unique(simple_models_preds[f'{model}_prediction']))
#    temp = simple_models_preds.copy()
#    # Replace snoRNA_pseudogene predictions as "other" to keep it binary
#    temp[f'{model}_prediction'] = temp[f'{model}_prediction'].replace('snoRNA_pseudogene', 'other')
#    precision, recall, f1_score, _ = precision_recall_fscore_support(
#                                    simple_models_preds['target'], temp[f'{model}_prediction'], 
#                                    pos_label='expressed_CD_snoRNA', average='binary')
#
#    temp_df = pd.DataFrame([['accuracy', accuracy, model],
#                            ['precision', precision, model], 
#                            ['recall', recall, model], 
#                            ['f1_score', f1_score, model]], 
#                            columns=['score_name', 'score_value', 'predictor'])
#    print(model, ' on snoRNA pseudogenes in test set')
#    print('accuracy: ', accuracy)
#    print('precision: ', precision)
#    print('recall: ', recall)
#    print('f1_score: ', f1_score)
#    dfs.append(temp_df)
#
#final_df = pd.concat(dfs)
#
#
## Compute the precision, recall and f1_score based on the snoRNA_pseudogene class
#dfs_ = []
#for model in pseudo_df_dict.keys():
#    # Based here on the snoRNA_pseudogene as TP
#    model_df = pseudo_df_dict[model]
#    # Replace expressed_CD_snoRNA predictions as "other" to keep it binary
#    model_df[f'{model}_prediction'] = model_df[f'{model}_prediction'].replace('expressed_CD_snoRNA', 'other')
#    model_df['target'] = model_df['target'].replace('expressed_CD_snoRNA', 'other')
#
#    precision, recall, f1_score, _ = precision_recall_fscore_support(
#                                    model_df['target'], model_df[f'{model}_prediction'], 
#                                    pos_label='snoRNA_pseudogene', average='binary')
#
#    temp_df = pd.DataFrame([['precision_on_pseudosno', precision, model], 
#                            ['recall_on_pseudosno', recall, model], 
#                            ['f1_score_on_pseudosno', f1_score, model]], 
#                            columns=['score_name', 'score_value', 'predictor'])
#    print(model, ' on snoRNA pseudogenes in test set')
#    print('precision: ', precision)
#    print('recall: ', recall)
#    print('f1_score: ', f1_score)
#    dfs_.append(temp_df)
#
#
#final_df_pseudosno = pd.concat([final_df] + dfs_)
#ft.lineplot(final_df_pseudosno, 'score_name', 'score_value', 'predictor', 
#            'Metrics', 'Metrics value', '', simple_models_color_dict, output_simple_models)
#