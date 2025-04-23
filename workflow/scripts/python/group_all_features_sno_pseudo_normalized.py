#!/usr/bin/python3
import pandas as pd
import collections as coll

target_output = [snakemake.output.target_tuning, snakemake.output.target_training, snakemake.output.target_test]
scaled_df_output = [snakemake.output.sno_pseudo_tuning, snakemake.output.sno_pseudo_training, snakemake.output.sno_pseudo_test]
# Load the already split examples into respectives sets for the first snoBIRD model
# Then keep only expressed CD and snoRNA pseudogene examples and no data augmented examples
tuning = pd.read_csv(snakemake.input.tuning, sep='\t')
tuning = tuning[tuning['target'] != 'other']
tuning = tuning[~tuning['gene_id'].str.contains('_A[0-9]$|_B[0-9]$')]
training = pd.read_csv(snakemake.input.training, sep='\t')
training = training[training['target'] != 'other']
training = training[~training['gene_id'].str.contains('_A[0-9]$|_B[0-9]$')]
test = pd.read_csv(snakemake.input.test, sep='\t')
test = test[test['target'] != 'other']
test = test[~test['gene_id'].str.contains('_A[0-9]$|_B[0-9]$')]

# Load features df
mfe = pd.read_csv(snakemake.input.structure_mfe, sep='\t')
terminal_stem = pd.read_csv(snakemake.input.terminal_stem, sep='\t')
box_score = pd.read_csv(snakemake.input.box_score, sep='\t')
all_features = box_score.merge(mfe[['gene_id', 'sno_stability', 'normalized_sno_stability',
                                    ]], how='left', on='gene_id')
all_features = all_features.merge(terminal_stem[['gene_id', 'terminal_stem_stability', 
                    'terminal_stem_length_score']], how='left', on='gene_id')
# Add dataset in which the examples are from
all_features.loc[all_features['gene_id'].isin(tuning['gene_id']), 'dataset'] = 'Tuning'
all_features.loc[all_features['gene_id'].isin(training['gene_id']), 'dataset'] = 'Training'
all_features.loc[all_features['gene_id'].isin(test['gene_id']), 'dataset'] = 'Test'

# Add predicted sno length
print(all_features.columns)
all_features['predicted_sno_len'] = all_features['predicted_sequence'].apply(lambda x: len(x))
all_features.to_csv(snakemake.output.all_examples, sep='\t', index=False)
relevant_features = ['score_c', 'score_d', 'score_c_prime', 'score_d_prime', 
                    'predicted_sno_len', 'sno_stability', 'terminal_stem_stability', 
                    'terminal_stem_length_score']

# Merge relevant features to the 3 sets
# Due to duplicate removal, some sno are present only in their augmented form (_A1 or _B1)
# 2 snoRNAs are like that and are therefore removed from the set (['Oa1850b_O_ana_Schmitz_2008', 'ENSG00000238519'])
tuning = tuning.merge(all_features[relevant_features + ['gene_id']], on='gene_id', how='left')
training = training.merge(all_features[relevant_features + ['gene_id']], on='gene_id', how='left')
test = test.merge(all_features[relevant_features + ['gene_id']], on='gene_id', how='left')



# Normalize (standardize) the datasets separately
norm_df = []
for i, df in enumerate([tuning, training, test]):
    target_df = df.copy()
    target_df['target'] = target_df['target'].replace({'snoRNA_pseudogene':0, 'expressed_CD_snoRNA': 1})
    target_df = target_df[['gene_id', 'target']]
    target_df.to_csv(target_output[i], sep='\t', index=False)


    # Scale/normalize the added features only using standardization
    df_num = df.copy()
    cols = ['score_c', 'score_d', 'score_c_prime', 'score_d_prime', 'sno_stability', 
            'terminal_stem_stability', 'terminal_stem_length_score', 'predicted_sno_len']
    df_num = df_num[cols]
    col_values = [df[['gene_id', 'target']]]
    for col in cols:
        mean = df_num[col].mean()
        std = df_num[col].std()
        if std != 0:
            col_values.append(pd.DataFrame((df[col] - mean) / std).rename(
                                columns={col: f'{col}_norm'}))
        else:  # to deal with column that has all the same value, thus a std=0
            col_values.append(pd.DataFrame(df[col]).rename(  # we don't scale, but these values will either be all 0 or all 1
                            columns={col: f'{col}_norm'}))  
    scaled_df = pd.concat(col_values, axis=1)
    norm_df.append(scaled_df)
    scaled_df.to_csv(scaled_df_output[i], index=False, sep='\t')  # save scaled output

all_norm_df = pd.concat(norm_df)
all_norm_df.to_csv(snakemake.output.all_examples_norm, sep='\t', index=False)



