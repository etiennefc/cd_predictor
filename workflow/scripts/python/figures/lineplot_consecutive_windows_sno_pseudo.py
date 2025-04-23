#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft
import collections as coll

df_pred = pd.read_csv(snakemake.input.preds, sep='\t')
color_dict = snakemake.params.color_conf_val

# Add confusion value columns based on pseudogene class
df_pred.loc[(df_pred['y_true'] == 0) & (df_pred['y_pred'] == 0), 'conf_val_pseudo'] = 'TP'
df_pred.loc[(df_pred['y_true'] == 0) & (df_pred['y_pred'] == 1), 'conf_val_pseudo'] = 'FN'
df_pred.loc[(df_pred['y_true'] == 1) & (df_pred['y_pred'] == 0), 'conf_val_pseudo'] = 'FP'
df_pred.loc[(df_pred['y_true'] == 1) & (df_pred['y_pred'] == 1), 'conf_val_pseudo'] = 'TN'

# Create simplified gene_id and suffix columns
df_pred['base_gene_id'] = df_pred['gene_id'].str.replace(r'(_[AB]\d+)?$', '', regex=True)
df_pred['suffix'] = df_pred['gene_id'].str.extract(r'_([AB]\d+)?$', expand=False).fillna('0')

# Create df copy to keep the prediction at the center window only for each example
df_pred2 = df_pred.copy()
df_pred2 = df_pred2[df_pred2['base_gene_id'] == df_pred2['gene_id']]


# Pivot df to have the prediction for each window per example
pivot_df = df_pred.pivot(index='base_gene_id', columns='suffix', values='y_pred').reset_index()

# Merge with the confusion value of the prediction of the center window
pivot_df = pivot_df.merge(df_pred2[['base_gene_id', 'conf_val_pseudo']], how='left', on='base_gene_id')
cols_ = [f'B{i}' for i in range(15, 0, -1)] + ['0'] + [f'A{i}' for i in range(1, 15+1)]
cols_2 = [f'B{i}' for i in range(5, 0, -1)] + ['0'] + [f'A{i}' for i in range(1, 5+1)]
pivot_df = pivot_df.drop_duplicates(['base_gene_id'] + cols_)



# Iterate over the type of predictions (TP/TN/FP/FN)
col_dict = {'TP': cols_, 'FN': cols_, 'TN': cols_2, 'FP': cols_2}
for v in ['TP', 'TN', 'FP', 'FN']:
	print(color_dict[v])
	output_ = [path for path in snakemake.output if v in path][0]
	conf_val_df = pivot_df[pivot_df['conf_val_pseudo'] == v]
	chosen_col = col_dict[v]

	# Groupby/count the number of example that have the same pattern of 
	# prediction across windows (ex: all 0, all 1, 0/1 alternance, etc.)
	overlap_df = conf_val_df.groupby(chosen_col).size().reset_index(name='overlap_count')

	# Regulate the alpha (transparency) based on the nb of examples that have the same prediction 
	# pattern across windows (the more common example, the more opaque the line will be)
	overlap_df['alpha'] = overlap_df['overlap_count'] / overlap_df['overlap_count'].sum() + 0.05
	overlap_df['fake_hue'] = 'test'  # fake_hue column needed for the plot generation


	fig, ax = plt.subplots(1, 1, figsize=(50, 14))
	plt.rcParams['svg.fonttype'] = 'none'

	for i, row in overlap_df.iterrows():
	    alph = row['alpha']
	    lw = 50 * alph
	    plt.xticks(fontsize=40)
	    plt.yticks(fontsize=40) 
	    pd.plotting.parallel_coordinates(overlap_df.iloc[[i]], 'fake_hue', chosen_col, 
					axvlines=False, ax=ax, alpha=alph, lw=lw, color=color_dict[v])
	    plt.legend().set_visible(False)
	plt.title(f'SnoBIRD 2nd model prediction distribution across windows of ({len(conf_val_df)}) {v} test set examples w/r to pseudogenes', 
				fontsize=40)
	plt.xlabel('Window position around example', fontsize=40)
	plt.ylabel('Prediction output', fontsize=40)
	plt.savefig(output_)
