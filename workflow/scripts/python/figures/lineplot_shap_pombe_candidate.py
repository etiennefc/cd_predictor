#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import functions as ft

shap_df = pd.read_csv(snakemake.input.shap_df, sep='\t')
shap_df = shap_df[shap_df['gene_id'] == 'CD_531'].drop(
            columns=['gene_id', 'predicted_label', 'probability'])

# Create lineplot
xticklabels = [i.split('_pos_')[1] if i not in ['CLS', 'SEP'] else i for i in shap_df.columns]
ft.shap_lineplot(range(len(xticklabels)), shap_df.values.tolist()[0], xticklabels, 'Sequence token position', 
                'Average SHAP values', 
                'Distribution of avg SHAP values over\nCD_531 candidate sequence in S. pombe', 
                snakemake.output.lineplot)