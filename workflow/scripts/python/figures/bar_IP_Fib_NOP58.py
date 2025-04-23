#!/usr/bin/python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats

# Load dfs
fib = pd.read_csv(snakemake.input.fib, sep='\t')
nop58 = pd.read_csv(snakemake.input.nop58, sep='\t')

print(fib)

results_fib, results_nop58 = [], []

# Perform unpaired t-test for each gene for Fib IP
for _, row in fib.iterrows():
    gene = row["gene_id"]
    control_values = list(row.filter(like="Flag").values)
    cond1_values = list(row.filter(like="Fibrillarin").values)
    # Test if normally distributed:
    stat_control, p_control = stats.shapiro(control_values)
    stat_cond1, p_cond1 = stats.shapiro(cond1_values)
    # Test if variance are significantly different
    if (p_control < 0.05) | (p_cond1 < 0.05):  #not normal
        var_stat, var_p_value = stats.levene(control_values, cond1_values)
        if var_p_value < 0.05:  # different variances
            # Independent t-test (assumes samples are from different distributions)
            t_stat, p_value = stats.ttest_ind(control_values, cond1_values, equal_var=False)
        else: # variance are not different
            t_stat, p_value = stats.ttest_ind(control_values, cond1_values, equal_var=True)
      
    else: # both normally distributed
        var_stat, var_p_value = stats.bartlett(control_values, cond1_values)
        if var_p_value < 0.05:  # different variances
            # Independent t-test (assumes samples are from different distributions)
            t_stat, p_value = stats.ttest_ind(control_values, cond1_values, equal_var=False)
        else: # variance are not different
            t_stat, p_value = stats.ttest_ind(control_values, cond1_values, equal_var=True)
    
    
    results_fib.append({"gene": gene, "t_stat": t_stat, "p_value": p_value})

df_stats_fib = pd.DataFrame(results_fib)
print(df_stats_fib)

# Perform unpaired t-test for each gene for NOP58 IP
for _, row in nop58.iterrows():
    gene = row["gene_id"]
    control_values = list(row.filter(like="Untag").values)
    cond1_values = list(row.filter(like="NOP58-TAP").values)
    # Test if normally distributed:
    stat_control, p_control = stats.shapiro(control_values)
    stat_cond1, p_cond1 = stats.shapiro(cond1_values)
    # Test if variance are significantly different
    if (p_control < 0.05) | (p_cond1 < 0.05):  #not normal
        var_stat, var_p_value = stats.levene(control_values, cond1_values)
        if var_p_value < 0.05:  # different variances
            # Independent t-test (assumes samples are from different distributions)
            t_stat, p_value = stats.ttest_ind(control_values, cond1_values, equal_var=False)
        else: # variance are not different
            t_stat, p_value = stats.ttest_ind(control_values, cond1_values, equal_var=True)
      
    else: # both normally distributed
        var_stat, var_p_value = stats.bartlett(control_values, cond1_values)
        if var_p_value < 0.05:  # different variances
            # Independent t-test (assumes samples are from different distributions)
            t_stat, p_value = stats.ttest_ind(control_values, cond1_values, equal_var=False)
        else: # variance are not different
            t_stat, p_value = stats.ttest_ind(control_values, cond1_values, equal_var=True)
    
    
    results_nop58.append({"gene": gene, "t_stat": t_stat, "p_value": p_value})

df_stats_nop58 = pd.DataFrame(results_nop58)
print(nop58)
print(df_stats_nop58)


#Create bar plot for fib
#Melt df in right format for barplot
df_long = fib.melt(id_vars="gene_id", var_name="Fibrillarin", value_name="fold_enrichment")

# Extract condition names (control FLAG vs. Fib IP)
df_long["group"] = df_long["Fibrillarin"].str.extract(r"(Flag|Fibrillarin)")

# Create figure
plt.rcParams['svg.fonttype'] = 'none'
plt.figure(figsize=(20, 10))

# Barplot: Shows mean expression per gene
sns.barplot(data=df_long, x="gene_id", y="fold_enrichment", hue="group", errorbar="sd", 
        capsize=0.2, palette=["lightgray", "darkred"], alpha=0.6)

# Stripplot: Adds dots for individual replicates
sns.stripplot(data=df_long, x="gene_id", y="fold_enrichment", hue="group", dodge=True, 
            jitter=True, size=4, marker="o", edgecolor="black", linewidth=0.2, 
            palette=["black", "black"], legend=False)

# Adjust legend and title
plt.legend(loc="upper right", fontsize=25)
plt.text(1, 6, 'N=4', fontsize=25)
plt.xlabel("Gene", fontsize=30)
plt.ylabel("Fold enrichment", fontsize=30)
plt.xticks(fontsize=25, rotation=25, ha='right')
plt.yticks(fontsize=25)
plt.ylim(0, None)

# Save plot
plt.savefig(snakemake.output.bar_fib, bbox_inches='tight', dpi=500)


#Create bar plot for fib
#Melt df in right format for barplot
df_long = nop58.melt(id_vars="gene_id", var_name="NOP58-TAP", value_name="fold_enrichment")

# Extract condition names (control Untag vs. NOP58-TAP IP)
df_long["group"] = df_long["NOP58-TAP"].str.extract(r"(Untag|NOP58-TAP)")

# Create figure
plt.rcParams['svg.fonttype'] = 'none'
plt.figure(figsize=(20, 10))

# Barplot: Shows mean expression per gene
sns.barplot(data=df_long, x="gene_id", y="fold_enrichment", hue="group", errorbar="sd", 
            capsize=0.2, palette=["lightgray", "darkgreen"], alpha=0.6)

# Stripplot: Adds dots for individual replicates
sns.stripplot(data=df_long, x="gene_id", y="fold_enrichment", hue="group", dodge=True, 
        jitter=True, size=4, marker="o", edgecolor="black", linewidth=0.2, 
        palette=["black", "black"], legend=False)

# Adjust legend and title
plt.legend(loc="upper right", fontsize=25)
plt.text(1, 200, 'N=5', fontsize=25)
plt.xlabel("Gene", fontsize=30)
plt.ylabel("Fold enrichment", fontsize=30)
plt.xticks(fontsize=25, rotation=25, ha='right')
plt.yticks(fontsize=25)
plt.ylim(0, None)

# Save plot
plt.savefig(snakemake.output.bar_nop58, bbox_inches='tight', dpi=500)
