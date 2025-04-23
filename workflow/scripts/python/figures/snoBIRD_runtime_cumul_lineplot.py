#!/usr/bin/python3
import pandas as pd
import os
import functions as ft
import subprocess as sp

output = snakemake.output.lineplot
log_dir = snakemake.input.logs
color_dict_all = snakemake.params.color_dict
color_dict = {}
for s in ['tetrahymena_thermophila', 'drosophila_melanogaster', 
        'homo_sapiens_chr1', 'gallus_gallus', 'macaca_mulatta']:
        if 'homo_sapiens' in s:
            color = color_dict_all['homo_sapiens']
        else:
            color = color_dict_all[s]
        color_dict[s] = color


months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 
        'Nov', 'Dec']
_31_days = [1, 3, 5, 7, 8, 10, 12]

# GPU type on which snoBIRD has been run depending on the species
gpu_dict = {'gallus_gallus': 'V100', 'macaca_mulatta': 'V100', 
            'drosophila_melanogaster': 'V100', 
            'tetrahymena_thermophila': 'V100', 
            'homo_sapiens_chr1': 'V100'}


def get_days_month(month_num):
    # Get the total number of days per month depending on the month
    if month_num in _31_days:
        return 31
    else:
        return 30
month_dict = {m: i+1 for i, m in enumerate(months)}

# Iterate over each log file to find rule name, start and end time
rule_time = {}
for sp_ in gpu_dict.keys():
    rule_time_sp = {}
    sp_dir = [f for f in log_dir if sp_ in f][0]
    time_genome_pred = 0
    for file_ in os.listdir(sp_dir):
        start = sp.run(['grep', '-B1', 'rule ', f'{sp_dir}/{file_}'], 
                        capture_output=True, text=True).stdout.strip()
        start = start.replace('  ', ' ')
        rule_name = start.split('rule ')[-1].replace(':', '')
        end = sp.run(['grep', '-B1', 'Finished job', f'{sp_dir}/{file_}'], 
                        capture_output=True, text=True).stdout.strip()
        end = end.replace('  ', ' ')
        if rule_name != 'all':
            # Get start and end time
            weekday_s, month_s, day_s, time_s, year_s = start.split(
                                                ']')[0].replace('[', '').split(' ')
            weekday_e, month_e, day_e, time_e, year_e = end.split(
                                                ']')[0].replace('[', '').split(' ')
            month_s, month_e = month_dict[month_s], month_dict[month_e]
            day_s, day_e = int(day_s), int(day_e)
            if month_s <= month_e:
                day_diff = (day_e - day_s) * 24 * 60  # in min
            else:
                day_diff = (day_e + (get_days_month(month_s) - day_s)) * 24 * 60  # in min
            hour_s, min_s, sec_s = map(int, time_s.split(':'))
            hour_e, min_e, sec_e = map(int, time_e.split(':'))
            min_s = min_s + (sec_s / 60)
            min_e = min_e + (sec_e / 60)
            hour_diff = (hour_e - hour_s) * 60  # in min
            min_diff = min_e - min_s
            time_diff = day_diff + hour_diff + min_diff
            if rule_name != 'genome_prediction':
                rule_time_sp[rule_name] = time_diff
            else:  # genome_prediciton is parallelized, so take only the longest time
                if time_diff > time_genome_pred:
                    time_genome_pred = time_diff 
                    rule_time_sp[rule_name] = time_genome_pred
    rule_time[sp_] = rule_time_sp



print(rule_time)

#total_time = sum([v for k,v in rule_time.items()])
# Get time values per rule in the right order
l = []
rules = ['[start]', 'split_chr', 'genome_prediction', 'merge_filter_windows', 
        'sno_pseudo_prediction', 'shap_snoBIRD', 'find_sno_limits_shap', 
        'filter_sno_pseudo_predictions_with_features']
for sp_, dictio in rule_time.items():
    vals = ['snoBIRD', sp_]
    for r in rules:
        if r in dictio.keys():
            vals.append(dictio[r])
        else:
            vals.append(0)
    l.append(vals)


df = pd.DataFrame(l, columns=['tool', 'species_genome']+rules)
df['start'] = 0
pivot_df = pd.melt(df, id_vars=['species_genome'], value_vars=rules, 
            var_name='rule_name', value_name='time')
print(df)
print(pivot_df)

ordered_species = list(color_dict.keys())
# Filter by species and rule name
final_dfs = [pd.DataFrame()] * len(ordered_species)
for i, group in enumerate(pivot_df.groupby('species_genome')):
    species_name = group[0]
    temp_df = group[1]
    order_ = ordered_species.index(species_name)
    total_time = float(group[1]['time'].sum())
    temp_df['cumsum'] = temp_df['time'].cumsum()
    temp_df['time_proportion'] = temp_df['cumsum'] / total_time * 100
    final_dfs[order_] = temp_df

final_df = pd.concat(final_dfs)
#final_df = final_df[final_df['species_genome'] == 'macaca_mulatta']
print(final_df)
print(final_df.columns)

# Create lineplot
ft.lineplot(final_df, 'rule_name', 'time_proportion', 'species_genome', 
    'Rule order in SnoBIRD pipeline', 'Proportion of total runtime (%)', '', 
    color_dict, output)


#df.to_csv(output_runtime, sep='\t', index=False)
