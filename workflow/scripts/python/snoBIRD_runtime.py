#!/usr/bin/python3
import pandas as pd
import os
import subprocess as sp
import re

output_runtime = snakemake.output.runtime
species = snakemake.wildcards.species
log_dir = snakemake.input.logs
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 
        'Nov', 'Dec']
_31_days = [1, 3, 5, 7, 8, 10, 12]

# GPU type on which snoBIRD has been run depending on the species
gpu_dict = {'schizosaccharomyces_pombe': 'A100', 'homo_sapiens': 'A100', 
            'gallus_gallus': 'V100', 'macaca_mulatta': 'V100', 
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
time_genome_pred = 0
for file_ in os.listdir(log_dir):
    start = sp.run(['grep', '-B1', 'rule ', f'{log_dir}/{file_}'], 
                    capture_output=True, text=True).stdout.strip()
    start = start.replace('  ', ' ')
    rule_name = start.split('rule ')[-1].replace(':', '')
    end = sp.run(['grep', '-B1', 'Finished job', f'{log_dir}/{file_}'], 
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
            rule_time[rule_name] = time_diff
        else:  # genome_prediciton is parallelized, so take only the longest time
            if time_diff > time_genome_pred:
                time_genome_pred = time_diff 
                rule_time[rule_name] = time_genome_pred

# SnoBIRD in human was tested only with step_size =1 on A100 GPUs (in order 
# to not run it another time uselessly with s=5), so divide genome_prediction 
# time by 5
if species == 'homo_sapiens': 
    rule_time['genome_prediction'] = rule_time['genome_prediction'] / 5

print(rule_time)

total_time = sum([v for k,v in rule_time.items()])
print(f'Total time SnoBIRD on {species} genome: {total_time} min')

df = pd.DataFrame([['snoBIRD', species, total_time, gpu_dict[species]]], columns=[
                'tool', 'species_genome', 'runtime_min', 'gpu_type'])
df.to_csv(output_runtime, sep='\t', index=False)
