#!/usr/bin/python3
import pandas as pd
from glob import glob
import subprocess as sp
import re

output_runtime = snakemake.output.runtime
log_snoreport_pombe = snakemake.input.log_snoreport_pombe
log_snoscan_pombe = snakemake.input.log_snoscan_pombe
log_infernal_pombe = snakemake.input.log_infernal_pombe
log_snoreport_human = snakemake.input.log_snoreport_human
log_snoscan_human = snakemake.input.log_snoscan_human
log_infernal_human = snakemake.input.log_infernal_human
dir_dict = {'schizosaccharomyces_pombe': {'snoreport2': log_snoreport_pombe, 
            'infernal_rfam': log_infernal_pombe, 'snoscan': log_snoscan_pombe}, 
            'homo_sapiens': {'snoreport2': log_snoreport_human, 
            'infernal_rfam': log_infernal_human, 'snoscan': log_snoscan_human}}
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 
        'Nov', 'Dec']
_31_days = [1, 3, 5, 7, 8, 10, 12]

def get_days_month(month_num):
    # Get the total number of days per month depending on the month
    if month_num in _31_days:
        return 31
    else:
        return 30
month_dict = {m: i+1 for i, m in enumerate(months)}

# Iterate over each log file to find rule name, start and end time
tool_time_pombe, tool_time_human = {}, {}
for species in ['schizosaccharomyces_pombe', 'homo_sapiens']:
    for tool in ['snoreport2', 'infernal_rfam', 'snoscan']:
        temp_time = 0
        for file_ in glob(dir_dict[species][tool]+'/*out'):
            start = sp.run(['grep', '-B1', 'rule ', file_], 
                            capture_output=True, text=True).stdout.strip()
            start = start.replace('  ', ' ')
            rule_name = start.split('rule ')[-1].replace(':', '')
            tool_name = rule_name.split('_')[-2]
            if tool_name == 'infernal':
                tool_name = 'infernal_rfam'
            end = sp.run(['grep', '-B1', 'Finished job', file_], 
                            capture_output=True, text=True).stdout.strip()
            end = end.replace('  ', ' ')
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
            # These tools are not parallelized, so we must add the time that 
            # they took on each chromosome
            temp_time += time_diff
        if species == 'schizosaccharomyces_pombe':
            tool_time_pombe[tool_name] = temp_time
        elif species == 'homo_sapiens':
            tool_time_human[tool_name] = temp_time

print(tool_time_pombe)
print(tool_time_human)            


df_val = []
for tool, total_time in tool_time_pombe.items():
    df_val.append([tool, 'schizosaccharomyces_pombe', total_time, 'Narval_2_cpu'])
for tool, total_time in tool_time_human.items():
    df_val.append([tool, 'homo_sapiens', total_time, 'Narval_2_cpu'])

df = pd.DataFrame(df_val, columns=[
                'tool', 'species_genome', 'runtime_min', 'gpu_type'])
df.to_csv(output_runtime, sep='\t', index=False)
