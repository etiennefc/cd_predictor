#!/usr/bin/python3
import pandas as pd
from glob import glob
import subprocess as sp
import re

output_runtime = snakemake.output.runtime
species = snakemake.wildcards.species
max_time = snakemake.params.max_time
log_snoreport = snakemake.input.log_snoreport
log_snoscan = snakemake.input.log_snoscan
log_infernal = snakemake.input.log_infernal
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
tool_time = {}
files = glob(log_snoreport+'/*.out') + glob(log_snoscan+'/*.out') + glob(
        log_infernal+'/*.out')


for file_ in files:
    start = sp.run(['grep', '-B1', 'rule ', file_], 
                    capture_output=True, text=True).stdout.strip()
    start = start.replace('  ', ' ')
    rule_name = start.split('rule ')[-1].replace(':', '')
    tool_name = rule_name.split('_')[-2]
    if tool_name == 'infernal':
        tool_name = 'infernal_rfam'
    end = sp.run(['grep', '-B1', 'Finished job', file_], 
                    capture_output=True, text=True).stdout.strip()
    if end == '':
        if sp.run(['grep', 'CANCELLED', file_], 
            capture_output=True, text=True).stdout.strip() != '':
            print(f'{tool_name} running time maxed out at 3600 min (60h) on {species}')
            tool_time[tool_name] = max_time 
    else:
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
        tool_time[tool_name] = time_diff
            


df_val = []
for tool, total_time in tool_time.items():
    print(f'Total time {tool} on {species} genome: {total_time} min')
    df_val.append([tool, species, total_time, 'Beluga_2_cpu'])

df = pd.DataFrame(df_val, columns=[
                'tool', 'species_genome', 'runtime_min', 'gpu_type'])
df.to_csv(output_runtime, sep='\t', index=False)
