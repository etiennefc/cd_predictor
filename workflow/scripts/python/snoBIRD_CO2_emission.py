#!/usr/bin/python3
import subprocess as sp
import json  
import pandas as pd
import os 

cookies_beluga = snakemake.input.cookies_beluga
cookies_narval = snakemake.input.cookies_narval
jobs_beluga = snakemake.input.jobs_list_beluga
jobs_narval = snakemake.input.jobs_list_narval
url_beluga = snakemake.params.url_beluga
url_narval = snakemake.params.url_narval

beluga_co2 = 0
# Iterate over all jobs that were run for SnoBIRD on Beluga cluster
with open(jobs_beluga, 'r') as f:
    for line in f:
        if 'fcouture' in line:
            job_id = line.split('fcouture ')[1].split(' ')[0]
            sp.call(
                f'wget --load-cookies {cookies_beluga} {url_beluga}{job_id}/value/cost.json -O beluga_{job_id}', shell=True)
            with open(f'beluga_{job_id}', 'r') as temp_f:
                dictio = json.load(temp_f)
                if 'co2_emissions_kg' in dictio.keys():
                    beluga_co2 += float(dictio['co2_emissions_kg'])




print(f'Total CO2 production (kg) from optimizing and running SnoBIRD on Beluga: {beluga_co2}')



narval_co2 = 0
# Iterate over all jobs that were run for SnoBIRD on Narval cluster
with open(jobs_narval, 'r') as f:
    for line in f:
        if 'fcouture' in line:
            job_id = line.split('fcouture ')[1].split(' ')[0]
            sp.call(
                f'wget --load-cookies {cookies_narval} {url_narval}{job_id}/value/cost.json -O narval_{job_id}', shell=True)
            with open(f'narval_{job_id}', 'r') as temp_f:
                if os.path.getsize(f'narval_{job_id}') > 0:
                    dictio = json.load(temp_f)
                    if 'co2_emissions_kg' in dictio.keys():
                        narval_co2 += float(dictio['co2_emissions_kg'])

print(f'Total CO2 production (kg) from optimizing and running SnoBIRD on Narval: {narval_co2}')


# Create final df
df = pd.DataFrame([[beluga_co2, narval_co2, beluga_co2+narval_co2, 'between_Augst_2023_and_end_of_January_2025']], 
            columns=['Beluga_total_CO2_kg', 'Narval_total_CO2_kg', 'Total_CO2', 'time_lapse'])
print(df)
df.to_csv(snakemake.output.table, sep='\t', index=False)
sp.call('rm narval_* beluga_*', shell=True)


