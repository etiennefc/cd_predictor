#!/usr/bin/python3
import pandas as pd
import collections as coll 
from sklearn.utils import shuffle, resample
from sklearn.model_selection import train_test_split
seed = snakemake.params.random_seed

# Load dfs
sno_literature = pd.read_csv(snakemake.input.sno_literature, sep='\t')
sno_literature = sno_literature[['gene_id', 'species_name']]
tgirt_dfs = []
for path in snakemake.input.sno_tgirt:
    sno_tgirt = pd.read_csv(path, sep='\t')
    sno_tgirt = sno_tgirt[['gene_id', 'species_name']]
    tgirt_dfs.append(sno_tgirt)
sno_species_df = pd.concat(tgirt_dfs + [sno_literature])

# Simplify species name and merge it to sno_rfam df
sno_species_df['species_name'] = sno_species_df['species_name'].str.replace('mus_', 'M_')
sno_species_df['species_name'] = sno_species_df['species_name'].str.replace('homo_', 'H_')
sno_species_df['species_name'] = sno_species_df['species_name'].str.replace('saccharomyces_', 'S_')
sno_rfam = pd.read_csv(snakemake.input.sno_rfam, sep='\t')
sno_rfam = sno_rfam.merge(sno_species_df, how='left', on='gene_id')

# Remove U3 (RF00012) snoRNA family, as they are way larger and 
# different (other types of boxes) than other types of C/D snoRNAs
sno_rfam = sno_rfam[sno_rfam['rfam_family_id'] != 'RF00012']

# For the snoRNAs with a Rfam id, separate between those with 1 vs 2 to 10 members vs >10 members per family
sno_w_rfam = shuffle(sno_rfam[~sno_rfam['rfam_family_id'].isna()], 
                    random_state=seed).reset_index(drop=True)

sno_families_1_member, sno_big_families, sno_families = {}, {}, {}
for i, group in enumerate(sno_w_rfam.groupby('rfam_family_id')):
    if len(group[1]) == 1:  # 1 member per family
        sno_families_1_member[group[0]] = group[1]
    elif len(group[1]) > 10:  # big families
        sno_big_families[group[0]] = group[1]
    else: # 2 to 10 members per family
        sno_families[group[0]] = group[1]

# For the snoRNAs without a Rfam family or with a rfam id but only 1 member per family, 
# shuffle and split them randomly within the 3 sets (10 %, 70% and 20 % for tuning, 
# training and test set respectively)
sno_no_or_1_rfam = sno_rfam[(sno_rfam['rfam_family_id'].isna()) | 
                            (sno_rfam['rfam_family_id'].isin(sno_families_1_member.keys()))]
sno_no_or_1_rfam = shuffle(sno_no_or_1_rfam, random_state=seed).reset_index(drop=True) 

sno_no_or_1_rfam_rest, sno_no_or_1_rfam_tuning = train_test_split(sno_no_or_1_rfam, train_size=0.9,  
                                                test_size=0.1, random_state=seed)
sno_no_or_1_rfam_train, sno_no_or_1_rfam_test = train_test_split(sno_no_or_1_rfam_rest, train_size=0.78, # 78% of the remaining 90% gives ~ 70%
                                                test_size=0.22, random_state=seed) # 22% of the remaining 90% gives ~ 20%

# For snoRNAs with > 10 members per family, reduce the number to 10 snoRNAs to avoid overrepresentation of some families
# (maximize the number of different species snoRNAs then randomly choose for the remaining snoRNAs)
def largest_factor(num):
    """ Returns the maximal number if times that num is contained within 10."""
    for i in range(10, 0, -1):  # go backward from num-1 to 1
        if i * num <= 10:
            return i * num
            
reduced_big_families = []
for id, df in sno_big_families.items():
    temp_df = df.copy()
    species_dict = dict(coll.Counter(df.species_name)) # nb sno per species (one family per dict)
    species_nb_per_family = len(species_dict.keys())
    largest = largest_factor(species_nb_per_family) # largest number of sno that can be equally chosen from all species
    remaining = 10 - largest  # remaining number of sno that will be randomly chosen within the rest
    i = 0  # i is the number of chosen sno (max(i) = largest)
    j = 0  # j is the row index in a species sno_df 
    sno_ids, sno_rows = [], []
    for n in range(1, largest+1):  # iterate until largest (to ensure that the maximum nb of sno is chosen)
        for name in species_dict.keys():  # iterate across the species present in that rfam sno family
            species_df = df[df['species_name'] == name]
            if (i < largest) & (j <= len(species_df) - 1):  # i < largest and row index must be present within given df 
                random_rows = species_df.sample(n=len(species_df),  # sample randomly all rows in species_df
                                random_state=seed).reset_index(drop=True)
                sno_id = random_rows.loc[j, 'gene_id']  # get sno_id at index j
                if sno_id not in sno_ids:  # if sno not already chosen
                    sno_row = list(random_rows.iloc[j, :].reset_index(drop=True)) # select that snoRNA
                    sno_rows.append(sno_row)
                    sno_ids.append(sno_id)
                    i += 1
        j += 1
    
    # Create list of 1 df from appended sno rows
    sno_rows_df_list = [pd.DataFrame((sno_rows), 
                        columns=['gene_id', 'chr', 'strand', 'start', 'end', 
                                'sequence', 'rfam_family_id', 'species_name'])]
    # Complete with randomly chosen snoRNAs if 10 snoRNAs are not already chosen per family
    if remaining > 0:
        remaining_rows = temp_df[~temp_df['gene_id'].isin(sno_ids)].sample(n=remaining,
                                random_state=seed).reset_index(drop=True)
        sno_rows_df_list.append(remaining_rows)
    reduced_big_families.append(sno_rows_df_list)

# Concat the rows per big family into one df per family, then add them to remaining sno_families dict
big_fam_final_df = [pd.concat(sub) for sub in reduced_big_families]
remaining_sno_dict = sno_families.copy()
for big_fam_df in big_fam_final_df:
    rfam_id = pd.unique(big_fam_df['rfam_family_id'])[0]
    remaining_sno_dict[rfam_id] = big_fam_df


# For snoRNAs with Rfam id and betwen 2 to 10 members per family
remaining_sno_nb = sum([len(v) for k, v in remaining_sno_dict.items()])
selected_sno_nb = remaining_sno_nb + len(sno_no_or_1_rfam_tuning) + len(sno_no_or_1_rfam_train) + len(sno_no_or_1_rfam_test)

# Count the remaining number of sno needed in each set
tuning_nb, train_nb, test_nb = round(selected_sno_nb * 0.1), round(selected_sno_nb * 0.7), round(selected_sno_nb * 0.2)
tuning_nb_remaining = tuning_nb - len(sno_no_or_1_rfam_tuning)  # 52 C/D
train_nb_remaining = train_nb - len(sno_no_or_1_rfam_train)  # 365 C/D
test_nb_remaining = test_nb - len(sno_no_or_1_rfam_test)  # 105 C/D

## Distribute the families pseudo-randomly with respect to the proportion of snoRNAs per set
fam_len_dict = {k:len(v) for k,v in remaining_sno_dict.items()}

tuning, train, test = [], [], []
# To get to the 52 remaining C/D to add in the tuning set,
# randomly choose a family of 10, 10, 10, 7, 4, 3, 2, 2, 2, 2 snoRNAs (52 snoRNAs)
# This combination of 10, 10, 10, 7, 3, 3, ... 2 was manually picked to ensure a total of 52
def sample_cd(dictio, nb_per_family, given_list, big_df, n_samples, rs):
    """ Select from dictio all families of n (nb_per_family) snoRNAs, retrieve all 
        snoRNAs  of that family from big_df and pick randomly n_samples (i.e n different 
        families of size nb_per_family). Ensure that the family was not already selected 
        in given_list from another dataset"""
    selected_fams = [id for id, val in dictio.items() if (len(val) == nb_per_family) & (id not in given_list)]
    df_ = big_df[big_df['rfam_family_id'].isin(selected_fams)].drop_duplicates(subset=['rfam_family_id'])
    sno_selected = resample(df_, n_samples=n_samples, random_state=rs)
    return sno_selected.rfam_family_id.tolist()

tuning_occurences = [3, 1, 1, 1, 4]
tuning_ids = []
for i, number in enumerate([10, 7, 4, 3, 2]):
    ids = sample_cd(remaining_sno_dict, number, tuning_ids, sno_rfam, tuning_occurences[i], seed)
    tuning += ids
    for id in ids:
        tuning_ids.append(id)

filtered_sno = [df.gene_id.tolist() for id, df in remaining_sno_dict.items()]
filtered_sno = [item for sublist in filtered_sno for item in sublist]
filtered_df = sno_rfam[sno_rfam['gene_id'].isin(filtered_sno)]
tuning_df = filtered_df[filtered_df['rfam_family_id'].isin(tuning)]  # 52 C/D

# For test set, randomly choose a family of 10, 10, 10, 8, 7, 7, 6, 6, 6, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2 snoRNAs (105 snoRNAs)
# This combination of 8,8, 7, ... 2 was manually picked to ensure a total of 105
test_occurences = [3, 1, 2, 3, 2, 3, 3, 2]
for i, number in enumerate([10, 8, 7, 6, 5, 4, 3, 2]):
    ids = sample_cd(remaining_sno_dict, number, tuning_ids, sno_rfam, test_occurences[i], seed)
    test += ids
    for id in ids:
        tuning_ids.append(id)

test_df = filtered_df[filtered_df['rfam_family_id'].isin(test)]  # 105 C/D

# For training set, select the remaining snoRNAs not in the test nor tuning sets
train_df = filtered_df[~filtered_df['rfam_family_id'].isin(test+tuning)]  # 365 C/D

# Concat the sets composed of families of 2-10 members to their respective set 
# composed of families with 0 or 1 rfam id
final_tuning = shuffle(pd.concat([sno_no_or_1_rfam_tuning, tuning_df]), 
                        random_state=seed).reset_index(drop=True)
final_train = shuffle(pd.concat([sno_no_or_1_rfam_train, train_df]), 
                        random_state=seed).reset_index(drop=True)
final_test = shuffle(pd.concat([sno_no_or_1_rfam_test, test_df]), 
                        random_state=seed).reset_index(drop=True)       

final_tuning.to_csv(snakemake.output.tuning, sep='\t', index=False)
final_train.to_csv(snakemake.output.training, sep='\t', index=False)
final_test.to_csv(snakemake.output.test, sep='\t', index=False)
