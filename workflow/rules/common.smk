""" Common functions used across the workflow"""

def get_species_genome(wildcards):
    # Get the fasta of the genome of a given species
    # The wildcard here is sno_fasta
    species_name = str(wildcards).split('_cd_')[0]
    protists = ["leishmania_major", "dictyostelium_discoideum", "giardia_lamblia"]
    if species_name == 'tetrahymena_thermophila':
        path = rules.download_tetrahymena_genome.output.genome
    elif species_name == 'ostreococcus_tauri':
        path = rules.download_o_tauri_genome.output.genome
    elif species_name in protists:
        path = expand(rules.download_other_protist_genome.output.genome, species=species_name)
    elif species_name in ['homo_sapiens', 'mus_musculus']:
        path = rules.download_mammal_genome.output.genome
    elif species_name in ['saccharomyces_cerevisiae', 'schizosaccharomyces_pombe']:
        path = rules.download_yeast_genome.output.genome
    else:
        path = expand(rules.download_genome.output.genome, species=species_name)
    return path
    print(path)

def get_species_gtf(species):
    # Get the gtf of the genome of a given species
    species = str(species)
    if 'saccharomyces' in species:
        path = rules.download_yeast_gtf.output.gtf
    elif species == 'homo_sapiens':
        path = rules.download_human_gtf.output.gtf
    elif species == 'mus_musculus':
        path = rules.download_mouse_gtf.output.gtf
    return path

def join_list(l, subl, remove=False):
    # From a list l, return a string of all items in subl joined by '|'
    small_list = [a for a in l if a in subl]
    # If we want to remove (instead of only keeping) items of subl from l
    if remove==True:
        small_list = [a for a in l if a not in subl]
    return "{}".format("|".join(small_list))
