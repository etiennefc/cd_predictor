rule get_three_sets_initial:
    """ Concat the positives and negative examples into the 
        tuning (10%), training (70%) and test (20%) sets. This is 
        the initial split based on sequence only without any 
        snoRNA pseudogenes in the negative examples. It is a 
        stratified split (same proportion of pos/neg examples 
        across the 3 sets, i.e. ~ 20:1 (negatives:positives))."""
    input:
        positives = rules.tuning_train_test_split_rfam.output,
        negatives = rules.get_all_initial_negatives.output
    output:
        tuning = 'data/references/positives_and_negatives/initial/initial_tuning_set.tsv',
        training = 'data/references/positives_and_negatives/initial/initial_training_set.tsv',
        test = 'data/references/positives_and_negatives/initial/initial_test_set.tsv'
    params:
        random_state = 42,
        short_name_dict = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_three_sets_initial.py"

rule get_three_sets_initial_fixed_length:
    """ Concat the positives and negative examples into the 
        tuning (10%), training (70%) and test (20%) sets. This is 
        the initial split based on sequence only without any 
        snoRNA pseudogenes in the negative examples. It is a 
        stratified split (same proportion of pos/neg examples 
        across the 3 sets, i.e. ~ 20:1 (negatives:positives))."""
    input:
        positives = rules.tuning_train_test_split_rfam_fixed_length.output,
        negatives = rules.get_all_initial_negatives_fixed_length.output
    output:
        tuning = 'data/references/positives_and_negatives/initial/initial_tuning_set_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/positives_and_negatives/initial/initial_training_set_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/positives_and_negatives/initial/initial_test_set_fixed_length_{fixed_length}nt.tsv'
    params:
        random_state = 42,
        short_name_dict = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_three_sets_initial_fixed_length.py"


# rule to use not only sequence but other specified features as 
# input (box_score, structure stability, terminal_stem_stability)