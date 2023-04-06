rule get_three_sets_initial:
    """ Merge the positives and negative examples into the 
        tuning (10%), training (70%) and test (20%) sets. This is 
        the initial split based on sequence only without any 
        snoRNA pseudogenes in the negative examples. It is a 
        stratified split (same proprotion of pos/neg examples 
        across the 3 sets)."""
    input:
        tune_pos = rules.tuning_train_test_split_rfam.output.tuning,
        train_pos = rules.tuning_train_test_split_rfam.output.training,
        test_pos = rules.tuning_train_test_split_rfam.output.test,
        tune_neg = rules.get_final_negatives.output.
    


# rule same as initial but with snoRNA pseudogenes as negatives

# rule to use not only sequence but other specified features as 
# input (box_score, structure stability, terminal_stem_stability)