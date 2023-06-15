rule box_score:
    """ Compute a box score for all examples (score of 0, 
        means all the C/D/C'/D' boxes are perfect, and 
        with each mutation wih regards to the consensus, 
        the score increases)."""
    input:
        positives_fa = rules.tuning_train_test_split_rfam_fixed_length.output.all_positives,
        negatives_fa = rules.get_all_initial_negatives_fixed_length.output
    output:
        positives = 'data/references/box_score/all_positives_box_score_{fixed_length}nt.tsv',
        negatives = 'data/references/box_score/all_negatives_and_pseudosno_box_score_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/box_score.py"


rule structure_stability:
    """ Compute the structure stability"""


# terminal_stem stability