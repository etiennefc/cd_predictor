rule training_transformer_w_features:
    """ Train the Transformer with the sequence and added features. Must be connected to internet to load the pretrained model for the first time"""
    input:
        X_train = 'data/references/positives_and_negatives/added_features/added_features_training_set_fixed_length_{fixed_length}nt.tsv',
        y_train = 'data/references/positives_and_negatives/added_features/added_features_training_target_{fixed_length}nt.tsv'
    output:
        model = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_fold_{fold_num}.pt',
        fold_loss = 'results/predictions/transformer/{fixed_length}/transformer_trained_fold_w_features_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/transformer/{fixed_length}/transformer_trained_fold_w_features_{fold_num}_f1_score_per_epoch.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/training_transformer_w_features.py"

rule learning_curve_avg_f1_score_training_transformer:
    """ Create average learning curve (of avg f1-score across 3 classes) 
        across 10 folds on training set."""
    input:
        f1_score_tsv = glob.glob('/home/etienne/Narval/scratch/cd_predictor/workflow/results/predictions/transformer/211/transformer_trained_fold_*f1_score_per_epoch.tsv')
    output:
        learning_curve = 'results/figures/lineplot/transformer/211nt/transformer_training_f1_score_avg_across_fold.svg'
    params:
        num_epoch = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_transformer.py"

rule learning_curve_avg_f1_score_training_transformer_w_features:
    """ Create average learning curve (of avg f1-score across 3 classes) 
        across 10 folds on training set for transformer trained also with 
        4 numerical features."""
    input:
        f1_score_tsv = glob.glob('/home/etienne/Narval/scratch/cd_predictor/workflow/results/predictions/transformer/211/transformer_trained_w_features_*f1_score_per_epoch.tsv')
    output:
        learning_curve = 'results/figures/lineplot/transformer/211nt/transformer_training_f1_score_avg_across_fold_w_features.svg'
    params:
        num_epoch = 25
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_transformer.py"

rule learning_curve_avg_f1_score_training_transformer_2_classes:
    """ Create average learning curve (of avg f1-score across 2 classes (other vs sno (sno|pseudosno))) 
        across 10 folds on training set for transformer trained w sequence only."""
    input:
        f1_score_tsv = glob.glob('/home/etienne/Narval/scratch/cd_predictor/workflow/results/predictions/transformer/211/transformer_2_classes_25_epochs_4e-5_16/transformer_2_classes_t*f1_score_per_epoch.tsv')
    output:
        learning_curve = 'results/figures/lineplot/transformer/211nt/transformer_2_classes_training_f1_score_avg_across_fold.svg'
    params:
        num_epoch = 25
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_transformer.py"

rule learning_curve_avg_f1_score_training_transformer_w_features_2_classes:
    """ Create average learning curve (of avg f1-score across 2 classes (other vs sno (sno|pseudosno))) 
        across 10 folds on training set for transformer trained w sequence and 4 intrinsic features."""
    input:
        f1_score_tsv = glob.glob('/home/etienne/Narval/scratch/cd_predictor/workflow/results/predictions/transformer/211/transformer_trained_w_features_2_classes_fold_*f1_score_per_epoch.tsv')
    output:
        learning_curve = 'results/figures/lineplot/transformer/211nt/transformer_w_features_2_classes_training_f1_score_avg_across_fold.svg'
    params:
        num_epoch = 200
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_transformer.py"
