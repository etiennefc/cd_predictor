rule training_transformer:
    """ Train the Transformer with the best hyperparams. Must be connected to internet to load the pretrained model for the first time"""
    input:
        X_train = 'data/references/positives_and_negatives/added_features/added_features_training_set_fixed_length_{fixed_length}nt.tsv',
        y_train = 'data/references/positives_and_negatives/added_features/added_features_training_target_{fixed_length}nt.tsv'
    output:
        model = 'results/predictions/transformer/{fixed_length}/transformer_trained_fold_{fold_num}.pt',
        fold_loss = 'results/predictions/transformer/{fixed_length}/transformer_trained_fold_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/transformer/{fixed_length}/transformer_trained_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/transformer_training_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_training.sh "
	"{params.pretrained_model} {wildcards.fold_num} "
	"{params.random_state} "
	"{input.X_train} {input.y_train} "
	"{output.model} {output.fold_loss} {output.fold_f1_score} &> {log}"


rule hypertuning_transformer_2_classes:
    """ Get the best hyperparameters of the transformer with 2 classes using Grid Search."""
    input:
        X_tuning = 'data/references/positives_and_negatives/added_features/added_features_tuning_set_fixed_length_{fixed_length}nt.tsv',
        y_tuning = 'data/references/positives_and_negatives/added_features/added_features_tuning_target_{fixed_length}nt.tsv'
    output:
        best_hyperparams = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_best_hyperparams.tsv'
    params:
        hyperparams_space = config['hyperparameter_space_transformer_2_classes'],
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/hypertuning_transformer_2_classes_{fixed_length}nt.log"
    shell:
        "bash scripts/bash/hypertuning_transformer_2_classes.sh "
        "{params.random_state} {params.pretrained_model} "
        "{input.X_tuning} {input.y_tuning} "
        "{output.best_hyperparams} &> {log}"


rule training_transformer_2_classes:
    """ Train the Transformer with the best hyperparams. ONLY 2 classes to predict: sno/pseudosno (1) or other (0). Must be connected to internet to load the pretrained model for the first time. ORGIINAL hyperparams: lr=2e-5, batch_size=100"""
    input:
        X_train = 'data/references/positives_and_negatives/added_features/added_features_training_set_fixed_length_{fixed_length}nt.tsv',
        y_train = 'data/references/positives_and_negatives/added_features/added_features_training_target_{fixed_length}nt.tsv'
    output:
        model = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_trained_fold_{fold_num}.pt',
        fold_loss = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_trained_fold_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_trained_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/transformer_2_classes_training_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_2_classes_training.sh "
        "{params.pretrained_model} {wildcards.fold_num} "
        "{params.random_state} "
        "{input.X_train} {input.y_train} "
        "{output.model} {output.fold_loss} {output.fold_f1_score} &> {log}"


rule training_transformer_w_features:
    """ Train the Transformer with sequence and numerical features. Must be connected to internet to load the pretrained model for the first time"""
    input:
        X_train = 'data/references/positives_and_negatives/added_features/added_features_training_set_fixed_length_{fixed_length}nt.tsv',
        y_train = 'data/references/positives_and_negatives/added_features/added_features_training_target_{fixed_length}nt.tsv'
    output:
        model = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_fold_{fold_num}.pt',
        fold_loss = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_fold_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/transformer_w_features/transformer_training_w_features_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_training_w_features.sh "
        "{params.pretrained_model} {wildcards.fold_num} "
        "{params.random_state} "
        "{input.X_train} {input.y_train} "
        "{output.model} {output.fold_loss} {output.fold_f1_score} &> {log}"


rule training_transformer_w_features_2_classes:
    """ Train the Transformer with sequence and numerical features. Must be connected to internet to load the pretrained model for the first time"""
    input:
        X_train = 'data/references/positives_and_negatives/added_features/added_features_training_set_fixed_length_{fixed_length}nt.tsv',
        y_train = 'data/references/positives_and_negatives/added_features/added_features_training_target_{fixed_length}nt.tsv'
    output:
        model = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_2_classes_fold_{fold_num}.pt',
        fold_loss = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_2_classes_fold_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/transformer/{fixed_length}/transformer_trained_w_features_2_classes_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/transformer_w_features_2_classes/transformer_training_w_features_2_classes_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_training_w_features_2_classes.sh "
        "{params.pretrained_model} {wildcards.fold_num} "
        "{params.random_state} "
        "{input.X_train} {input.y_train} "
        "{output.model} {output.fold_loss} {output.fold_f1_score} &> {log}"


rule test_transformer_2_classes:
    """ Test the performance of the trained transform on the actual test
        set. Don't forget to use the best hyperparams (and specify
        the model architecture before loading the weights/parameters
        learned during training). This transformer predicts only 2 classes (other vs sno (sno|pseudosno))."""
    input:
        X_test = 'data/references/positives_and_negatives/added_features/added_features_test_set_fixed_length_211nt.tsv',
        y_test = 'data/references/positives_and_negatives/added_features/added_features_test_target_211nt.tsv',
        #best_hyperparams = rules.hypertuning_gru_added_features_half_normalized_simplified.output.best_hyperparams,
        model = rules.training_transformer_2_classes.output.model 
    output:
        df_metrics_on_test = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_test_metrics_{fixed_length}nt_fold_{fold_num}.tsv',
        test_predictions = 'results/predictions/transformer/{fixed_length}/transformer_2_classes_test_predictions_{fixed_length}nt_fold_{fold_num}.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        python_script = "scripts/python/test_transformer_2_classes.py"
    log:
        "logs/test_transformer_2_classes_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_test_2_classes.sh "
        "{params.pretrained_model} "
        "{input.X_test} {input.y_test} "
        "{input.model} {output.df_metrics_on_test} {output.test_predictions} "
        "{params.python_script} &> {log}"


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

