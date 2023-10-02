rule scale_data_sno_pseudo_sno_simple_models:
    """ Fit scaler (standardization) on training data only and apply scaler to tuning, training and test set separately. This ensures no data leakage between the sets since we fit the scaler only on the training set."""
    input:
        X_tuning = 'data/references/positives_and_negatives/added_features/added_features_tuning_set_fixed_length_{fixed_length}nt.tsv',
        X_train = 'data/references/positives_and_negatives/added_features/added_features_training_set_fixed_length_{fixed_length}nt.tsv',
        X_test = 'data/references/positives_and_negatives/added_features/added_features_test_set_fixed_length_211nt.tsv'
    output:
        X_tuning_scaled = 'data/references/positives_and_negatives/added_features/added_features_tuning_set_fixed_length_{fixed_length}nt_scaled_sno_pseudo.tsv',
        X_train_scaled = 'data/references/positives_and_negatives/added_features/added_features_training_set_fixed_length_{fixed_length}nt_scaled_sno_pseudo.tsv',
        X_test_scaled = 'data/references/positives_and_negatives/added_features/added_features_test_set_fixed_length_{fixed_length}nt_scaled_sno_pseudo.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/scale_data_sno_pseudo_sno_simple_models.py"



rule hypertuning_sno_pseudo_simple_models: 
    """ Get the best hyperparameters for simple models (logreg, svc, rf, KNN, and gbm on the classification task sno vs pseudosno."""
    input:
        X_tuning_scaled = rules.scale_data_sno_pseudo_sno_simple_models.output.X_tuning_scaled,
        y_tuning = 'data/references/positives_and_negatives/added_features/added_features_tuning_target_{fixed_length}nt.tsv'
    output:
        best_hyperparams = 'results/predictions/sno_pseudo/simple_models/{fixed_length}/{simple_models}_best_hyperparams.tsv'
    params:
        hyperparams_space = lambda wildcards: config['hyperparameter_space'][wildcards.simple_models],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hypertuning_sno_pseudo_simple_models.py"


rule training_sno_pseudo_simple_models:
    """ Train a classifier to predict if it is an expressed C/D or a C/D pseudogene in the given sequence."""
    input:
        X_train_scaled = rules.scale_data_sno_pseudo_sno_simple_models.output.X_train_scaled,
        y_train = 'data/references/positives_and_negatives/added_features/added_features_training_target_{fixed_length}nt.tsv',
        best_hyperparameters = rules.hypertuning_sno_pseudo_simple_models.output.best_hyperparams
    output:
        model = 'results/predictions/sno_pseudo/simple_models/{fixed_length}/{simple_models}_fold_{fold_num}.sav',
        train_metrics = 'results/predictions/sno_pseudo/simple_models/{fixed_length}/{simple_models}_train_metrics_fold_{fold_num}.tsv'
    params:
        random_state = 42,
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_sno_pseudo_simple_models.py"

rule test_sno_pseudo_simple_models:
    """ Test the performance of the sno_pseudo classifier on the test set, but only on examples predicted as C/D by the first classifier (transformer). ***Add specific threshold to filter out FP based on prediction confidence? (ex: if probability is around 0.5 insteadd of a clear 0 or 1 for logreg, consider it as FP?)***"""
    input:
        X_test_scaled = rules.scale_data_sno_pseudo_sno_simple_models.output.X_test_scaled,
        y_test = 'data/references/positives_and_negatives/added_features/added_features_test_target_211nt.tsv',
        transformer_CD_preds = expand('results/predictions/transformer/{fixed_length}/transformer_2_classes_LR_schedule_4e-5_1e-5_16_25_epochs/transformer_2_classes_LR_schedule_test_predictions_{fixed_length}nt_fold_9.tsv', fixed_length=211),
        best_hyperparams = rules.hypertuning_sno_pseudo_simple_models.output.best_hyperparams,
        model = rules.training_sno_pseudo_simple_models.output.model
    output:
        df_metrics_on_test = 'results/predictions/sno_pseudo/simple_models/{fixed_length}/{simple_models}_test_metrics_{fixed_length}_fold_{fold_num}.tsv',
        test_predictions = 'results/predictions/sno_pseudo/simple_models/{fixed_length}/{simple_models}_test_predictions_{fixed_length}_fold_{fold_num}.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/test_sno_pseudo_simple_models.py"
