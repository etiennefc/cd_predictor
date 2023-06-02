rule hypertuning_gru:
    """ Define the best hyperparameters of a gated recurrent units (GRU), 
        a type of recurrent neural network. By default, we use the 
        Tree-structured Parzen estimator (TPE), which is a bayesian optimization 
        method that is in general more efficent and effective than grid search.
        We use a 3-fold cross-validation to evaluate the f1-score to find the 
        hyperparams combination which maximizes the f1-score (which accounts for 
        all 3 classes equally)."""
    input:
        X_tuning = rules.onehot_encode_normalize_initial_fixed_length.output.normalized_tuning,
        y_tuning = rules.onehot_encode_normalize_initial_fixed_length.output.target_tuning
    output:
        best_hyperparams = 'results/predictions/gru/{fixed_length}nt/gru_best_params_{fixed_length}nt.tsv'
    params:
        hyperparams_search_space = config['hyperparameter_space_GRU'],
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/hypertuning_gru.py"

rule training_gru:
    """ Train the GRU with the best hyperparams. Use a 10-fold CV 
        approach to evaluate the average performance of the model. 
        Save the model trained on each fold and its performance 
        metrics. Save also the average metrics across the 10 folds."""
    input:
        X_train = rules.onehot_encode_normalize_initial_fixed_length.output.normalized_training,
        y_train = rules.onehot_encode_normalize_initial_fixed_length.output.target_training,
        best_hyperparams = rules.hypertuning_gru.output.best_hyperparams
    output:
        trained_model = expand('results/predictions/gru/{fixed_length}nt/gru_trained_{fixed_length}nt_fold_{fold_num}.pt', 
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        training_metrics = expand('results/predictions/gru/{fixed_length}nt/gru_training_metrics_{fixed_length}nt_fold_{fold_num}.tsv',
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True),
        avg_train_metrics = 'results/predictions/gru/{fixed_length}nt/gru_avg_training_metrics_{fixed_length}nt_across_folds.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_gru.py"