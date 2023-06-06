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
        training_metrics_per_fold = 'results/predictions/gru/{fixed_length}nt/gru_training_metrics_{fixed_length}nt_last_epoch_per_fold_w_avg.tsv',
        learning_curves = expand('results/figures/lineplot/gru/{fixed_length}nt/gru_training_f1_score_{fixed_length}nt_fold_{fold_num}.svg',
                                        fold_num=[str(i) for i in range(1,11)], allow_missing=True)                                        
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/training_gru.py"

#rule test_gru:
#    """ Test the performance of the trained GRU on the actual test 
#        set. Don't forget to use the best hyperparams (and specify
#        the model architecture before loading the weights/parameters
#        learned during training)."""