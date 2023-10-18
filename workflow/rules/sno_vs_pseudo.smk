rule hypertuning_transformer_sno_pseudo:
    """ Get the best hyperparameters of the transformer sno_pseudo predictor using Grid Search."""
    input:
        X_tuning = rules.get_three_sets_initial_fixed_length.output.tuning,
        y_tuning = rules.get_three_sets_initial_fixed_length.output.tuning_target
    output:
        best_hyperparams = 'results/predictions/sno_pseudo/transformer/{fixed_length}/transformer_sno_pseudo_best_hyperparams.tsv'
    params:
        hyperparams_space = config['hyperparameter_space_transformer_2_classes'],
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/hypertuning_transformer_sno_pseudo_{fixed_length}nt.log"
    shell:
        "bash scripts/bash/hypertuning_transformer_sno_pseudo.sh "
        "{params.random_state} {params.pretrained_model} "
        "{input.X_tuning} {input.y_tuning} "
        "{output.best_hyperparams} &> {log}"

rule training_transformer_sno_pseudo:
    """ Train the Transformer with the best hyperparams to predict sno (2) vs pseudosno (1). 
        Must be connected to internet to load the pretrained model for the first time. """
    input:
        X_train = rules.get_three_sets_initial_fixed_length.output.training,
        y_train = rules.get_three_sets_initial_fixed_length.output.training_target,
        best_hyperparams = rules.hypertuning_transformer_sno_pseudo.output.best_hyperparams
    output:
        model = 'results/predictions/sno_pseudo/transformer/{fixed_length}/transformer_sno_pseudo_trained_fold_{fold_num}.pt',
        fold_loss = 'results/predictions/sno_pseudo/transformer/{fixed_length}/transformer_sno_pseudo_trained_fold_{fold_num}_loss_per_epoch.tsv',
        fold_f1_score = 'results/predictions/sno_pseudo/transformer/{fixed_length}/transformer_sno_pseudo_trained_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6"
    log:
        "logs/transformer_sno_pseudo_training_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_sno_pseudo_training.sh "
        "{params.pretrained_model} {wildcards.fold_num} "
        "{params.random_state} "
        "{input.X_train} {input.y_train} {input.best_hyperparams} "
        "{output.model} {output.fold_loss} {output.fold_f1_score} &> {log}"

rule test_transformer_sno_pseudo:
    """ Test the performance of the trained transformer on the actual test
        set. Don't forget to use the best hyperparams (and specify
        the model architecture before loading the weights/parameters
        learned during training). This transformer predicts C/D sno (2) vs pseudosno (1)."""
    input:
        X_test = rules.get_three_sets_initial_fixed_length.output.test,
        y_test = rules.get_three_sets_initial_fixed_length.output.test_target,
        best_hyperparams = rules.hypertuning_transformer_sno_pseudo.output.best_hyperparams,
        model = rules.training_transformer_sno_pseudo.output.model
    output:
        df_metrics_on_test = 'results/predictions/sno_pseudo/transformer/{fixed_length}/transformer_sno_pseudo_test_metrics_{fixed_length}nt_fold_{fold_num}.tsv',
        test_predictions = 'results/predictions/sno_pseudo/transformer/{fixed_length}/transformer_sno_pseudo_test_predictions_{fixed_length}nt_fold_{fold_num}.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        python_script = "scripts/python/test_transformer_sno_pseudo.py"
    log:
        "logs/test_transformer_sno_pseudo_{fixed_length}nt_{fold_num}.log"
    shell:
        "bash scripts/bash/transformer_test_2_classes.sh "
        "{params.pretrained_model} "
        "{input.X_test} {input.y_test} {input.best_hyperparams} "
        "{input.model} {output.df_metrics_on_test} {output.test_predictions} "
        "{params.python_script} &> {log}"

rule learning_curve_avg_f1_score_training_transformer_sno_pseudo:
    """ Create average learning curve (of avg f1-score across 2 classes (other vs sno (sno|pseudosno))) 
        across 10 folds on training set for transformer trained w sequence only."""
    input:
        f1_before_train = glob.glob('/home/etienne/Narval/scratch/cd_predictor/workflow/results/predictions/transformer/198/transformer_2_classes_Before_t*f1_score_per_epoch.tsv'),
        f1_score_tsv = glob.glob('/home/etienne/Narval/scratch/cd_predictor/workflow/results/predictions/sno_pseudo/transformer/198/transformer_sno_pseudo_*f1_score_per_epoch.tsv')
    output:
        learning_curve = 'results/figures/lineplot/transformer/198nt/sno_pseudo/transformer_sno_pseudo_training_f1_score_avg_across_fold.svg'
    params:
        num_epoch = 30
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_transformer.py"