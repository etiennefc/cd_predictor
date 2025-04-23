rule SHAP_sno_limits_FP:
    """ Get SHAP values of the first model on predicted C/D (test FP predictions) 
        and define start/end and box positions, as well as feature values (box score, 
        terminal stem, normalized sno mfe). Return also prediction of first and 
        second snoBIRD model."""
    input:
        model = 'data/references/best_data_aug_snoBIRD/transformer_2_classes_LR_schedule_trained_fold_8.pt',
        FP_first_model = 'results/predictions/snoBIRD/transformer/194/3e-5_3e-6_32_4_data_aug_1_ratio/transformer_2_classes_LR_schedule_test_predictions_194nt_fold_8.tsv',
        seq_FP = 'data/references/positives_and_negatives/data_augmentation/test_set_1_ratio_fixed_length_194nt.tsv',
        second_model = 'data/references/best_data_aug_sno_pseudo/transformer_sno_pseudo_trained_fold_9.pt'
    output:
        df = 'results/shap/snoBIRD/FP_first_model_sno_limits_and_preds.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 194,
        python_script = "scripts/python/SHAP_sno_limits.py",
        batch_size = 1024
    conda:
        "../envs/python_new2.yaml"
    shell:
        "python3 {params.python_script} "
        "{input.model} {input.FP_first_model} "
        "{input.seq_FP} {input.second_model} "
        "{params.pretrained_model} {params.fixed_length} "
        "{params.batch_size} {output.df}"

rule SHAP_sno_limits_pombe_w_second_filters:
    """ From the SHAP values of the first model on predicted C/D (S. pombe predictions), 
        define start/end and box positions, as well as feature values (box score, 
        terminal stem, normalized sno mfe). Then, using the second model and 
        these features, predict if predicted snoRNA is expressed or pseudogene."""
    input:
        SHAP_pombe = 'results/shap/snoBIRD/all_pombe_preds_shap_values.tsv',
        pombe_pred_seq = 'results/predictions/snoBIRD/schizosaccharomyces_pombe/filtered_center_preds_step5.bed',
        genome = 'data/references/genome_fa/schizosaccharomyces_pombe_genome.fa',
        second_model_pombe_pred = 'results/predictions/snoBIRD/schizosaccharomyces_pombe/filtered_center_preds_step5_sno_pseudo.tsv'
    output:
        df = 'results/shap/snoBIRD/predicted_cd_pombe_sno_limits.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 194,
        python_script = "scripts/python/SHAP_sno_limits_pombe_w_second_filters.py",
        batch_size = 1024
    conda:
        "../envs/python_new2.yaml"
    shell:
        "python3 {params.python_script} "
        "{input.SHAP_pombe} {input.genome} "
        "{input.pombe_pred_seq} {input.second_model_pombe_pred} "
        "{params.pretrained_model} {params.fixed_length} "
        "{params.batch_size} {output.df}"


#rule run_second_snoBIRD:
#    """ Get SHAP values of the first snoBIRD model, then run the second model of snoBIRD and apply feature-based 
#        filters (based on sno start/end/C/D boxes dfined with SHAP values) on given examples. For instance, try the model on the 
#        FP of the first model and on the predicted C/D in S. pombe."""
#    input:
#        FP_first_model = 'results/shap/snoBIRD/all_cd_shap_values_sno_pseudo.tsv',
#        location_FP = 'data/references/positives_and_negatives/data_augmentation/test_set_1_ratio_fixed_length_194nt.tsv',
#        all_examples = 'data/references/sno_pseudo/all_cd_pseudo_all_features_{fixed_length}nt.tsv'
#    output:
#        preds = 'results/predictions/sno_pseudo/filters/filtered_preds_sno_pseudo_test_{fixed_length}.tsv',
#        metrics = 'results/predictions/sno_pseudo/filters/metrics_filtered_preds_sno_pseudo_test_{fixed_length}.tsv'
#    conda:
#        "../envs/python_new3.yaml"
#    script:
#        "../scripts/python/run_second_snoBIRD.py"