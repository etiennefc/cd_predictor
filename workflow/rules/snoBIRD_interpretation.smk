rule attention_heatmaps_pombe:
    """ Create heatmaps for each of the 12 attention heads 
        in snoBIRD for given input sequences."""
    input:
        model = 'data/references/best_data_aug_snoBIRD/transformer_2_classes_LR_schedule_trained_fold_8.pt'
    output:
        heatmaps = 'results/figures/heatmap/attention_test.svg'
    params:
        fixed_length = 194,
        pretrained_model = "zhihan1996/DNA_bert_6"
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/attention_heatmaps_pombe.py"
'''
rule shap_snoBIRD:
    """ Compute average SHAP values per nt for all positive examples. TO RUN ON GPU"""
    input:
        model = 'data/references/best_data_aug_snoBIRD/transformer_2_classes_LR_schedule_trained_fold_8.pt',
        #dataset = ''
    output:
        shap_df = 'results/shap/snoBIRD/all_cd_shap_values.tsv'
    params:
        fixed_length = 194,
        pretrained_model = "zhihan1996/DNA_bert_6"
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/shap_snoBIRD.py"

rule shap_snoBIRD_pombe:
    """ Compute average SHAP values per nt for specific examples in S. pombe"""
    input:
        model = 'data/references/best_data_aug_snoBIRD/transformer_2_classes_LR_schedule_trained_fold_8.pt',
        #dataset = ''
    output:
        shap_df = 'results/shap/snoBIRD/pombe_shap_values.tsv'
    params:
        fixed_length = 194,
        pretrained_model = "zhihan1996/DNA_bert_6"
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/shap_snoBIRD.py"
'''
rule shap_heatmap_all_positives:
    """ Create a heatmap with the SHAP values of all positive examples"""
    input:
        #shap_df = rules.shap_snoBIRD.output.shap_df
        shap_df = 'results/shap/snoBIRD/all_cd_shap_values.tsv',
        cd_df = 'data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_194nt.tsv'
    output:
        heatmap_clustered = 'results/figures/heatmap/shap_all_cd_sno_clustered.png',
        heatmap = 'results/figures/heatmap/shap_all_cd_sno_length_filtered.svg',
        heatmap_species = 'results/figures/heatmap/shap_all_cd_sno_species_filtered.png'
    params:
        species_dict = config['colors']['species'],
        sp_name = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/shap_heatmap_all_positives.py"

rule shap_heatmap_all_positives_num_eval:
    """ Create a heatmap with the SHAP values of all positive examples 
        depending on the num_evals parameter in the SHAP computation (higher 
        is usually more accurate as more permutations are done)."""
    input:
        shap_df = 'results/shap/snoBIRD/shap_values_all_cd_{numeval}.tsv',  # obtained from shap_snoBIRD on Narval
        cd_df = 'data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_194nt.tsv'
    output:
        heatmap_clustered = 'results/figures/heatmap/shap_all_cd_sno_clustered_{numeval}.png',
        heatmap = 'results/figures/heatmap/shap_all_cd_sno_length_filtered_{numeval}.png',
        heatmap_species = 'results/figures/heatmap/shap_all_cd_sno_species_filtered_{numeval}.png'
    params:
        species_dict = config['colors']['species'],
        sp_name = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/shap_heatmap_all_positives.py"

rule find_sno_limits_shap:
    """ Based on the SHAP values, find the C and/or the D box to delimit the snoRNA start 
        and end, as well as the C' and D' boxes."""
    input:
        #shap_df = rules.shap_snoBIRD.output.shap_df
        shap_df = 'results/shap/snoBIRD/all_cd_shap_values.tsv',
        cd_df = 'data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_194nt.tsv'
    output:
        df = 'results/shap/snoBIRD/all_cd_predicted_sno_limits.tsv',
        jointplot = 'results/figures/jointplot/len_diff_overlap.svg',
        jointplot_biotype = 'results/figures/jointplot/len_diff_overlap_gene_biotype.svg',
        jointplot_species = 'results/figures/jointplot/len_diff_overlap_gene_species.svg',
        density_box_score = 'results/figures/density/box_score.svg'
    params:
        species_short_name = config['species_short_name'],
        species_colors = config['colors']['species']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/find_sno_limits_shap.py"

rule find_sno_limits_shap_num_eval:
    """ Based on the SHAP values (depending on the max_evals parameters in 
        the explainer used for SHAP computations), find the C and/or the D box 
        to delimit the snoRNA start and end, as well as the C' and D' boxes."""
    input:
        shap_df = 'results/shap/snoBIRD/shap_values_all_cd_{numeval}.tsv',  # obtained from shap_snoBIRD on Narval
        cd_df = 'data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_194nt.tsv'
    output:
        df = 'results/shap/snoBIRD/all_cd_predicted_sno_limits_{numeval}.tsv',
        jointplot = 'results/figures/jointplot/len_diff_overlap_{numeval}.svg',
        jointplot_biotype = 'results/figures/jointplot/len_diff_overlap_gene_biotype_{numeval}.svg',
        jointplot_species = 'results/figures/jointplot/len_diff_overlap_gene_species_{numeval}.svg',
        density_box_score = 'results/figures/density/box_score_{numeval}.svg'
    params:
        species_short_name = config['species_short_name'],
        species_colors = config['colors']['species']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/find_sno_limits_shap.py"

rule box_logo:
    """ Based on the predicted C and D (and prime) boxes, 
        generate logo of the boxes per species."""
    input:
        box_df = rules.find_sno_limits_shap.output.df,
        tuning = 'data/references/positives_and_negatives/data_augmentation/tuning_set_1_ratio_fixed_length_194nt.tsv',
        training = 'data/references/positives_and_negatives/data_augmentation/training_set_1_ratio_fixed_length_194nt.tsv',
        test = 'data/references/positives_and_negatives/data_augmentation/test_set_1_ratio_fixed_length_194nt.tsv'
    output:
        logo_c = 'results/figures/logo/c_logo_{species}.svg',
        logo_d = 'results/figures/logo/d_logo_{species}.svg',
        logo_c_prime = 'results/figures/logo/c_prime_logo_{species}.svg',
        logo_d_prime = 'results/figures/logo/d_prime_logo_{species}.svg'
    params:
        species_short_name = config['species_short_name']
    conda:
        "../envs/logomaker.yaml"
    script:
        "../scripts/python/figures/box_logo.py"

rule find_terminal_stem_stability:
    """ Based on the predicted start/end of snoRNAs (based on the predicted C and D boxes), 
        find the terminal stem stability."""
    input:
        box_df = rules.find_sno_limits_shap.output.df
    output:
        fasta_terminal_stem = 'results/shap/snoBIRD/all_cd_predicted_terminal_stem.fa',
        terminal_stem_df = 'results/shap/snoBIRD/all_cd_predicted_terminal_stem.tsv',
        density_biotype = 'results/figures/density/terminal_stem_stability_length_score_gene_biotype.svg'
    params:
        species_short_name = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/find_terminal_stem_stability.py"

rule find_sno_stability:
    """ Based on the predicted start/end of snoRNAs (based on the predicted C and D boxes), 
        find the snoRNA whole structure stability."""
    input:
        box_df = rules.find_sno_limits_shap.output.df
    output:
        fasta_mfe = 'results/shap/snoBIRD/all_cd_predicted_structure_stability.fa',
        mfe_df = 'results/shap/snoBIRD/all_cd_predicted_structure_stability.tsv',
        density_biotype = 'results/figures/density/structure_stability_gene_biotype.svg'
    params:
        species_short_name = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/find_sno_stability.py"

rule group_all_features_sno_pseudo_normalized:
    """ Regroup the box scores, sno length, terminal stem stability/score 
        and structure stability into one df. Normalize these values (standardize) 
        and return the tuning, training and test sets for the second snoBIRD model 
        (svm, rf, logreg or gbm) which could use these features."""
    input:
        tuning = 'data/references/positives_and_negatives/data_augmentation/tuning_set_1_ratio_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/positives_and_negatives/data_augmentation/training_set_1_ratio_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/positives_and_negatives/data_augmentation/test_set_1_ratio_fixed_length_{fixed_length}nt.tsv',
        structure_mfe = rules.find_sno_stability.output.mfe_df,
        terminal_stem = rules.find_terminal_stem_stability.output.terminal_stem_df,
        box_score = rules.find_sno_limits_shap.output.df
    output:
        sno_pseudo_tuning = 'data/references/sno_pseudo/sno_pseudo_normalized_features_tuning_set_{fixed_length}nt.tsv',
        sno_pseudo_training = 'data/references/sno_pseudo/sno_pseudo_normalized_features_training_set_{fixed_length}nt.tsv',
        sno_pseudo_test = 'data/references/sno_pseudo/sno_pseudo_normalized_features_test_set_{fixed_length}nt.tsv',
        target_tuning = 'data/references/sno_pseudo/sno_pseudo_normalized_features_tuning_target_{fixed_length}nt.tsv',
        target_training = 'data/references/sno_pseudo/sno_pseudo_normalized_features_training_target_{fixed_length}nt.tsv',
        target_test = 'data/references/sno_pseudo/sno_pseudo_normalized_features_test_target_{fixed_length}nt.tsv',
        all_examples = 'data/references/sno_pseudo/all_cd_pseudo_all_features_{fixed_length}nt.tsv',
        all_examples_norm = 'data/references/sno_pseudo/all_cd_pseudo_all_features_normalized_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/group_all_features_sno_pseudo_normalized.py"

rule find_terminal_stem_stability_numeval:
    """ Based on the predicted start/end of snoRNAs (based on the predicted C and D boxes), 
        find the terminal stem stability."""
    input:
        box_df = rules.find_sno_limits_shap_num_eval.output.df
    output:
        fasta_terminal_stem = 'results/shap/snoBIRD/all_cd_predicted_terminal_stem_{numeval}.fa',
        terminal_stem_df = 'results/shap/snoBIRD/all_cd_predicted_terminal_stem_{numeval}.tsv',
        density_biotype = 'results/figures/density/terminal_stem_stability_length_score_gene_biotype_{numeval}.svg'
    params:
        species_short_name = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/find_terminal_stem_stability.py"

rule find_sno_stability_numeval:
    """ Based on the predicted start/end of snoRNAs (based on the predicted C and D boxes), 
        find the snoRNA whole structure stability."""
    input:
        box_df = rules.find_sno_limits_shap_num_eval.output.df
    output:
        fasta_mfe = 'results/shap/snoBIRD/all_cd_predicted_structure_stability_{numeval}.fa',
        mfe_df = 'results/shap/snoBIRD/all_cd_predicted_structure_stability_{numeval}.tsv',
        density_biotype = 'results/figures/density/structure_stability_gene_biotype_{numeval}.svg'
    params:
        species_short_name = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/find_sno_stability.py"

rule group_all_features_sno_pseudo_normalized_numeval:
    """ Regroup the box scores, sno length, terminal stem stability/score 
        and structure stability into one df. Normalize these values (standardize) 
        and return the tuning, training and test sets for the second snoBIRD model 
        (svm, rf, logreg or gbm) which could use these features."""
    input:
        tuning = 'data/references/positives_and_negatives/data_augmentation/tuning_set_1_ratio_fixed_length_{fixed_length}nt.tsv',
        training = 'data/references/positives_and_negatives/data_augmentation/training_set_1_ratio_fixed_length_{fixed_length}nt.tsv',
        test = 'data/references/positives_and_negatives/data_augmentation/test_set_1_ratio_fixed_length_{fixed_length}nt.tsv',
        structure_mfe = rules.find_sno_stability_numeval.output.mfe_df,
        terminal_stem = rules.find_terminal_stem_stability_numeval.output.terminal_stem_df,
        box_score = rules.find_sno_limits_shap_num_eval.output.df
    output:
        sno_pseudo_tuning = 'data/references/sno_pseudo/sno_pseudo_normalized_features_{numeval}_tuning_set_{fixed_length}nt.tsv',
        sno_pseudo_training = 'data/references/sno_pseudo/sno_pseudo_normalized_features_{numeval}_training_set_{fixed_length}nt.tsv',
        sno_pseudo_test = 'data/references/sno_pseudo/sno_pseudo_normalized_features_{numeval}_test_set_{fixed_length}nt.tsv',
        target_tuning = 'data/references/sno_pseudo/sno_pseudo_normalized_features_{numeval}_tuning_target_{fixed_length}nt.tsv',
        target_training = 'data/references/sno_pseudo/sno_pseudo_normalized_features_{numeval}_training_target_{fixed_length}nt.tsv',
        target_test = 'data/references/sno_pseudo/sno_pseudo_normalized_features_{numeval}_test_target_{fixed_length}nt.tsv',
        all_examples = 'data/references/sno_pseudo/all_cd_pseudo_{numeval}_all_features_{fixed_length}nt.tsv',
        all_examples_norm = 'data/references/sno_pseudo/all_cd_pseudo_{numeval}_all_features_normalized_{fixed_length}nt.tsv'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/group_all_features_sno_pseudo_normalized.py"

rule pairplot_features_scatterplot_TPM:
    """ Based on the collected intrinsic features, do a pairplot for each species 
        with a hue of gene_biotype to show potential correlation between features. 
        Also do a scatterplot to show the relationship between the abundance in TPM
        and the different features of mouse vs human snoRNAs."""
    input:
        feature_df = rules.group_all_features_sno_pseudo_normalized.output.all_examples,
        human_sno_tpm_df = 'data/references/tgirt_seq_output/homo_sapiens_expressed_snoRNAs.tsv',
        human_pseudosno_tpm_df = 'data/references/tgirt_seq_output/homo_sapiens_pseudogene_snoRNAs.tsv',
        mouse_sno_tpm_df = 'data/references/tgirt_seq_output/mus_musculus_expressed_snoRNAs.tsv',
        mouse_pseudosno_tpm_df = 'data/references/tgirt_seq_output/mus_musculus_pseudogene_snoRNAs.tsv'
    output:
        pairplot_human = 'results/figures/pairplot/human_features_sno_pseudo_{fixed_length}.svg',
        pairplot_mouse = 'results/figures/pairplot/mouse_features_sno_pseudo_{fixed_length}.svg',
        pairplot_droso = 'results/figures/pairplot/droso_features_sno_pseudo_{fixed_length}.svg',
        scatter_tpm_human = 'results/figures/scatter/human_tpm_features_sno_pseudo_{fixed_length}.svg',
        scatter_tpm_mouse = 'results/figures/scatter/mouse_tpm_features_sno_pseudo_{fixed_length}.svg'
    params:
        color_biotype = config['colors']['target']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/pairplot_features_scatterplot_TPM.py"

rule density_features_TP_FN_FP_TN:
    """ Based on the second model of snoBIRD, create density plots of different intrinsic 
        features based on the type of prediction (True/False Positives/Negatives)."""
    input:
        feature_df = rules.group_all_features_sno_pseudo_normalized.output.all_examples,
        preds = 'results/predictions/sno_pseudo/transformer/{fixed_length}/3e-5_3e-6_16_4_data_aug_equal_ratio_lossfn_1/transformer_sno_pseudo_test_predictions_194nt_fold_9.tsv'
    output:
        density = 'results/figures/density/sno_pseudo_conf_val_terminal_stem_stability_{fixed_length}.svg'
    params:
        color_conf_val = config['colors']['confusion_value']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/density_features_TP_FN_FP_TN.py"

rule lineplot_consecutive_windows_sno_pseudo:
    """ Based on the second model of snoBIRD, a line plot of all test set examples of a certain 
        prediction type (TP, FP, FN, TN based on the pseudogene class) to see how they are 
        predicted (0 or 1) depending on the window in the data augmented examples 
        (_B15, _B14, ..., _A14, _A15)."""
    input:
        preds = 'results/predictions/sno_pseudo/transformer/{fixed_length}/3e-5_3e-6_16_4_data_aug_equal_ratio_lossfn_1/transformer_sno_pseudo_test_predictions_194nt_fold_9.tsv'
    output:
        lineplot_TP = 'results/figures/lineplot/sno_pseudo/{fixed_length}/consecutive_windows_test_set_TP.svg',
        lineplot_FP = 'results/figures/lineplot/sno_pseudo/{fixed_length}/consecutive_windows_test_set_FP.svg',
        lineplot_TN = 'results/figures/lineplot/sno_pseudo/{fixed_length}/consecutive_windows_test_set_TN.svg',
        lineplot_FN = 'results/figures/lineplot/sno_pseudo/{fixed_length}/consecutive_windows_test_set_FN.svg'
    params:
        color_conf_val = config['colors']['confusion_value']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/lineplot_consecutive_windows_sno_pseudo.py"

rule lineplot_shap_sno_pseudo:
    """ Based on the second model of snoBIRD, create a line plot of the shap values (across nt in the window) 
        per sno prediction type (FN, FP, TP, TN)."""
    input:
        shap_df = 'results/shap/snoBIRD/all_cd_shap_values_sno_pseudo.tsv',
        cd_df = 'data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_194nt.tsv'
    output:
        lineplot_TP = directory('results/figures/lineplot/sno_pseudo/shap/TP/'),
        lineplot_FP = directory('results/figures/lineplot/sno_pseudo/shap/FP/'),
        lineplot_TN = directory('results/figures/lineplot/sno_pseudo/shap/TN/'),
        lineplot_FN = directory('results/figures/lineplot/sno_pseudo/shap/FN/')
    params:
        fixed_length = 194
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/lineplot_shap_sno_pseudo.py"

rule shap_heatmap_sno_pseudo:
    """ Create a heatmap with the SHAP values of all positive examples of the second snoBIRD model (sno vs pseudo)."""
    input:
        #shap_df = rules.shap_snoBIRD.output.shap_df
        shap_df = 'results/shap/snoBIRD/all_cd_shap_values_sno_pseudo.tsv',
        cd_df = 'data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_194nt.tsv'
    output:
        heatmap_clustered = 'results/figures/heatmap/shap_sno_pseudo_clustered.png',
        heatmap = 'results/figures/heatmap/shap_sno_pseudo_length_filtered.png',
        heatmap_species = 'results/figures/heatmap/shap_sno_pseudo_species_filtered.png'
    params:
        species_dict = config['colors']['species'],
        sp_name = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/shap_heatmap_sno_pseudo.py"

rule filters_second_snobird:
    """ Try different filter orders to see the get the best performance 
        possible (highest precision on pseudo while keeping a hgh recall 
        on expressed C/D)."""
    input:
        shap_df = 'results/shap/snoBIRD/all_cd_shap_values_sno_pseudo.tsv',
        all_examples = 'data/references/sno_pseudo/all_cd_pseudo_all_features_{fixed_length}nt.tsv'
    output:
        preds = 'results/predictions/sno_pseudo/filters/filtered_preds_sno_pseudo_test_{fixed_length}.tsv',
        metrics = 'results/predictions/sno_pseudo/filters/metrics_filtered_preds_sno_pseudo_test_{fixed_length}.tsv'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/filters_second_snobird.py"

rule filters_second_snobird_numeval:
    """ Try different filter orders to see the get the best performance 
        possible (highest precision on pseudo while keeping a hgh recall 
        on expressed C/D)."""
    input:
        shap_df = 'results/shap/snoBIRD/sno_pseudo_predictions.tsv',  # 2nd model output
        all_examples = 'data/references/sno_pseudo/all_cd_pseudo_{numeval}_all_features_{fixed_length}nt.tsv'
    output:
        preds = 'results/predictions/sno_pseudo/filters/filtered_preds_sno_pseudo_{numeval}_test_{fixed_length}.tsv',
        metrics = 'results/predictions/sno_pseudo/filters/metrics_filtered_preds_sno_pseudo_{numeval}_test_{fixed_length}.tsv'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/filters_second_snobird_numeval.py"