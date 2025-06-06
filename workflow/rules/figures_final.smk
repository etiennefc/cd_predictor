rule donut_positives_negatives_initial:
    """ Create donut charts showing the proportion of species
        for the positive examples and the proportion of 
        negative types and species for the negatives."""
    input:
        sets_first_model = expand(
            'data/references/positives_and_negatives/data_augmentation/{set_name}_set_1_ratio_fixed_length_194nt.tsv', 
            set_name=['tuning', 'training', 'test']),
        sets_second_model = expand(
            'data/references/positives_and_negatives/data_augmentation_equal_ratio/{set_name}_set_equal_ratio_fixed_length_194nt.tsv', 
            set_name=['tuning', 'training', 'test'])
    output:
        pie_species = 'results/figures/pie/expressed_CD_per_species.svg',
        pie_neg_type = 'results/figures/pie/negatives_per_biotype.svg',
        pie_neg_species = 'results/figures/pie/negatives_per_species.svg',
        pie_pseudo_species = 'results/figures/pie/CD_pseudogene_per_species.svg'
    params:
        species_colors = config['colors']['species'],
        biotype_colors = config['colors']['biotypes'],
        species_name = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/donut_positives_negatives_initial.py"

rule bar_cd_expression_species_initial:
    """ Create a stacked bar chart per species with TGIRT-Seq dataset (that we 
        had access to at the start of this study) to show 
        the proportion of expressed CD vs pseudogenes."""
    input:
        tgirt_dir = 'data/references/tgirt_seq_output/'
    output:
        bar = 'results/figures/bar/cd_expression_per_species.svg'
    params:
        colors = config['colors']['target'],
        species_short_name = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/bar_cd_expression_species.py"



rule learning_curve_avg_f1_score_training_transformer_first_model:
    """ Create average learning curve (of avg f1-score across 2 classes (other vs sno (sno|pseudosno))) 
        across 10 folds on training set for transformer trained w sequence only."""
    input:
        fake_dep = expand(rules.download_f1_score.output.f1, stage=['Before', 'LR_schedule'], 
                    fold_num=[1,2,3,4,5,6,7,8,9,10])
    output:
        learning_curve = 'results/figures/lineplot/transformer/194nt/transformer_2_classes_training_f1_score_avg_across_fold.svg'
    params:
        num_epoch = 4,
        f1_before_train = 'results/predictions/transformer/194/3e-5_3e-6_32_4_data_aug_1_ratio/',
        f1_score_tsv = 'results/predictions/transformer/194/3e-5_3e-6_32_4_data_aug_1_ratio/'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/learning_curve_avg_f1_score_training_transformer_first_model.py"

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

rule shap_heatmap_sno_pseudo_final:
    """ Create a heatmap with the SHAP values of all positive examples of the second snoBIRD model (sno vs pseudo)."""
    input:
        shap_df = 'results/shap/snoBIRD/all_cd_shap_values_sno_pseudo.tsv',
        cd_df = 'data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_194nt.tsv'
    output:
        heatmap_clustered = 'results/figures/heatmap/shap_sno_pseudo_clustered_final.png',
        heatmap = 'results/figures/heatmap/shap_sno_pseudo_length_filtered_final.png',
        heatmap_species = 'results/figures/heatmap/shap_sno_pseudo_species_filtered_final.png'
    params:
        species_dict = config['colors']['species'],
        sp_name = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/shap_heatmap_sno_pseudo.py"

rule lineplot_shap_pombe_candidate:
    """ Create a lineplot using the avg SHAP values of SnoBIRD's first model 
        for the given CD_531 prediction in S. pombe."""
    input:
        shap_df = 'results/shap/snoBIRD/all_pombe_preds_shap_values.tsv'  # obtained after running SnoBIRD on S. pombe on HPC (in results/intermediate/predictions/SHAP/shap_values_all_predictions.tsv)
        #shap_df = '../../snoBIRD/workflow/results/intermediate/predictions/first_model/SHAP/all_cd_predicted_sno_limits_pombe.tsv'
    output:
        lineplot = 'results/figures/lineplot/s_pombe_CD_531_shap_lineplot.svg'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/lineplot_shap_pombe_candidate.py"


rule bar_IP_Fib_NOP58:
    """ Create a bar plot showing the fold enrichment for Fibrillarin 
        (vs FLAG control) and NOP58-TAP (vs  untagged control)."""
    input:
        fib = 'results/IP_cedric/Fib_fold_enrichment.tsv',
        nop58 = 'results/IP_cedric/NOP58_TAP_fold_enrichment.tsv'
    output:
        bar_fib = 'results/figures/bar/fibrillarin_IP_enrichment.svg',
        bar_nop58 = 'results/figures/bar/NOP58_IP_enrichment.svg'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/bar_IP_Fib_NOP58.py"

rule clustal_omega_sf3b3:
    """ Multi-align with clustal omega sequences from snoBIRD predicted snoRNAs
        in the SF3B3 gene across droso, macaque, chicken and zebrafish. Return 
        the phylogenetic tree and multi-alignment"""
    input:
        fasta = 'data/references/sf3b3_snoRNAs_across_species.fa',  # manually taken from each prediction file (TO PUT ON ZENODO!)
        phylo_tree_nwk = 'data/references/sf3b3_snoRNAs_across_species.nwk'  # obtained from clustalo web server w default parameters
    output:
        tree_fig = 'results/figures/tree/sf3b3_snoRNAs_multi_species_phylogenetic_tree.svg'
    params:
        num_epoch = 4
    conda:
        "../envs/clustalo_python.yaml"
    script:
        "../scripts/python/figures/clustal_omega_sf3b3.py"

