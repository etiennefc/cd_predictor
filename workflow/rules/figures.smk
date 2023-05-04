rule length_selection_positives_density:
    """ Create a density plot showing the length 
        distribution of positive examples. Also 
        show different thresholds based on 
        percentile/quartiles to choose a cutoff. 
        The chosen percentile is 95% (< 183 nt)"""
    input:
        positives = rules.get_sno_sequences.output.df
    output:
        density = 'results/figures/density/length_selection_positives.svg'
    params:
        percent_colors = config["colors"]['percent_colors']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/length_selection_positives_density.py"

rule roc_curve_cd_predictors:
    """ Create a ROC curve of the existing CD predictors.*******NOT COMPLETED****"""
    input:
        snoreport = rules.filter_snoreport_predictions.output.predictions_tsv,
        snoscan = rules.filter_snoscan_predictions.output.predictions_tsv,
        infernal_rfam = rules.filter_rfam_infernal_predictions.output.predictions_tsv
    output:
        roc = 'results/figures/roc/existing_cd_predictors_{fixed_length}.svg'
    params:
        colors = config['colors']['predictors']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/roc_curve_cd_predictors.py"
        

#rule pr_curve_cd_predictor:

rule metrics_lineplot_predictors:
    """ Compute the precision, recall, accuracy and f1-score 
        on the test set and also true negative rate (i.e specificity) 
        on the C/D pseudogenes negatives (i.e not expressed human C/D) 
        for the models and display it as a dot plot."""
    input:
        snoreport = rules.filter_snoreport_predictions.output.predictions_tsv,
        snoscan = rules.filter_snoscan_predictions.output.predictions_tsv,
        infernal_rfam = rules.filter_rfam_infernal_predictions.output.predictions_tsv,
        pseudosno_preds = rules.filter_cd_predictors_pseudosno.output.pseudosno_preds
    output:
        dotplot = 'results/figures/lineplot/metrics_existing_cd_predictors_{fixed_length}.svg'
    params:
        colors = config['colors']['predictors']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/metrics_lineplot_predictors.py"

rule FP_FN_initial_analyses_pie:
    """ Create pie charts showing the proportion of species
        or type of negatives predicted as false positives/negatives 
        (FP/FN) by existing CD_predictors. """
    input:
        snoreport = rules.filter_snoreport_predictions.output.predictions_tsv,
        snoscan = rules.filter_snoscan_predictions.output.predictions_tsv,
        infernal_rfam = rules.filter_rfam_infernal_predictions.output.predictions_tsv
    output:
        pie_species = 'results/figures/pie/{error}_per_species_existing_cd_predictors_{fixed_length}.svg',
        pie_neg_type = 'results/figures/pie/{error}_per_negative_type_existing_cd_predictors_{fixed_length}.svg'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/FP_FN_initial_analyses_pie.py"


























