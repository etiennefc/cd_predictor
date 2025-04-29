
rule pred_overlap_bedgraph_multiple_filters:
    """ Find the overlap between the snoBIRD predictions and 
        the bedgraph of expression to see which predictions are expressed."""
    input:
        snoBIRD_preds = 'results/predictions/snoBIRD/final/{species}/',  # zenodo
        known_cd = "data/references/sno_type_df/{species}_snotype.tsv",  # zenodo
        known_cd_tetrahymena = 'data/references/sno_type_df/tetrahymena_thermophila_snotype_good_genome.tsv',  # zenodo
        bedgraph = 'results/bedgraph_TGIRT/{species}/',  # zenodo
        chr_size = 'data/references/chr_size/{species}_chr_size.tsv',
        gtf = 'data/references/gtf/{species}.gtf'
    output:
        coverage = 'results/predictions/snoBIRD/bedgraph_overlap/multiple_filters/{species}_TGIRT_coverage.tsv',
        batch_script_known_cd = 'results/igv_scripts/snoBIRD_{species}_known_cd.batch',
        batch_script_new_cd = 'results/igv_scripts/snoBIRD_{species}_new_cd.batch',
        known_cd_df = 'results/predictions/snoBIRD/bedgraph_overlap/multiple_filters/{species}_overlap_known_cd_w_snoBIRD_preds.tsv',
        new_cd_df = 'results/predictions/snoBIRD/bedgraph_overlap/multiple_filters/{species}_potential_new_cd_snoBIRD_preds.tsv'
    params:
        fixed_length = 194,
        screenshot_dir = "~/lab_32GB/Desktop/Etienne/cd_predictor/workflow/results/figures/screenshots/{species}/"
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/pred_overlap_bedgraph_multiple_filters.py"

rule snoBIRD_candidate_table:
    """ Create a table of all interesting candidates that SnoBIRD predicted 
        across the different species."""
    input:
        snoBIRD = expand('results/predictions/snoBIRD/final/{species}/', 
                        species=['homo_sapiens', 'tetrahymena_thermophila', 
                            'macaca_mulatta', 'schizosaccharomyces_pombe', 
                            'gallus_gallus', 'danio_rerio']),
        new_cd_df = expand('results/predictions/{cd_predictors}/bedgraph_overlap/multiple_filters/{species}_potential_new_cd_{cd_predictors}_preds.tsv', 
                        cd_predictors='snoBIRD', species=['homo_sapiens', 'tetrahymena_thermophila', 
                            'macaca_mulatta', 'schizosaccharomyces_pombe', 
                            'gallus_gallus', 'danio_rerio'])
    output:
        final_table= 'results/predictions/snoBIRD/final_candidates_cross_species.tsv'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/snoBIRD_candidate_table.py"

rule venn_annotated_snoBIRD_predictions:
    """ Find the overlap between snoBIRD predictions and annotated snoRNAs 
        (w/r to their expression). """
    input:
        snoBIRD_preds = expand(rules.pred_overlap_bedgraph_multiple_filters.output.coverage, 
                            species=['homo_sapiens', 'tetrahymena_thermophila', 
                            'macaca_mulatta', 'schizosaccharomyces_pombe', 'gallus_gallus', 
                            'drosophila_melanogaster', 'danio_rerio', 'arabidopsis_thaliana', 
                            'plasmodium_falciparum'])
    output:
        venn = 'results/figures/venn/overlap_expressed_annot_snoBIRD_preds_all_species.svg',
        donut = 'results/figures/donut/overlap_annot_snoBIRD_preds_all_species.svg'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/venn_annotated_snoBIRD_predictions.py"


rule bar_expressed_pseudo_species:
    """ Create a stacked bar plot to show the % of expressed vs 
        snoRNA_pseudogene predictions across the different species for which we 
        used SnoBIRD on."""
    input:
        snoBIRD_preds = expand('results/predictions/snoBIRD/final/{species}/', 
                            species=['homo_sapiens', 'tetrahymena_thermophila', 
                            'macaca_mulatta', 'schizosaccharomyces_pombe', 'gallus_gallus', 
                            'drosophila_melanogaster', 'danio_rerio', 
                            'plasmodium_falciparum']),
    output:
        bar = 'results/figures/bar/expressed_pseudo_snoBIRD_all_species.svg'
    params:
        color_dict = config['colors']['target'],
        species_short_name = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/bar_expressed_pseudo_species.py"

rule filter_other_tools_predictions:
    """ Filter the output of snoscan, snoreport and infernal on S. pombe and 
        H. sapiens to get a tsv in the end. These results were obtained after 
        running the tools on a HPC cluster (Narval)."""
    input:
        snoscan_preds = expand('results/predictions/snoscan/{species}/', # The user should run snoscan on these genomes
                species=['schizosaccharomyces_pombe', 'homo_sapiens']),
        snoreport_preds = expand('results/predictions/snoreport2/{species}/', # The user should run snoreport2 on these genomes
                species=['schizosaccharomyces_pombe', 'homo_sapiens']),
        infernal_preds = expand('results/predictions/infernal_rfam/{species}/', # The user should run infernal on these genomes 
                species=['schizosaccharomyces_pombe', 'homo_sapiens']),
        chr_size = expand('data/references/chr_size/{species}_chr_size.tsv', 
                species=['schizosaccharomyces_pombe', 'homo_sapiens'])
    output:
        snoscan_pombe = 'results/predictions/snoscan/snoscan_schizosaccharomyces_pombe_filtered.tsv',
        snoreport_pombe = 'results/predictions/snoreport2/snoreport2_schizosaccharomyces_pombe_filtered.tsv',
        infernal_pombe = 'results/predictions/infernal_rfam/infernal_rfam_schizosaccharomyces_pombe_filtered.tsv',
        snoscan_human = 'results/predictions/snoscan/snoscan_homo_sapiens_filtered.tsv',
        snoreport_human = 'results/predictions/snoreport2/snoreport2_homo_sapiens_filtered.tsv',
        infernal_human = 'results/predictions/infernal_rfam/infernal_rfam_homo_sapiens_filtered.tsv'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/filter_other_tools_predictions.py"

rule pred_overlap_bedgraph_multiple_filters_other_tools:
    """ Find the overlap between the other tools' predictions and 
        the bedgraph of expression to see which predictions are expressed."""
    input:
        preds = 'results/predictions/{cd_predictors}/{cd_predictors}_{species}_filtered.tsv', 
        known_cd = "data/references/sno_type_df/{species}_snotype.tsv",
        known_cd_tetrahymena = 'data/references/sno_type_df/tetrahymena_thermophila_snotype_good_genome.tsv',
        bedgraph = 'results/bedgraph_TGIRT/{species}/',
        chr_size = 'data/references/chr_size/{species}_chr_size.tsv',
        gtf = 'data/references/gtf/{species}.gtf'
    output:
        coverage = 'results/predictions/{cd_predictors}/bedgraph_overlap/multiple_filters/{species}_TGIRT_coverage.tsv',
        batch_script_known_cd = 'results/igv_scripts/{cd_predictors}_{species}_known_cd.batch',
        batch_script_new_cd = 'results/igv_scripts/{cd_predictors}_{species}_new_cd.batch',
        known_cd_df = 'results/predictions/{cd_predictors}/bedgraph_overlap/multiple_filters/{species}_overlap_known_cd_w_{cd_predictors}_preds.tsv',
        new_cd_df = 'results/predictions/{cd_predictors}/bedgraph_overlap/multiple_filters/{species}_potential_new_cd_{cd_predictors}_preds.tsv'
    params:
        fixed_length = 194,
        screenshot_dir = "~/lab_32GB/Desktop/Etienne/cd_predictor/workflow/results/figures/screenshots/{cd_predictors}/{species}/"
    wildcard_constraints:
        species=join_list(config['species_tgirt'], ["schizosaccharomyces_pombe",    
              "homo_sapiens"]),
        tool = join_list(config['cd_predictors'], ['snoscan', 'snoreport2', 'infernal_rfam'])
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/pred_overlap_bedgraph_multiple_filters_other_tools.py"

rule venn_annotated_other_tools_predictions:
    """ Find the overlap between other tools predictions and annotated snoRNAs 
        (w/r to their expression). """
    input:
        preds = expand(rules.pred_overlap_bedgraph_multiple_filters_other_tools.output.coverage, 
                            species=['homo_sapiens', 'schizosaccharomyces_pombe'], allow_missing=True)
    output:
        venn = 'results/figures/venn/overlap_expressed_annot_{cd_predictors}_preds_all_species.svg',
        donut = 'results/figures/donut/overlap_annot_{cd_predictors}_preds_all_species.svg'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/venn_annotated_other_tools_predictions.py"

rule upset_test_set:
    """ Create an upset plot for the predictions between tools on the test set 
        (one for TP and one for TN). Create also donut plot per confusion value 
        and per tool to highlight potential prediction biases in species and/or biotype"""
    input:
        #snoBIRD = 'results/predictions/snoBIRD/fixed_length_194nt/3e-5_3e-6_32_4_data_aug_1_ratio/transformer_2_classes_LR_schedule_test_predictions_194nt_fold_8.tsv',
        snoBIRD = expand(rules.download_f1_score.output.f1, stage=['LR_schedule'], 
                    fold_num=[8]),
        other_tool = expand('results/predictions/{cd_predictors}/fixed_length_194nt/test_predictions.tsv', 
                    cd_predictors=['infernal_rfam', 'snoscan', 'snoreport2'])
    output:
        TP_FN_upset = 'results/figures/upset/TP_FN_test_set_all_tools.svg',
        TN_FP_upset = 'results/figures/upset/TN_FP_test_set_all_tools.svg',
        donut_conf_value_species = 'results/figures/donut/conf_value_test_set_all_tools_species_hue.svg',
        donut_conf_value_biotype = 'results/figures/donut/conf_value_test_set_all_tools_biotype_hue.svg'
    params:
        biotype_colors = config['colors']['biotypes'],
        species_colors = config['colors']['species'],
        models_colors = config['colors']['predictors']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/upset_test_set.py"

rule snoBIRD_test_set_preds_features:
    """ Create different plots to show the features of expressed CD and snoRNA 
        pseudogenes accurately predicted by SnoBIRD in the test set 
        (box scores, stability, terminal stem stability) (all species in test 
        set combined)."""
    input:
        snoBIRD_test_set = 'results/predictions/snoBIRD/fixed_length_194nt/3e-5_3e-6_32_4_data_aug_1_ratio/test_set_snoBIRD_predictions.tsv',
        test_set = 'data/references/positives_and_negatives/data_augmentation/test_set_1_ratio_fixed_length_194nt.tsv'
    output:
        box_score = 'results/figures/density/box_score_sno_pseudo_test_set.svg',
        c_score = 'results/figures/density/c_score_sno_pseudo_test_set.svg',
        d_score = 'results/figures/density/d_score_sno_pseudo_test_set.svg',
        c_prime_score = 'results/figures/density/c_prime_score_sno_pseudo_test_set.svg',
        d_prime_score = 'results/figures/density/d_prime_score_sno_pseudo_test_set.svg',
        terminal_stability = 'results/figures/density/terminal_stem_stability_sno_pseudo_test_set.svg',
        structure_stability = 'results/figures/density/global_stability_sno_pseudo_test_set.svg',
    params:
        target_colors = config['colors']['target']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/snoBIRD_test_set_preds_features.py"

rule snoBIRD_species_preds_features:  
    """ Create different plots to show the features of expressed CD and snoRNA 
        pseudogenes predicted by SnoBIRD across different species genomes 
        (box scores, stability, terminal stem stability). """
    input:
        snoBIRD = expand('results/predictions/snoBIRD/final/{species}/', 
                        species=['homo_sapiens', 'tetrahymena_thermophila', 
                            'macaca_mulatta', 'schizosaccharomyces_pombe', 'gallus_gallus', 
                            'drosophila_melanogaster', 'danio_rerio', 'plasmodium_falciparum'])
    output:
        boxes_and_stability = 'results/figures/density/boxes_scores_and_stabilities_sno_pseudo_all_predicted_species.svg'
    params:
        target_colors = config['colors']['target']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/snoBIRD_species_preds_features.py"

rule genomic_element_overlap_preds:
    """ Find the overlap between intron, exon, intergenic with the predictions 
        of SnoBIRD and the other tools on S. pombe and H. sapiens. Find it also 
        for SnoBIRD with other species for which we predicted on."""
    input:
        snoBIRD_preds = expand('results/predictions/snoBIRD/final/{species}/', 
                            species=['homo_sapiens', 'tetrahymena_thermophila', 
                            'macaca_mulatta', 'schizosaccharomyces_pombe', 'gallus_gallus', 
                            'drosophila_melanogaster', 'danio_rerio', 'arabidopsis_thaliana', 
                            'plasmodium_falciparum']),
        gtf = expand('data/references/gtf/{species}.gtf', 
                            species=['homo_sapiens', 'tetrahymena_thermophila', 
                            'macaca_mulatta', 'schizosaccharomyces_pombe', 'gallus_gallus', 
                            'drosophila_melanogaster', 'danio_rerio', 'arabidopsis_thaliana', 
                            'plasmodium_falciparum']),
        chr_size = expand('data/references/chr_size/{species}_chr_size.tsv', 
                            species=['homo_sapiens', 'tetrahymena_thermophila', 
                            'macaca_mulatta', 'schizosaccharomyces_pombe', 'gallus_gallus', 
                            'drosophila_melanogaster', 'danio_rerio', 'arabidopsis_thaliana', 
                            'plasmodium_falciparum']),
        snoscan_preds = expand('results/predictions/snoscan/snoscan_{species}_filtered.tsv', 
                species=['schizosaccharomyces_pombe', 'homo_sapiens']),
        snoreport_preds = expand('results/predictions/snoreport2/snoreport2_{species}_filtered.tsv', 
                species=['schizosaccharomyces_pombe', 'homo_sapiens']),
        infernal_preds = expand('results/predictions/infernal_rfam/infernal_rfam_{species}_filtered.tsv', 
                species=['schizosaccharomyces_pombe', 'homo_sapiens']),
    output:
        bar_human_pombe = 'results/figures/bar/overlap_preds_genomic_element_human_pombe.svg',
        bar_other_species = 'results/figures/bar/overlap_preds_genomic_element_other_species.svg'
    params:
        color_dict = config['colors']['genomic_element']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/genomic_element_overlap_preds.py"



rule snoBIRD_runtime:
    """ Based on log files, get the total time elapsed for the snoBIRD pipeline
        for different species."""
    input:
        logs = "logs/snoBIRD/{species}/"  # logs obtained after running snoBIRD on HPC cluster
    output:
        runtime = "results/runtime/snoBIRD/{species}_runtime.tsv"
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/snoBIRD_runtime.py"

rule other_tools_runtime:
    """ Based on log files, get the total time elapsed for the other tool 
        prediction (snoreport, snoscan, infernal_rfam) for different species."""
    input:
        log_snoreport = "logs/snoreport2/{species}/",  # log obtained after running snoreport2 on HPC cluster,
        log_snoscan = "logs/snoscan/{species}/",  # log obtained after running snoreport2 on HPC cluster,
        log_infernal = "logs/infernal_rfam/{species}/",  # log obtained after running snoreport2 on HPC cluster,
    output:
        runtime = "results/runtime/other_tools/{species}_runtime.tsv"
    params:
        max_time = 3600  # minutes
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/other_tools_runtime.py"

rule other_tools_runtime_pombe_human:
    """ Based on log files, get the total time elapsed for the other tool 
        prediction (snoreport, snoscan, infernal_rfam) for S. pombe and human (on Narval)."""
    input:
        log_snoreport_pombe = "logs/snoreport2/schizosaccharomyces_pombe/",  # log obtained after running snoreport2 on HPC cluster,
        log_snoscan_pombe = "logs/snoscan/schizosaccharomyces_pombe/",  # log obtained after running snoreport2 on HPC cluster,
        log_infernal_pombe = "logs/infernal_rfam/schizosaccharomyces_pombe/",  # log obtained after running snoreport2 on HPC cluster,
        log_snoreport_human = "logs/snoreport2/homo_sapiens/",  # log obtained after running snoreport2 on HPC cluster,
        log_snoscan_human = "logs/snoscan/homo_sapiens/",  # log obtained after running snoreport2 on HPC cluster,
        log_infernal_human = "logs/infernal_rfam/homo_sapiens/",  # log obtained after running snoreport2 on HPC cluster,
    output:
        runtime = 'results/runtime/other_tools/schizosaccharomyces_pombe_homo_sapiens.tsv'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/other_tools_runtime_pombe_human.py"

rule runtime_scatter:
    """ Create a scatter plot of runtime of snoBIRD and other tools 
        prediction (snoreport, snoscan, infernal_rfam) for different species on V100 GPUs."""
    input:
        snoBIRD = expand(rules.snoBIRD_runtime.output.runtime, 
                species=['tetrahymena_thermophila', 
                        'drosophila_melanogaster', 'macaca_mulatta', 'gallus_gallus', 
                        'homo_sapiens_chr1']),
        other_tool = expand(rules.other_tools_runtime.output.runtime,
                species=['tetrahymena_thermophila', 
                        'drosophila_melanogaster', 'macaca_mulatta', 'gallus_gallus', 
                        'homo_sapiens_chr1']),
        genome_dir = "data/references/genome_fa/"
    output:
        scatter = "results/figures/scatter/runtime_tools_species.svg"
    params:
        colors = config['colors']['predictors']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/runtime_scatter.py"

rule runtime_bar:
    """ Create a bar plot of runtime of snoBIRD and other tools 
        prediction (snoreport, snoscan, infernal_rfam) for S. pombe and human on A100 GPU."""
    input:
        snoBIRD = expand(rules.snoBIRD_runtime.output.runtime, 
                species=['homo_sapiens', 'schizosaccharomyces_pombe']),
        other_tool = rules.other_tools_runtime_pombe_human.output.runtime,
        genome_dir = "data/references/genome_fa/"
    output:
        bar_pombe = "results/figures/bar/runtime_tools_pombe_A100.svg",
        bar_human = "results/figures/bar/runtime_tools_human_A100.svg"
    params:
        colors = config['colors']['predictors']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/runtime_bar.py"

rule snoBIRD_runtime_cumul_lineplot:
    """ Based on log files, get the time elapsed per rule (step) in the snoBIRD 
        pipeline for different species. Create a lineplot showing the 
        cumulative time per step depending on the genome size."""
    input:
        logs = expand("logs/snoBIRD/{species}/", species=[
                'tetrahymena_thermophila', 'drosophila_melanogaster', 
                'macaca_mulatta', 'gallus_gallus', 'homo_sapiens_chr1'])  # logs obtained after running snoBIRD on HPC cluster
    output:
        lineplot = "results/figures/lineplot/snoBIRD_cumul_runtime_per_species.svg"
    params:
        color_dict = config['colors']['species']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/snoBIRD_runtime_cumul_lineplot.py"

rule step_size_analysis:
    """ Find which step_size is optimal to limit FP and FN in S. pombe and 
        human."""
    input:
        preds_s1 = "results/predictions/snoBIRD/step_size_s1/{species}/",  # all positive windows obtained after running snoBIRD with step_size=1
        known_cd = "data/references/sno_type_df/{species}_snotype.tsv",  # from get_s_pombe_sno_type or concat sno and pseudo CD from whole dataset
        input_fasta_dir = "data/references/genome_fa/{species}_chunks/" # from split_chr on HPC
    output:
        step_size_plot = "results/figures/lineplot/step_size_analysis_{species}.svg",
        df = "results/step_size/preds_per_step_size_{species}.tsv"
    params:
        fixed_length = 194,
        prob_threshold = 0.999,
        min_consecutive_windows = 10
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/step_size_analysis.py"

rule snoBIRD_CO2_emission:
    """ Based on DRAC job summary pages, compute the total CO2 emission of the project."""
    input:
        cookies_beluga = 'data/cookies_beluga.txt', # for DRAC authentification obtained from Get cookies.txt Chrome extension
        cookies_narval = 'data/cookies_narval.txt', # for DRAC authentification obtained from Get cookies.txt Chrome extension
        jobs_list_narval = 'data/all_jobs_narval.tsv',  # obtained with sacct --user=fcouture --starttime=2023-08-01 --endtime=2025-01-31 | grep -vE 'ba\+|ex\+'
        jobs_list_beluga = 'data/all_jobs_beluga.tsv'  # obtained with sacct --user=fcouture --starttime=2023-08-01 --endtime=2025-01-31 | grep -vE 'ba\+|ex\+'
    output:
        table = "results/CO2_emission.tsv"
    params:
        url_narval = 'https://portail.narval.calculquebec.ca/secure/jobstats/fcouture/',
        url_beluga = 'https://portail.beluga.calculquebec.ca/secure/jobstats/fcouture/'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/snoBIRD_CO2_emission.py"

#rule pred_overlap_bedgraph:
#    """ Find the overlap between the snoBIRD (and other tools) predictions and 
#        the bedgraph of expression (one representative per species or max/avg?) 
#        to see which predictions are expressed."""
#    input:
#        snoBIRD_preds = 'results/predictions/snoBIRD/final/{species}/',  # TO PUT ON ZENODO
#        #snocan_preds = '',
#        #snoreport_preds = '',
#        #infernal_preds = '',
#        known_cd = "data/references/sno_type_df/{species}_snotype.tsv",
#        known_cd_tetrahymena = 'data/references/sno_type_df/tetrahymena_thermophila_snotype_good_genome.tsv',
#        bedgraph = 'results/bedgraph_TGIRT/{species}/'
#    output:
#        coverage = 'results/predictions/snoBIRD/bedgraph_overlap/{species}_TGIRT_coverage.tsv'
#    params:
#        fixed_length = 194
#    conda:
#        "../envs/python_new2.yaml"
#    script:
#        "../scripts/python/pred_overlap_bedgraph.py"  # filter only for avg coverage >= 5 reads
#
#rule genome_stats:
#    """ Compute different stats on the genome of species with TGIRT-Seq to explain 
#        the number of SnoBIRD predictions. These stats include genome length, 
#        %GC, intron length, intron % of the genome, intron number"""
#    input:
#        snoBIRD_preds = expand('results/predictions/snoBIRD/final/{species}/', 
#                            species=['homo_sapiens', 'tetrahymena_thermophila', 
#                            'macaca_mulatta', 'schizosaccharomyces_pombe', 'gallus_gallus', 
#                            'drosophila_melanogaster', 'danio_rerio', 'arabidopsis_thaliana', 
#                            'plasmodium_falciparum']),
#        genome = expand('data/references/genome_fa/{species}_genome.fa', 
#                            species=['homo_sapiens', 'tetrahymena_thermophila', 
#                            'macaca_mulatta', 'schizosaccharomyces_pombe', 'gallus_gallus', 
#                            'drosophila_melanogaster', 'danio_rerio', 'arabidopsis_thaliana', 
#                            'plasmodium_falciparum']),
#        gtf = expand('data/references/gtf/{species}.gtf', 
#                            species=['homo_sapiens', 'tetrahymena_thermophila', 
#                            'macaca_mulatta', 'schizosaccharomyces_pombe', 'gallus_gallus', 
#                            'drosophila_melanogaster', 'danio_rerio', 'arabidopsis_thaliana', 
#                            'plasmodium_falciparum']),
#        chr_size = expand('data/references/chr_size/{species}_chr_size.tsv', 
#                            species=['homo_sapiens', 'tetrahymena_thermophila', 
#                            'macaca_mulatta', 'schizosaccharomyces_pombe', 'gallus_gallus', 
#                            'drosophila_melanogaster', 'danio_rerio', 'arabidopsis_thaliana', 
#                            'plasmodium_falciparum'])
#    output:
#        multi_plot = 'results/figures/genome_stats_snoBIRD_preds.svg'
#    params:
#        color_dict = config['colors']['species']
#    conda:
#        "../envs/python_new2.yaml"
#    script:
#        "../scripts/python/figures/genome_stats.py"

#rule create_windows_size_reduction_analysis:
#    """ From the positive examples in the test set, create windows with 
#        increasing number (from the edge to the center) of Ns instead of the 
#        actual sequence to see the window size reduction effect on SnoBIRD 
#        predictions."""
#    input:
#        test_set = 'data/references/positives_and_negatives/data_augmentation/test_set_1_ratio_fixed_length_194nt.tsv',
#        snoBIRD_test_set = 'results/predictions/snoBIRD/fixed_length_194nt/3e-5_3e-6_32_4_data_aug_1_ratio/test_set_snoBIRD_predictions.tsv'
#    output:
#        windows = 'data/references/window_size_reduction_analysis/test_set_reduced_windows_with_Ns.fa'
#    params:
#        fixed_length = 194
#    conda:
#        "../envs/python_new2.yaml"
#    script:
#        "../scripts/python/create_windows_size_reduction_analysis.py"
#
#rule window_size_reduction_analysis:
#    """ From the SnoBIRD predictions on all the test set reduced windows with increasing 
#        Ns (SnoBIRD was run locally), create a line plot showing the % of 
#        windows predicted as snoRNAs when we replace the input window sequence 
#        by increasing number of surrounding Ns. Also overlap on the plot a cumulative line 
#        showing the total of snoRNAs that are of length 194 nt up to 0nt."""
#    input:
#        positives_info = 'data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_194nt.tsv',
#        snoBIRD_N_preds = 'results/predictions/snoBIRD/window_size_reduction_analysis/filtered_positive_windows_test_set_reduced_windows_with_Ns.bed'
#    output:
#        lineplot = 'results/figures/lineplot/window_size_reduction_analysis_test_set_with_Ns.svg'
#    params:
#        fixed_length = 194
#    conda:
#        "../envs/python_new2.yaml"
#    script:
#        "../scripts/python/figures/window_size_reduction_analysis.py"


#rule upset_human_pombe:
#    """ Create an upset plot for the predictions between tools for Human and S. pombe
#        (also showing stacked bars of expressed predictions or not). 
#        Show only merged overlapping predictions per tool. """
#    input:
#        snoBIRD_preds = expand('results/predictions/snoBIRD/final/{species}/', 
#                            species=['homo_sapiens', 'schizosaccharomyces_pombe']),
#        other_preds = rules.filter_other_tools_predictions.output
#    output:
#        human_upset = 'results/figures/upset/human_preds_overlap_all_tools.svg',
#        pombe_upset = 'results/figures/upset/pombe_preds_overlap_all_tools.svg'
#    conda:
#        "../envs/python_new2.yaml"
#    script:
#        "../scripts/python/figures/upset_human_pombe.py"
