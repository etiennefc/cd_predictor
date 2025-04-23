rule window_step_analysis:
    """ Find what is the optimal step size based on the number of consecutive 
        positive windows predicted as C/D snoRNAs based on the correctly 
        predicted snoRNAs (expressed and pseudogenes) by the first snoBIRD 
        model."""
    input:
        snoBIRD = 'data/references/best_data_aug_snoBIRD/transformer_2_classes_LR_schedule_trained_fold_8.pt',
        genome = 'data/references/genome_fa/',
        preds = 'results/predictions/snoBIRD/transformer/{fixed_length}/3e-5_3e-6_32_4_data_aug_1_ratio/transformer_2_classes_LR_schedule_test_predictions_194nt_fold_8.tsv',
        all_cd = 'data/references/positives/cd_rfam_filtered_all_sno_pseudo_fixed_length_{fixed_length}nt.tsv'
    output:
        window_preds = 'results/predictions/transformer/window_step_analysis/{fixed_length}/prediction_50nt_up_downstream_positives.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6",
        step_size = 1,
        sp_name_dict = config['species_short_name']
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/window_step_analysis.py"

rule window_prediction_dist:
    """ Create a lineplot that shows the proportion of predicted windows 
        as positives 50 nt before and after the actaul snoRNA window 
        (for accurately predicted CD snoRNAs in the test set)."""
    input:
        df = rules.window_step_analysis.output.window_preds
    output:
        lineplot = 'results/figures/lineplot/transformer/{fixed_length}nt/window_prediction_dist.svg',
        lineplot_strand = 'results/figures/lineplot/transformer/{fixed_length}nt/window_prediction_dist_strand.svg',
        density = 'results/figures/density/sno_length_per_strand_{fixed_length}nt.svg'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/window_prediction_dist.py"
'''
rule cerevisiae_window_prob_dist:
    """ Create a lineplot that shows the prediction probability for each window located 
        100 nt before and after the actual snoRNA window 
        (for all C/D snoRNAs in S. cerevisiae separately)."""
    input:
        all_cd = 'data/references/positives/cd_rfam_filtered_all_fixed_length_190nt.tsv',
        pred_dir = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/'
    output:
        lineplot = 'results/figures/lineplot/snoBIRD/190nt/{cd_cerevisiae}.svg'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/cerevisiae_window_prob_dist.py"

rule cerevisiae_window_prob_dist_all:
    """ Create a lineplot that shows the prediction probability for each window located 
        100 nt before and after the actual snoRNA window 
        (for all C/D snoRNAs in S. cerevisiae in one graph)."""
    input:
        all_cd = 'data/references/positives/cd_rfam_filtered_all_fixed_length_190nt.tsv',
        pred_dir = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/'
    output:
        lineplot = 'results/figures/lineplot/snoBIRD/190nt/all_cd_S_cerevisiae.svg'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/cerevisiae_window_prob_dist_all.py"

rule cerevisiae_filter_windows:
    """ Filter windows per score, consecutive blocks of nt and merge afterwards. Creates test_test.bed"""
    input:
        all_cd = 'data/references/positives/cd_rfam_filtered_all_fixed_length_190nt.tsv',
        #pred_dir = 'results/predictions/snoBIRD/saccharomyces_cerevisiae_OLD/'
        #pred_dir = 'results/predictions/snoBIRD/drosophila_melanogaster/',
        #pred_dir = 'results/predictions/snoBIRD/saccharomyces_cerevisiae_EQUAL/'
        pred_dir = 'results/predictions/snoBIRD/saccharomyces_cerevisiae_data_aug_1_ratio_best_hyperparams/'
        #fa = 'data/references/genome_fa/drosophila_melanogaster/split/',
        #cd_bed_sklias = 'droso_cd_bed_slias_2R.bed'
    output:
        df = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/filtered_preds_all.tsv',
        density_block_length = 'results/figures/density/snoBIRD/pred_block_length.svg'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/cerevisiae_filter_windows.py"
'''
rule pombe_filter_windows:
    """ Filter windows per score, consecutive blocks of nt and merge afterwards.
        Step_size=5 is the best option for S. pombe."""
    input:
        all_cd = 'data/references/positives/filtered_CD_snoRNAs.bed',
        pred_dir = 'results/predictions/snoBIRD/schizosaccharomyces_pombe_step5/',
        genome = 'data/references/genome_fa/schizosaccharomyces_pombe_genome.fa'
    output:
        filtered_preds = 'results/predictions/snoBIRD/schizosaccharomyces_pombe/filtered_preds_step5.bed',
        center_preds = 'results/predictions/snoBIRD/schizosaccharomyces_pombe/filtered_center_preds_step5.bed',
        density_block_length = 'results/figures/density/snoBIRD/pred_block_length_s_pombe.svg',
        bed_overlap_sno = 'results/predictions/snoBIRD/schizosaccharomyces_pombe/overlap_snoBIRD_annotated_CD.bed'
    params:
        fixed_length = 194
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/pombe_filter_windows.py"

rule igv_batch_screenshot:
    """ From a bed, create a batch script which will be fed to IGV 
        to automate screenshot creation at each locus in the bed file."""
    input:
        snoBIRD_bed = 'test_test.bed',  # created by pombe_filter_windows
    output:
        batch_script = 'results/igv_scripts/snoBIRD_s_pombe.batch',

    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/igv_batch_screenshot.py"

rule cerevisiae_FP_overlap:
    """ Find with what genomic elemnts the "positives" (FP and TP) of the 
        first model overlap. (Mostly exon, intron, intergenic, tRNA?)."""
    input:
        pos_bed = 'test_test.bed',
        gtf = 'data/references/gtf/saccharomyces_cerevisiae.gtf'
    output:
        df = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/FP_location_099996_190_no_merge_overlap.tsv'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/cerevisiae_FP_overlap.py"
'''
rule human_cd_windows:
    """ Based on the expressed C/D in human, create/get all windows 
        (of length 190 nt) 100 nt up/downstream of these snoRNAs."""
    input:
        cd_expressed = '', 
        genome = 'data/references/genome_fa/homo_sapiens_genome.fa',
        chr_size = 'data/references/chr_size/homo_sapiens/homo_sapiens_chr_size.tsv'
    output:
        df = 'results/predictions/snoBIRD/homo_sapiens/sno_cd_all_windows.tsv'
    params:
        fixed_length = 190
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/human_cd_windows.py"

rule cerevisiae_sno_pseudo_windows:
    """ Based on the blocks predicted as positives by the first step of snoBIRD, 
        create windows out of these blocks (to then apply the second step to predict if it is a sno or pseudosno. 
        Start off with score>0.99996; block>190 before merging)"""
    input:
        pos_blocks = 'test_test.bed',  # choose the optimal filters and change the name of test_test.bed in rule cerevisiae_filter_windows
        genome = 'data/references/genome_fa/saccharomyces_cerevisiae_genome.fa',
        chr_size = 'data/references/chr_size/saccharomyces_cerevisiae/saccharomyces_cerevisiae_chr_size.tsv'
    output:
        df = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/sno_pseudo_windows.tsv'
    params:
        fixed_length = 190
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/figures/cerevisiae_sno_pseudo_windows.py"

rule predict_cerevisiae_sno_pseudo:
    """ Based on the blocks predicted as positives by the first step of snoBIRD, 
        apply the second step on the generated windows to predict if it is a sno or pseudosno. 
        Start off with score>0.99996; block>190 before merging"""
    input:
        sno_pseudo_windows = rules.cerevisiae_sno_pseudo_windows.output.df,
        snoBIRD_pseudo = 'data/references/trained_transformer_sno_pseudo_4e-5_4e-6_32_25_epochs/transformer_sno_pseudo_trained_fold_7.pt'
    output:
        df = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/sno_pseudo_preds.tsv'
    params:
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 190,
        step_size = 1
    conda:
        "../envs/python_new3.yaml"
    script:
        "../scripts/python/figures/predict_cerevisiae_sno_pseudo.py"

rule filter_predictions_cerevisiae_sno_pseudo:
    """ Filter the second step of predictions to maximize the recall of expressed C/D and minimize the FP"""
    input:
        all_cd = 'data/references/positives/cd_rfam_filtered_all_fixed_length_190nt.tsv',
        sno_pseudo_windows = rules.cerevisiae_sno_pseudo_windows.output.df,
        sno_pseudo_preds = rules.predict_cerevisiae_sno_pseudo.output.df
    output:
        df = 'results/predictions/snoBIRD/saccharomyces_cerevisiae/sno_pseudo_filtered_preds.tsv'
    conda:
        "../envs/python_new3.yaml"
    script:
        "../scripts/python/figures/filter_predictions_cerevisiae_sno_pseudo.py"

rule genome_windows:
    """ Create windows of 190 the Candida genome."""
    input:
        snoBIRD = 'data/references/trained_transformer_2_classes_4e-5_4e-7_16_30/transformer_2_classes_LR_schedule_trained_fold_9.pt',
        genome = 'data/references/genome_fa/candida_albicans_genome.fa'
    output:
        windows = 'test_windows.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 190,
        step_size = 1
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/genome_windows.py"

rule genome_windows_separate_chrom:
    """ Create windows of 190 the Candida genome."""
    input:
        snoBIRD = 'data/references/trained_transformer_2_classes_4e-5_4e-7_16_30/transformer_2_classes_LR_schedule_trained_fold_9.pt',
        genome = 'data/references/genome_fa/candida_albicans/{chr}.fa'  # rule separate_chrom
    output:
        windows = 'test_windows_{chr}.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 190,
        step_size = 1
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/genome_windows_separate_chrom.py"

rule genome_windows_separate_chrom_onnx:
    """ Create windows of 190 the Candida genome."""
    input:
        snoBIRD = 'data/references/trained_transformer_2_classes_4e-5_4e-7_16_30/transformer_2_classes_LR_schedule_trained_fold_9.pt',
        genome = 'data/references/genome_fa/candida_albicans/{chr}.fa'  # rule separate_chrom
    output:
        windows = 'results/onnx/windows_onnx_{chr}.tsv'
    params:
        random_state = 42,
        pretrained_model = "zhihan1996/DNA_bert_6",
        fixed_length = 190,
        step_size = 1,
        strand = 'both'
    conda:
        "../envs/python_new3.yaml"
    script:
        "../scripts/python/genome_windows_separate_chrom_onnx.py"

rule test_snoscan_species:
    """ Predict with snoscan the presence of 
        C/D in human genome."""
    input:
        target_rDNA = rules.download_rDNA.output.rDNA_fa,
        genome = '../../snoBIRD/{species}_genome.fa'
    output:
        predictions = 'results/predictions/snoscan/fixed_length_{fixed_length}nt/predicted_cd_{species}.txt'
    conda: 
        "../envs/snoscan.yaml"
    script:
        "../scripts/python/test_snoscan_species.py"
'''