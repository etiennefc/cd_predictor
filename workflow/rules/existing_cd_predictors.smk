rule test_snoreport:
    """ Test the performance of snoreport2 on the test 
        set of fixed C/D length (by first converting the 
        test set to fasta). **WARNING**: Snoreport2 must 
        be installed manually on your computer before 
        using this rule and you might need to change the
        path in the command line below to access snoreport2."""
    input:
        test_set = rules.get_three_sets_initial_fixed_length.output.test
    output:
        predictions = 'results/predictions/snoreport2/fixed_length_{fixed_length}nt/predicted_cd.fa'
    shell:
        """awk 'NR>1 {{print ">"$1" "$2" "$3"\\n"$4}}' {input.test_set} > """
        """temp_snoreport_{wildcards.fixed_length}.fa && """
        """~/snoReport_2/snoreport_2 -i temp_snoreport_{wildcards.fixed_length}.fa """
        """-CD """
        """--positives """
        """-o {output.predictions} && """
        """rm temp_snoreport_{wildcards.fixed_length}.fa"""

rule filter_snoreport_predictions:
    """ Filter snoreport2 predictions to remove duplicates 
        and return a clear target prediction for each example."""
    input:
        predictions_fa = rules.test_snoreport.output.predictions,
        test_set = rules.get_three_sets_initial_fixed_length.output.test
    output:
        predictions_tsv = 'results/predictions/snoreport2/fixed_length_{fixed_length}nt/test_predictions.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_snoreport_predictions.py"

rule test_snoscan:
    """ Predict with snoscan the presence of 
        expressed C/D in the test set."""
    input:
        target_rDNA = rules.download_rDNA.output.rDNA_fa,
        test_set = rules.get_three_sets_initial_fixed_length.output.test
    output:
        predictions = 'results/predictions/snoscan/fixed_length_{fixed_length}nt/predicted_cd.txt'
    conda: 
        "../envs/snoscan.yaml"
    script:
        "../scripts/python/test_sno_scan.py"

rule filter_snoscan_predictions:
    """ Filter snoscan predictions to remove duplicates 
        and return a clear target prediction for each example."""
    input:
        predictions_fa = rules.test_snoscan.output.predictions,
        test_set = rules.get_three_sets_initial_fixed_length.output.test
    output:
        predictions_tsv = 'results/predictions/snoscan/fixed_length_{fixed_length}nt/test_predictions.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_snoscan_predictions.py"

rule test_rfam_infernal:
    """ Use Infernal and Rfam covariance 
        models to predict if testset examples 
        are part of a C/D Rfam family."""
    input:
        test_set = rules.get_three_sets_initial_fixed_length.output.test,
        rfam_cm = rules.download_rfam_covariance_models.output.rfam_cm
    output:
        infernal_tblout = 'results/predictions/infernal_rfam/fixed_length_{fixed_length}nt/predicted_cd.tblout',
        infernal_alignments = 'results/predictions/infernal_rfam/fixed_length_{fixed_length}nt/predicted_cd.txt'
    conda:
        "../envs/infernal.yaml"
    shell:
        """awk 'NR>1 {{print ">"$1" "$2" "$3"\\n"$4}}' {input.test_set} > infernal_temp.fa && """
        """cmpress -F {input.rfam_cm} && """
        """cmscan --cut_ga --rfam --nohmmonly -o {output.infernal_alignments} """
        """--tblout {output.infernal_tblout} {input.rfam_cm} infernal_temp.fa && """
        """rm infernal_temp.fa"""

rule filter_rfam_infernal_predictions:
    """ Filter infernal/rfam predictions to remove duplicates 
        and return a clear target prediction for each example."""
    input:
        predictions_table = rules.test_rfam_infernal.output.infernal_tblout,
        test_set = rules.get_three_sets_initial_fixed_length.output.test
    output:
        predictions_tsv = 'results/predictions/infernal_rfam/fixed_length_{fixed_length}nt/test_predictions.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_rfam_infernal_predictions.py"
