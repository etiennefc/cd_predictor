rule get_sno_sequences:
    """ Get all expressed C/D snoRNA sequences (from the literature and 
        TGIRT-seq) regrouped in a fasta file."""
    input:
        sno_literature = rules.merge_sno_location_species.output.df,
        sno_tgirt = expand(rules.get_expressed_snoRNAs_location.output.expressed_sno_df, 
                    species=['homo_sapiens', 'mus_musculus', 'saccharomyces_cerevisiae'])
    output:
        fa = 'data/references/all_expressed_cd_sequences.fa',
        df = 'data/references/all_expressed_cd_sequences_location.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_sno_sequences.py"

rule infernal:
    """ Use Infernal and Rfam covariance 
        models to identify the Rfam family 
        of all expressed C/D box snoRNAs."""
    input:
        sno_fasta = rules.get_sno_sequences.output.fa,
        rfam_cm = rules.download_rfam_covariance_models.output.rfam_cm
    output:
        infernal_tblout = 'data/references/infernal/sno_families.tblout',
        infernal_alignments = 'data/references/infernal/sno_families.txt'
    conda:
        "../envs/infernal.yaml"
    shell:
        "cmpress -F {input.rfam_cm} && "
        "cmscan --cut_ga --rfam --nohmmonly -o {output.infernal_alignments} "
        "--tblout {output.infernal_tblout} {input.rfam_cm} {input.sno_fasta}"

rule filter_infernal:
    """ Filter infernal output to return the Rfam family id per snoRNA."""
    input:
        infernal = rules.infernal.output.infernal_tblout,
        sno_df = rules.get_sno_sequences.output.df
    output:
        df = 'data/references/infernal/sno_families_filtered.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_infernal.py"


rule tuning_train_test_split_rfam:
    """ Split expressed C/D in 3 datasets (tuning (10%), training (70%) 
        and test set (20%)). SnoRNAs of a same Rfam clan (then Rfam family) 
        are all kept within the same set so that we limit overfitting. """
    input:
        sno_rfam = rules.filter_infernal.output.df,
        sno_literature = rules.merge_sno_location_species.output.df,
        sno_tgirt = expand(rules.get_expressed_snoRNAs_location.output.expressed_sno_df, 
                                species=['homo_sapiens', 'mus_musculus', 'saccharomyces_cerevisiae']),
        rfam_clans = rules.download_rfam_clans.output.df
    output:
        tuning = 'data/references/infernal/cd_rfam_filtered_tuning_set.tsv',
        training = 'data/references/infernal/cd_rfam_filtered_training_set.tsv',
        test = 'data/references/infernal/cd_rfam_filtered_test_set.tsv',
        all_positives = 'data/references/positives/cd_rfam_filtered_all.tsv'
    params:
        random_seed = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/tuning_train_test_split_rfam.py"

