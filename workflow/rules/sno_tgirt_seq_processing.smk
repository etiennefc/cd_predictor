### T. thermophila and S. pombe TGIRT-Seq datasets are 
### kept as external validation datasets 

rule get_chr_size_tgirt:
    """ Get the chr size from genome fasta. Species is defined in 
        sno_families.smk in the inputs of rule get_sno_sequences."""
    input:
        genome = get_species_genome
    output:
        chr_size = 'data/references/chr_size/{species}_chr_size.tsv'
    wildcard_constraints:
        species = join_list(config['species_tgirt'], ["homo_sapiens", "mus_musculus", 
                "saccharomyces_cerevisiae", "drosophila_melanogaster", "danio_rerio", 
                "tetrahymena_thermophila", "plasmodium_falciparum", "gallus_gallus", "macaca_mulatta"])
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input.genome} --fai-idx {wildcards.species}_temp.fai && "
        "cut -f1,2 {wildcards.species}_temp.fai > {output.chr_size} && "
        "rm {wildcards.species}_temp.fai"

rule get_expressed_snoRNAs_location:
    """ From the snoRNAs with TGIRT-Seq (human, mouse and 
        S. cerevisiae), separate the snoRNAs that are 
        expressed (> 1 TPM in at least one average sample) 
        from those that are not expressed (snoRNA 
        pseudogenes in the case of human snoRNAs). The snoRNA 
        type info comes from snoDB for human (legacy version), 
        rnacentral (for mouse from the Abundance determinants
        paper) and from the UMass yeast snoRNA database 
        (Samarksy et al., NAR, 1999). Then, get their 
        location/sequence from their respective gtf"""
    input: 
        tpm_df = 'data/references/tgirt_seq_output/{species}_merged_tpm_w_biotype.tsv',  # TO DO:  Add these datasets on Zenodo
        sno_type_df = lambda wildcards: config['sno_type_df'][wildcards.species],  # # TO DO:  Add to Zenodo
        gtf = get_species_gtf,
        genome = get_species_genome,
        chr_size = rules.get_chr_size_tgirt.output.chr_size
    output: 
        expressed_sno_df = 'data/references/tgirt_seq_output/{species}_expressed_snoRNAs.tsv'
    params:
        human_pseudosno = 'data/references/tgirt_seq_output/homo_sapiens_pseudogene_snoRNAs.tsv',
        mouse_pseudosno = 'data/references/tgirt_seq_output/mus_musculus_pseudogene_snoRNAs.tsv',
        extension = 15
    wildcard_constraints:
        species=join_list(config['species'], ["homo_sapiens", "mus_musculus", 
                "saccharomyces_cerevisiae"])
    conda:
        "../envs/python_new.yaml"
    script: 
        "../scripts/python/get_expressed_snoRNAs.py"

rule get_D_melanogaster_expressed_snoRNAs_location:
    """ From the (Sklias et al., NAR, 2024) paper (table S6), get the expressed C/D snoRNA 
        and snoRNA pseudogene based on their classification of expressed/not-expressed 
        (cannot reprocess their ovary/head/S2R datasets because they are not yet available).
        Get their location/sequence from the Drosophila gtf"""
    input: 
        tpm_df = 'data/references/tgirt_seq_output/{species}_expression_status_snoRNAs.tsv',  # TO DO:  Add these datasets on Zenodo
        gtf = get_species_gtf,
        genome = get_species_genome,
        chr_size = rules.get_chr_size_tgirt.output.chr_size
    output: 
        expressed_sno_df = 'data/references/tgirt_seq_output/{species}_expressed_snoRNAs.tsv',
        droso_pseudosno = 'data/references/tgirt_seq_output/{species}_pseudogene_snoRNAs.tsv'
    params:
        extension = 15
    wildcard_constraints:
        species=join_list(config['species'], ["drosophila_melanogaster"])
    conda:
        "../envs/python_new.yaml"
    script: 
        "../scripts/python/get_D_melanogaster_expressed_snoRNAs_location.py"

rule tgirt_seq_control:
    """ Check the quality of TGIRT-Seq samples (chicken, macaque and human SKOV). 
        Compare the size-selected SKOV vs fragmented to see if there are differences 
        in snoRNA TPM and TPM ranks. Create also a pie chart of TPM % per biotype."""
    input: 
        tpm_tgirt = expand('data/references/tgirt_seq_validation/{sp}_merged_tpm.tsv', sp=['homo_sapiens', 'macaca_mulatta', 'gallus_gallus']),
        gtf = expand('data/references/gtf/{sp}.gtf', sp=['homo_sapiens', 'macaca_mulatta', 'gallus_gallus']),
        human_tpm = 'data/references/tgirt_seq_output/homo_sapiens_merged_tpm_w_biotype.tsv'
    output: 
        tpm_skov = 'results/figures/scatterplot/TPM_SKOV_comparison_frag_size_selected.svg',
        tpm_skov_rank = 'results/figures/scatterplot/TPM_Rank_SKOV_comparison_frag_size_selected.svg',
        pie_biotype = 'results/figures/pie/biotype_proportion_TGIRT_chicken_monkey_skov.svg',
        df_skov = 'data/references/tgirt_seq_validation/homo_sapiens_SKOV_expressed_tpm_w_biotype.tsv',
        df_chicken = 'data/references/tgirt_seq_validation/gallus_gallus_expressed_tpm_w_biotype.tsv',
        df_macaque = 'data/references/tgirt_seq_validation/macaca_mulatta_expressed_tpm_w_biotype.tsv'
    params:
        tpm_tgirt = '_SEP_'.join(expand('data/references/tgirt_seq_validation/{sp}_merged_tpm.tsv', sp=['homo_sapiens', 'macaca_mulatta', 'gallus_gallus'])),
        gtf = '_SEP_'.join(expand('data/references/gtf/{sp}.gtf', sp=['homo_sapiens', 'macaca_mulatta', 'gallus_gallus'])),
    #conda:
    #    "../envs/python_new.yaml"
    script:
        "../scripts/python/figures/tgirt_seq_control.py"
    #shell:
    #    "source pyenv/bin/activate && " 
    #    "python3 scripts/python/figures/tgirt_seq_control.py {params.tpm_tgirt} "
    #    "{params.gtf} {input.human_tpm} {output.tpm_skov} {output.tpm_skov_rank} "
    #    "{output.pie_biotype}"