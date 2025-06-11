rule snoglobe_inputs:
    """ Create the formated inputs required by snoGloBe for rRNA target 
        prediction in the different validated species. These inputs are a gtf 
        containing the predicted snoRNAs and their target rRNA, as well as a 
        fasta of snoRNA sequences."""
    input:
        fake_input_gtf = 'data/references/gtf/{species}.gtf',
        rDNA_species = 'data/references/rDNA/rDNA_species.fa',
        rDNA_danio = 'data/references/rDNA/rDNA_danio_rerio.fa',
        sno_table = 'results/predictions/snoBIRD/final_candidates_cross_species_filtered.tsv'  # obtained by manual inspection of predictions in final_candidates_cross_species.tsv
    output:
        gtf = 'data/references/snoglobe/{species}_rRNA.gtf',
        sno_fasta = 'data/references/snoglobe/candidate_snoRNAs_{species}.fa',
        chr_fasta_dir = directory('data/references/snoglobe/chr_rRNA_{species}/'),
        target_ids = 'data/references/snoglobe/target_id_{species}.txt'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/snoglobe_inputs.py"


rule snoglobe_pred:
    """ Predict the potential rRNA targets of candidate snoRNAs predicted in 
        different species."""
    input:
        sno_fasta = rules.snoglobe_inputs.output.sno_fasta,
        gtf = rules.snoglobe_inputs.output.gtf,  
        chr_fasta_dir = rules.snoglobe_inputs.output.chr_fasta_dir,
        target_ids = rules.snoglobe_inputs.output.target_ids
    output:
        pred = 'results/snoglobe_preds/predictions_{species}.tsv'
    params:
        snoglobe_path = "git_repos/snoglobe/bin:$PATH"
    conda:
        "../envs/snoglobe.yaml"
    shell:
        "export PATH={params.snoglobe_path} && "
        "snoglobe {input.sno_fasta} {input.target_ids} {input.gtf} "
        "{input.chr_fasta_dir} {output.pred} "
        "--verbose "
        "--seq "
        "--merge"

rule filter_snoglobe_preds:
    """ Filter snoGloBe predictions to keep only those comprising interactions 
        upstream of the D or D' box identified by SnoBIRD. Keep only the 
        prediction with the highest probability per region."""
    input:
        snoglobe = expand(rules.snoglobe_pred.output.pred, 
            species=['danio_rerio', 'homo_sapiens', 'schizosaccharomyces_pombe', 
                'gallus_gallus', 'tetrahymena_thermophila', 'macaca_mulatta']),
        sno_table = 'results/predictions/snoBIRD/final_candidates_cross_species_filtered.tsv'
    output:
        filtered_preds = 'results/predictions/snoBIRD/final_candidates_cross_species_filtered_w_targets.tsv'
    conda:
        "../envs/python_new2.yaml"
    script:
        "../scripts/python/filter_snoglobe_preds.py"