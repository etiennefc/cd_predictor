# For the negatives, also include snoRNA pseudogenes (control for rfam_id?)

rule random_shuffle_sno:
    """ Shuffle expressed C/D sequences and use that as 
        negative examples from all species"""
    input:
        sno_sequences = rules.get_sno_sequences.output.df
    output:
        shuffled_sno_df = 'data/references/negatives/shuffle_sno/random_shuffle_all_expressed_cd.tsv'
    params:
        random_state = 42
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/random_shuffle_sno.py"

rule filter_rnacentral_tRNA_snRNA_pre_miRNA:
    """ Filter rnacentral ncRNA beds (keep non-overlapping entries 
        of tRNA, snRNA and pre-miRNAs) and return also sequence as 1 df.
        O. tauri is not on RNAcentral (thereby excluded here), S. pombe is on 
        RNAcentral but not included in our training set (thereby excluded here) and 
        T. thermophila is on RNAcentral (but relative to the previous incomplete/unassembled 
        genome version, thereby not usable/translatable to the latest assembled version of 
        its genome with which we work here, thereby excluded here)."""
    input:
        beds = expand(rules.download_rnacentral_ncRNA.output.bed, 
                species=[sp for sp in config['species'] + config['species_tgirt'] 
                            if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe', 'tetrahymena_thermophila']]),
        genomes = get_all_genomes('data/references/genome_fa/*.fa'),
        gallus_gallus_gff = rules.download_gallus_gallus_gff.output.gff
    output:
        df = 'data/references/rnacentral/filtered_tRNA_snRNA_pre_miRNA.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_rnacentral_tRNA_snRNA_pre_miRNA.py"

rule filter_split_rnacentral_HACA_snoRNA_fa:
    """ Filter the fasta of all eukaryote H/ACA box snoRNAs 
        present in RNAcentral to keep only the relevant species. 
        Split in .fa per species."""
    input:
        fa = rules.download_eukaryote_HACA_snoRNAs.output.HACA_fa
    output:
        fa = expand('data/references/HACA/{species}_HACA.fa', species=[sp for sp in 
                    config['species'] + config['species_tgirt'] if sp not in 
                    ['ostreococcus_tauri', 'schizosaccharomyces_pombe', 'tetrahymena_thermophila']])
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/filter_split_rnacentral_HACA_snoRNA_fa.py"

rule blat_haca_genome:
    """ Get the genomic location of a given H/ACA snoRNA sequence 
        in a given species genome using BLAT (H/ACA sequence collected 
        from RNAcentral)."""
    input:
        blat_fake_dependency = rules.download_blat.output.tmp_file,
        sno_sequences = "data/references/HACA/{species}_HACA.fa",
        genome = get_species_genome
    output:
        sno_location = "data/references/HACA/{species}_HACA_location.psl"
    params:
        blat_path = "data/references/blat/blat"
    shell:
        "{params.blat_path} {input.genome} {input.sno_sequences} "
        "{output.sno_location}"

rule format_blat_haca_output:
    """ Format BLAT output to keep only the match with the highest 
        number of matching nucleotide according to the original 
        snoRNA sequence. Update the snoRNA sequence based on the 
        actual location in the species genome (for merged snoRNAs 
        and potential mismatches mainly)."""
    input:
        blat = rules.blat_haca_genome.output.sno_location,
        genome = get_species_genome,
        chr_size = get_chr_size
    output:
        df = 'data/references/HACA/{species}_HACA_location_formated.tsv'
    params:
        dataset_attribute = config['dataset_attributes']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/format_blat_haca_output.py"

rule select_random_ncRNA:
    """ Select randomly chosen ncRNAs (H/ACA box snoRNAs, tRNAs, 
        snRNAs, pre-miRNA) from various species. ******TO IMPLEMENT**********"""
    input:
        t_sn_pre_miRNA = rules.filter_rnacentral_tRNA_snRNA_pre_miRNA.output.df,
        haca = rules.format_blat_haca_output.output.df
    output:
        random_ncRNA = 'data/references/negatives/random_ncRNA/random_ncRNA.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/select_random_ncRNA.py"

rule get_intergenic_intronic_regions:
    """ Select all intronic and intergenic regions in the genomes 
        of various species that do not overlap with expressed C/D or 
        other chosen random ncRNAs. Make sure that these random regions 
        are the same length distributions of the expressed C/D in each 
        species."""
    input:
        expressed_cd_all_sets = rules.tuning_train_test_split_rfam.output.all_positives,
        ncRNA = rules.filter_rnacentral_tRNA_snRNA_pre_miRNA.output.df,
        haca = expand(rules.format_blat_haca_output.output.df, 
                    species=[sp for sp in config['species']+config['species_tgirt'] 
                        if sp not in ['ostreococcus_tauri', 'schizosaccharomyces_pombe', 
                        'tetrahymena_thermophila']]),
        gtf = get_species_gtf,
        genome_size = glob.glob('data/references/chr_size/*.tsv')
    output:
        intronic_regions_bed = 'data/references/negatives/random_regions/all_intronic_regions_{species}.bed',
        intergenic_regions_bed = 'data/references/negatives/random_regions/all_intergenic_regions_{species}.bed'
    params:
        species_dict = config['species_short_name']
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/get_intergenic_intronic_regions.py"

rule select_random_intergenic_intronic_regions:
    """ Select random intronic and intergenic regions in the genomes 
        of various species. Make sure that these random regions 
        are the same length distributions of the expressed C/D in each 
        species."""
    input:
        expressed_cd_all_sets = rules.tuning_train_test_split_rfam.output.all_positives,
        intronic_regions = rules.get_intergenic_intronic_regions.output.intronic_regions_bed,
        intergenic_regions = rules.get_intergenic_intronic_regions.output.intergenic_regions_bed
    output:
        random_intronic_regions = 'data/references/negatives/random_regions/selected_intronic_regions_{species}.bed',
        random_intergenic_regions = 'data/references/negatives/random_regions/selected_intergenic_regions_{species}.bed'

#rule get_final_negatives:
 #   """ From all negative examples (other ncRNA sequences (H/ACA, 
  #      tRNA, snRNA, pre-miRNA), shuffle of C/D sequences, random 
   #     sequences in introns and intergenic regions, and potentially
    #    snoRNA pseudogene sequences?), select the wanted proportion 
     #   of each of these negatives relative to the number of positive 
     #   examples (expressed C/D)."""
    #input:
     #   random_ncRNA = rules.select_random_ncRNA.output.random_ncRNA,  # this previous rule could be included within this rule instead
      #  shuffle_sno = rules.random_shuffle_sno.output.shuffled_sno_df,
       # random_intronic_regions = rules.select_random_intergenic_intronic_regions.output.random_intronic_regions,
        #random_intergenic_regions = rules.select_random_intergenic_intronic_regions.output.random_intergenic_regions,
        #human_snoRNA_pseudogenes = rules.get_expressed_snoRNAs_location.params.human_pseudosno
