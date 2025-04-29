rule get_chr_size_literature:
    """ Get the chr size from genome fasta."""
    input:
        genome = get_species_genome
    output:
        chr_size = 'data/references/chr_size/{sno_fasta}_chr_size.tsv'
    wildcard_constraints:
        sno_fasta = join_list(config['sno_fasta'], config['sno_fasta'])
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input.genome} --fai-idx {wildcards.sno_fasta}_temp.fai && "
        "cut -f1,2 {wildcards.sno_fasta}_temp.fai > {output.chr_size} && "
        "rm {wildcards.sno_fasta}_temp.fai"

rule blat_sno_genome:
    """ Get the genomic location of a given snoRNA sequence 
        in a given species genome using BLAT (expressed C/D box 
        snoRNAs collected from the literature)."""
    input:
        blat_fake_dependency = rules.download_blat.output.tmp_file,
        sno_sequences = "data/sno_literature_processing/{sno_fasta}.fa",
    output:
        sno_location = "data/sno_literature_processing/{sno_fasta}_location.psl"
    params:
        genome = get_species_genome,
        blat_path = "data/references/blat/blat"
    shell:
        "{params.blat_path} {params.genome} {input.sno_sequences} "
        "{output.sno_location}"

#rule format_rnacentral_output:
#    """ Format RNAcentral output that contains all ncRNA per species to keep 
#        only C/D snoRNAs."""
#    input:
#        ncRNA_rnacentral_bed = rules.download_rnacentral_ncRNA.output.bed,
#        cd_id_rnacentral = "data/references/rnacentral/{species}_cd_id.bed"
#    output:
#        df = 'data/references/{species}_snotype.tsv'
#    wildcard_constraints:
#        species = "{}".format("|".join([]))
#    conda:
#        "../envs/python_new.yaml"
#    script:
#        "../scripts/python/format_rnacentral_output.py"  
#
#rule get_tetrahymena_sno_seq:
#    """ From tetrahymena C/D snoRNA from RNAcentral, get their sequence 
#        (based on Ensembl genome assembly). The sequences will be used to BLAT 
#        against the Tetrahymena genome version with 181 chromosomes"""
#    input:
#        df = expand(rules.format_rnacentral_output.output.df, species='tetrahymena_thermophila'),
#        genome_ensembl = 'data/references/genome_fa/Tetrahymena_thermophila.JCVI-TTA1-2.2.dna.toplevel.fa'
#    output:
#        fa = 'data/references/tetrahymena_rnacentral_cd.fa'
#    conda:
#        "../envs/python_new.yaml"
#    script:
#        "../scripts/python/get_tetrahymena_sno_seq.py"  
#    
#
#rule blat_sno_tetrahymena:
#    """ Get the genomic location (in genome assembly with 181 chr) of given snoRNA sequence s
#        in Tetrahymena using BLAT (snoRNA from RNAcentral)."""
#    input:
#        blat_fake_dependency = rules.download_blat.output.tmp_file,
#        sno_sequences = rules.get_tetrahymena_sno_seq.output.fa,
#        genome = "data/references/genome_fa/tetrahymena_thermophila_genome.fa"
#    output:
#        sno_location = "data/references/tetrahymena_thermophila_snotype_good_genome.tblout"
#    params:
#        blat_path = "data/references/blat/blat"
#    shell:
#        "{params.blat_path} {input.genome} {input.sno_sequences} "
#        "{output.sno_location}"
#
#rule format_blat_output_tetrahymena:
#    """ Format BLAT output to keep only the match with the highest 
#        number of matching nucleotide according to the original 
#        snoRNA sequence."""
#    input:
#        blat = rules.blat_sno_tetrahymena.output.sno_location
#    output:
#        df = 'data/references/tetrahymena_thermophila_snotype_good_genome.tsv'
#    conda:
#        "../envs/python_new.yaml"
#    script:
#        "../scripts/python/format_blat_output_tetrahymena.py"    

rule format_blat_output:
    """ Format BLAT output to keep only the match with the highest 
        number of matching nucleotide according to the original 
        snoRNA sequence. Update the snoRNA sequence based on the 
        actual location in the species genome (for merged snoRNAs 
        and potential mismatches mainly)."""
    input:
        blat = rules.blat_sno_genome.output.sno_location,
        chr_size = rules.get_chr_size_literature.output.chr_size
    output:
        df = 'data/sno_literature_processing/{sno_fasta}_location_formated.tsv'
    params:
        dataset_attribute = lambda wildcards: config['dataset_attributes'][wildcards.sno_fasta],
        extension = 15,
        genome = get_species_genome
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/format_blat_output.py"       

rule merge_sno_location_species:
    """ Merge all snoRNA info (from format_blat_output) from 
        all species in a single table."""
    input:
        dfs = expand(rules.format_blat_output.output.df, **config)
    output:
        df = 'data/references/sno_location_seq_all_species.tsv'
    conda:
        "../envs/python_new.yaml"
    script:
        "../scripts/python/merge_sno_location_species.py"