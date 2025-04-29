rule download_blat:
    """ Download BLAT to find coordinates 
        of a given sequence in a genome 
        (in fasta). Must enter sudo password at 
        a given point in the download."""
    output:
        tmp_file = "data/references/blat/blat_test.txt"
    params:
        blat = config['download']['blat']
    shell:
        "mkdir -p ./data/references/blat/ && "
        "rsync -aP {params.blat} ./data/references/blat/ &> {output.tmp_file}"

rule download_rfam_covariance_models:
    """ Download the Rfam library of covariance models that 
        will be used by Infernal."""
    output:
        rfam_cm = 'data/references/RFam.cm'
    params:
        link = config['download']['rfam_cm']
    shell:
        "wget {params.link} && "
        "gunzip Rfam.cm.gz && mv Rfam.cm {output.rfam_cm}"

rule download_genome:
    """ Download fasta of genome per species from Ensembl"""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["leishmania_major", 
              "dictyostelium_discoideum", "giardia_lamblia", "arabidopsis_thaliana", 
              "oryza_sativa", "aspergillus_fumigatus", "neurospora_crassa", "candida_albicans"], remove=True)
    params:
        link = "ftp://ftp.ensembl.org/pub/release-108/fasta/{species}/dna/*dna.toplevel.fa.gz",
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_mammal_genome:
    """ Download the reference genome (fasta file) of human and mouse
        from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species_tgirt'], ["mus_musculus",    
              "homo_sapiens"])
    params:
        link = "ftp://ftp.ensembl.org/pub/release-108/fasta/{species}/dna/*dna.primary_assembly.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_yeast_genome:
    """ Download the reference genome (fasta file) of S. pombe, 
        S. cerevisiae, N. crassa, C. albicans and A. fumigatus from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["saccharomyces_cerevisiae",
              "aspergillus_fumigatus", "neurospora_crassa", "candida_albicans"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-55/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_tetrahymena_genome:
    """ Download fasta of Tetrahymena thermophila genome per species from Zenodo/TGD"""
    output:
        genome = 'data/references/genome_fa/tetrahymena_thermophila_genome.fa'
    params:
        link = config['download']['t_thermophila_genome']
    shell:
        "wget -O {output.genome} {params.link}"

rule download_other_protist_genome:
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["leishmania_major", 
              "dictyostelium_discoideum", "giardia_lamblia"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/protists/release-55/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_plant_genome:
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species'], ["arabidopsis_thaliana", "oryza_sativa"])
    params:
        link = "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-55/fasta/{species}/dna/*dna.toplevel.fa.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.genome}"

rule download_o_tauri_genome:
    """ Download fasta of Ostreococcus tauri genome per species from Zenodo 
        (originally by concatenating the 20 chromosomes sequences from NCBI, 
        id ranging from NC_014426.2 to NC_014445.2)."""
    output:
        genome = 'data/references/genome_fa/ostreococcus_tauri_genome.fa'
    params:
        link = config['download']['o_tauri_genome']
    shell:
        "wget -O {output.genome} {params.link}"

rule download_mouse_gtf:
    """ Download the annotation of the mouse genome (gtf)
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/mus_musculus.gtf'
    params:
        link = config['download']['mouse_gtf']
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_yeast_gtf:
    """ Download the reference genome (fasta file) of S. cerevisiae, 
        S. pombe, A. fumigatus, N. crassa and C. albicans
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ["saccharomyces_cerevisiae", 
                "aspergillus_fumigatus", "neurospora_crassa", "candida_albicans"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-55/gtf/{species}/*5.gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_human_gtf:
    """ Download gtf of human genome from Zenodo. Remove trailing tabs."""
    output:
        gtf = 'data/references/gtf/homo_sapiens.gtf'
    params:
        link = config['download']['human_gtf']
    shell:
        "wget -O {output.gtf} {params.link} && "
        "sed -i 's/[\t]$//g' {output.gtf}"

rule download_tetrahymena_gtf:
    """ Download gtf of Tetrahymena thermophila genome from Zenodo. Remove trailing tabs."""
    output:
        gtf = 'data/references/gtf/tetrahymena_thermophila.gtf'
    params:
        link = config['download']['tetrahymena_gtf']
    shell:
        "wget -O {output.gtf} {params.link} && "
        "sed -i 's/[\t]$//g' {output.gtf}"

rule download_protist_gtf:
    """ Download the annotation (gtf) of different protists 
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ["leishmania_major", 
              "dictyostelium_discoideum", "giardia_lamblia", "plasmodium_falciparum"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/protists/release-55/gtf/{species}/*[^chr].gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"
    
rule download_plant_gtf:
    """ Download the annotation (gtf) of different plants
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ["oryza_sativa", "arabidopsis_thaliana"])
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/plants/release-55/gtf/{species}/*[^chr].gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_animal_gtf:
    """ Download the annotation (gtf) of different animals
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/{species}.gtf'
    wildcard_constraints:
        species=join_list(config['species'], ['macaca_mulatta', 'ornithorhynchus_anatinus', 'danio_rerio', 
                    'gallus_gallus', 'caenorhabditis_elegans', 'drosophila_melanogaster'])
    params:
        link = "ftp://ftp.ensembl.org/pub/release-108/gtf/{species}/*[^chrabinitio].gtf.gz"
    shell:
        "wget -O temp.gz {params.link} && "
        "gunzip temp.gz && "
        "mv temp {output.gtf}"

rule download_rnacentral_ncRNA:
    """ Download rnacentral bed file of all ncRNAs per species 
        (will be useful for tRNAs, snRNAs and pre-miRNA). 
        **Ostreococcus tauri is not present in rnacentral """
    output:
        bed = 'data/references/rnacentral/{species}.bed'
    wildcard_constraints:
        species=join_list(config['species']+config['species_tgirt'], 
                ["ostreococcus_tauri", "schizosaccharomyces_pombe"], remove=True)
    params:
        link = config['download']['rnacentral'] 
    shell:
        "wget -O temp_rnacentral_{wildcards.species}.gz {params.link}{wildcards.species}.*.bed.gz && "
        "gunzip temp_rnacentral_{wildcards.species}.gz && "
        "mv temp_rnacentral_{wildcards.species} {output.bed}"

rule download_gallus_gallus_gff:
    """ Download from rnacentral the latest version of Gallus gallus gff (bGalGal1/GRCg7b) to 
        convert the location of the old version (which is used in RNAcentral)"""
    output:
        gff = "data/references/rnacentral/gallus_gallus_bGalGal1_GRCg7b.gff3"  
    params:
        link = config['download']['gallus_gallus_gff']
    shell:
        "wget -O gff_gallus.gz {params.link} && "
        "gunzip gff_gallus.gz && mv gff_gallus {output.gff}"

rule download_eukaryote_HACA_snoRNAs:
    """ Download all H/ACA present in RNAcentral (v21) for all eukaryotes. It was done 
        manually by searching for snoRNAs and then selecting H/ACA snoRNAs (89 489 sequences). """
    output:
        HACA_fa = 'data/references/HACA/HACA_eukaryotes.fa'  # to put on Zenodo
    params:
        link = config['download']['haca_rnacentral']
    shell:
        "wget -O {output.HACA_fa} {params.link}"


rule download_rfam_clans:
    """ Download the clans in RFam, i.e. ~superfamilies (Release 14.9)"""
    output:
        df = 'data/references/rfam/fam_clans_families.tsv' 
    params:
        link = 'https://ftp.ebi.ac.uk/pub/databases/Rfam/14.9/database_files/clan_membership.txt.gz'  
    shell:
        "wget -O {output.df}.gz {params.link} && "
        "gunzip {output.df}.gz"

rule download_rDNA:
    """ Download ribosomal DNA fasta file per species 
        (it's needed to run snoscan). """ 
    output:
        rDNA_fa = 'data/references/rDNA/rDNA_species.fa'
    params:
        link = config['download']['rDNA'] 
    shell:
        "wget -O {output.rDNA_fa} {params.link}"

rule download_all_genes_pombe:
    """ Download from Pombase all genes and their description for S. pombe. 
        This is to retrieve all potential genes known as snoRNAs."""
    output:
        df = 'data/references/s_pombe/all_genes_s_pombe.tsv'
    params:
        link = config['download']['all_genes_pombe'] 
    shell:
        "wget -O {output.df} {params.link}"

rule download_CD_pombe_review:
    """ From Fafard-Couture et al. 2024 RNA Biology, get the 
        snoRNA type for snoRNAs in S. pombe."""
    output:
        df = 'data/references/s_pombe/cd_s_pombe_review.tsv'
    params:
        link = config['download']['cd_s_pombe_review']
    shell:
        "wget -O {output.df} {params.link}"

rule download_sno_type_info:
    """ Download snoRNA type per Rfam family. The snotype per Rfam 
        family was obtained manually by saving (right-click) the 
        Rfam webpage displaying all C/D (and then H/ACA) snoRNA 
        families and parsing the rfam family id only using the following command: 
        grep ">RF" ~/Downloads/HACA_RFAM_families_2024.html | sed -E 's/.*">//g; s/<.*//g'"""
    output:
        sno_type_rfam = 'data/references/snoRNA_type_rfam_families.tsv'
    params:
        link = config['download']['sno_type_rfam']
    shell:
        "wget -O {output.sno_type_rfam} {params.link}"

rule download_sno_litt:
    """ Download snoRNA fastas from different species obtained from literature 
        curation of expressed C/D box snoRNAs."""
    output:
        fa = "data/sno_literature_processing/{sno_fasta}.fa"
    params:
        link = "https://zenodo.org/records/15276930/files/"
    shell:
        "wget -O {output.fa} {params.link}{wildcards.sno_fasta}.fa"

rule download_CD_bed_rnacentral:
    """ Download bed of C/D snoRNAs in different species manually selected 
        from RNAcentral."""
    output:
        bed = "data/references/rnacentral/sno_type_df/{species}_cd_id.bed"
    params:
        link = "https://zenodo.org/records/15276930/files/"
    wildcard_constraints:
        species=join_list(config['species'], ['drosophila_melanogaster', 
            'gallus_gallus', 'macaca_mulatta', 'tetrahymena_thermophila'])
    shell:
        "wget -O {output.bed} {params.link}{wildcards.species}_cd_id.bed"

rule download_pombe_genome_negative_strand:
    """ Download genome sequence of - strand of S. pombe genome (manually 
        reverse transcribed) that is needed for snoreport2 to work."""
    output:
        fa = 'data/references/genome_fa/negative_strand/schizosaccharomyces_pombe_genome_negative_strand.fa'
    params:
        link = config['download']['pombe_negative_strand']
    shell:
        "wget -O {output.fa} {params.link}"

rule download_sno_type_df:
    """ Download snoRNA type info for different species."""
    output:
        human = 'data/references/sno_type_df/homo_sapiens_snotype.tsv',
        mouse =  'data/references/sno_type_df/mus_musculus_snotype.tsv',
        ttherm_gene = 'data/references/sno_type_df/tetrahymena_thermophila_snotype_good_genome.tsv',
        ttherm = 'data/references/sno_type_df/tetrahymena_thermophila_snotype.tsv',
        yeast = 'data/references/sno_type_df/saccharomyces_cerevisiae_snotype.tsv',
        pombe = 'data/references/sno_type_df/schizosaccharomyces_pombe_snotype.tsv',
        danio = 'data/references/sno_type_df/danio_rerio_snotype.tsv',
        macaca = 'data/references/sno_type_df/macaca_mulatta_snotype.tsv',
        gallus = 'data/references/sno_type_df/gallus_gallus_snotype.tsv'
    params:
        link = "https://zenodo.org/records/15276930/files/",
        link2 = "https://zenodo.org/records/15306359/files/",
        human = "homo_sapiens_snotype_snoDB.tsv",
        mouse = "mus_musculus_snotype_rnacentral",
        ttherm_gene = "tetrahymena_thermophila_snotype_good_genome.tsv",
        ttherm = "tetrahymena_thermophila_snotype.tsv",
        yeast = "saccharomyces_cerevisiae_snotype_umass.tsv",
        pombe = "schizosaccharomyces_pombe_snotype.tsv",  # from rnacentral
        danio = "danio_rerio_snotype.tsv",  # from rnacentral
        macaca = "macaca_mulatta_snotype.tsv",  # from rnacentral
        gallus = "gallus_gallus_snotype.tsv",  # from rnacentral
    shell:
        "wget -O {output.human} {params.link}{params.human} && "
        "wget -O {output.mouse} {params.link}{params.mouse} && "
        "wget -O {output.ttherm_gene} {params.link}{params.ttherm_gene} && "
        "wget -O {output.ttherm} {params.link}{params.ttherm} && "
        "wget -O {output.yeast} {params.link}{params.yeast} && "
        "wget -O {output.pombe} {params.link2}{params.pombe} && "
        "wget -O {output.danio} {params.link2}{params.danio} && "
        "wget -O {output.macaca} {params.link2}{params.macaca} && "
        "wget -O {output.gallus} {params.link2}{params.gallus}"

rule download_RIP_results:
    """ Download NOP58 and Fib RIP-qPCR results."""
    output:
        nop58 = 'results/IP_cedric/NOP58_TAP_fold_enrichment.tsv',
        fib = 'results/IP_cedric/Fib_fold_enrichment.tsv'
    params:
        link = 'https://zenodo.org/records/15276930/files/'
    shell:
        "wget -O {output.nop58} {params.link}NOP58_TAP_fold_enrichment.tsv && "
        "wget -O {output.fib} {params.link}Fib_fold_enrichment.tsv"

rule download_sf3b3_snoRNAs:
    """ Download Sf3b3 snoRNA sequence and Newick file alignment results."""
    output:
        fa = 'data/references/sf3b3_snoRNAs_across_species.fa',
        nwk = 'data/references/sf3b3_snoRNAs_across_species.nwk'
    params:
        link = 'https://zenodo.org/records/15276930/files/sf3b3_snoRNAs_across_species'
    shell:
        "wget -O {output.fa} {params.link}.fa && "
        "wget -O {output.nwk} {params.link}.nwk"

rule download_tgirt_seq_tpm:
    """ Download TGIRT-Seq abundance table (in TPM) of all genes in different 
        species after running the TGIRT-Seq analysis pipeline."""
    output:
        tsv = "data/references/tgirt_seq_output/{sp}.tsv"
    params:
        link = "https://zenodo.org/records/15276930/files/"
    wildcard_constraints:
        sp = "{}".format("|".join(['drosophila_melanogaster_expression_status_snoRNAs', 'homo_sapiens_merged_tpm_w_biotype', 
            'mus_musculus_merged_tpm_w_biotype', 'saccharomyces_cerevisiae_merged_tpm_w_biotype', 
            'gallus_gallus_expressed_tpm_w_biotype', 'macaca_mulatta_expressed_tpm_w_biotype']))
    shell:
        "wget -O {output.tsv} {params.link}{wildcards.sp}.tsv"

rule download_human_bedgraph_1:
    """ Download Human bedgraph obtained from TGIRT-Seq datasets."""
    output:
        bg = 'results/bedgraph_TGIRT/homo_sapiens/{species}.bam_{sense}.bedgraph'
    params:
        link = 'https://zenodo.org/records/15277346/files/'
    wildcard_constraints:
        species = "{}".format("|".join(config['bedgraph']['human_1'])),
        sense = "{}".format("|".join(['fwd', 'rev']))
    shell:
        "wget -O {output.bg} {params.link}{wildcards.sample_human1}.bam_{wildcards.sense}.bedgraph"

rule download_human_bedgraph_2:
    """ Download Human bedgraph obtained from TGIRT-Seq datasets."""
    output:
        bg = 'results/bedgraph_TGIRT/homo_sapiens/{sample_human2}.bam_{sense}.bedgraph'
    params:
        link = 'https://zenodo.org/records/15277470/files/'
    wildcard_constraints:
        sample_human2 = "{}".format("|".join(config['bedgraph']['human_2'])),
        sense = "{}".format("|".join(['fwd', 'rev']))
    shell:
        "wget -O {output.bg} {params.link}{wildcards.sample_human2}.bam_{wildcards.sense}.bedgraph"

rule download_danio_bedgraph:
    """ Download Danio rerio bedgraph obtained from TGIRT-Seq datasets."""
    output:
        bg = 'results/bedgraph_TGIRT/danio_rerio/{sample_danio}.bam_{sense}.bedgraph'
    params:
        link = 'https://zenodo.org/records/15277346/files/'
    wildcard_constraints:
        sample_danio = "{}".format("|".join(config['bedgraph']['danio'])),
        sense = "{}".format("|".join(['fwd', 'rev']))
    shell:
        "wget -O {output.bg} {params.link}{wildcards.sample_danio}.bam_{wildcards.sense}.bedgraph"

rule download_droso_bedgraph:
    """ Download D.melanogaster bedgraph obtained from TGIRT-Seq datasets."""
    output:
        bg = 'results/bedgraph_TGIRT/drosophila_melanogaster/{sample_droso}.bam_{sense}.bedgraph'
    params:
        link = 'https://zenodo.org/records/15277346/files/'
    wildcard_constraints:
        sample_droso = "{}".format("|".join(config['bedgraph']['droso'])),
        sense = "{}".format("|".join(['fwd', 'rev']))
    shell:
        "wget -O {output.bg} {params.link}{wildcards.sample_droso}.bam_{wildcards.sense}.bedgraph"

rule download_gallus_bedgraph:
    """ Download Gallus gallus bedgraph obtained from TGIRT-Seq datasets."""
    output:
        bg = 'results/bedgraph_TGIRT/gallus_gallus/{sample_gallus}.bam_{sense}.bedgraph'
    params:
        link = 'https://zenodo.org/records/15277346/files/'
    wildcard_constraints:
        sample_gallus = "{}".format("|".join(config['bedgraph']['gallus'])),
        sense = "{}".format("|".join(['fwd', 'rev']))
    shell:
        "wget -O {output.bg} {params.link}{wildcards.sample_gallus}.bam_{wildcards.sense}.bedgraph"

rule download_macaca_bedgraph:
    """ Download Macaca mulatta bedgraph obtained from TGIRT-Seq datasets."""
    output:
        bg = 'results/bedgraph_TGIRT/macaca_mulatta/{sample_macaca}.bam_{sense}.bedgraph'
    params:
        link = 'https://zenodo.org/records/15277346/files/'
    wildcard_constraints:
        sample_macaca = "{}".format("|".join(config['bedgraph']['macaca'])),
        sense = "{}".format("|".join(['fwd', 'rev']))
    shell:
        "wget -O {output.bg} {params.link}{wildcards.sample_macaca}.bam_{wildcards.sense}.bedgraph"

rule download_ttherm_bedgraph:
    """ Download T. thermophila bedgraph obtained from TGIRT-Seq datasets."""
    output:
        bg = 'results/bedgraph_TGIRT/tetrahymena_thermophila/{sample_ttherm}.bam_{sense}.bedgraph'
    params:
        link = 'https://zenodo.org/records/15277479/files/'
    wildcard_constraints:
        sample_ttherm = "{}".format("|".join(config['bedgraph']['tetrahymena'])),
        sense = "{}".format("|".join(['fwd', 'rev']))
    shell:
        "wget -O {output.bg} {params.link}{wildcards.sample_ttherm}.bam_{wildcards.sense}.bedgraph"

rule download_plasmodium_bedgraph:
    """ Download P. falciparum bedgraph obtained from TGIRT-Seq datasets."""
    output:
        bg = 'results/bedgraph_TGIRT/plasmodium_falciparum/{sample_plasmo}.bam_{sense}.bedgraph'
    params:
        link = 'https://zenodo.org/records/15277479/files/'
    wildcard_constraints:
        sample_plasmo = "{}".format("|".join(config['bedgraph']['plasmodium'])),
        sense = "{}".format("|".join(['fwd', 'rev']))
    shell:
        "wget -O {output.bg} {params.link}{wildcards.sample_plasmo}.bam_{wildcards.sense}.bedgraph"

rule download_pombe_bedgraph:
    """ Download S. pombe bedgraph obtained from TGIRT-Seq datasets."""
    output:
        bg = 'results/bedgraph_TGIRT/schizosaccharomyces_pombe/{sample_pombe}.bam_{sense}.bedgraph'
    params:
        link = 'https://zenodo.org/records/15277479/files/'
    wildcard_constraints:
        sample_pombe = "{}".format("|".join(config['bedgraph']['pombe'])),
        sense = "{}".format("|".join(['fwd', 'rev']))
    shell:
        "wget -O {output.bg} {params.link}{wildcards.sample_pombe}.bam_{wildcards.sense}.bedgraph"

rule download_f1_score:
    """ Download snoBIRD's f1 score before and after training (useful for 
        learning curve figure generation)."""
    output:
        f1 = 'results/predictions/transformer/194/3e-5_3e-6_32_4_data_aug_1_ratio/transformer_2_classes_{stage}_trained_fold_{fold_num}_f1_score_per_epoch.tsv'
    params:
        link = 'https://zenodo.org/records/15277470/files/'
    wildcard_constraints:
        stage = "{}".format("|".join(['Before', 'LR_schedule']))
    shell:
        "wget -O {output.f1} {params.link}transformer_2_classes_{wildcards.stage}_trained_fold_{wildcards.fold_num}_f1_score_per_epoch.tsv"


rule download_shap_first_second_model:
    """ Download SHAP values computed for SnoBIRD's first and second model on 
        the test set examples accureatly predicted as C/D box snoRNAs."""
    output:
        shap_first = 'results/shap/snoBIRD/all_cd_shap_values.tsv',
        shap_second = 'results/shap/snoBIRD/all_cd_shap_values_sno_pseudo.tsv'
    params:
        link = 'https://zenodo.org/records/15276930/files/'
    shell:
        "wget -O {output.shap_first} {params.link}all_cd_shap_values.tsv && "
        "wget -O {output.shap_second} {params.link}all_cd_shap_values_sno_pseudo.tsv"

rule download_all_jobs_cluster:
    """ Download all the job specs that were run during SnoBIRD's development 
        phase to estimate CO2 emission of the whole process."""
    output:
        narval = 'data/all_jobs_narval.tsv',
        beluga = 'data/all_jobs_beluga.tsv'
    params:
        link = 'https://zenodo.org/records/15277346/files/'
    shell:
        "wget -O {output.narval} {params.link}all_jobs_narval.tsv && "
        "wget -O {output.beluga} {params.link}all_jobs_beluga.tsv"

rule download_snoBIRD_preds:
    """ Download SnoBIRD's prediction in different species."""
    output:
        danio = 'results/predictions/snoBIRD/final/danio_rerio/danio_rerio_s5_cs100.tsv',
        droso = 'results/predictions/snoBIRD/final/drosophila_melanogaster/drosophila_melanogaster_s5_cs10.tsv',
        gallus = 'results/predictions/snoBIRD/final/gallus_gallus/gallus_gallus_s5_cs10.tsv',
        macaca = 'results/predictions/snoBIRD/final/macaca_mulatta/macaca_mulatta_s5_cs10.tsv',
        human = 'results/predictions/snoBIRD/final/homo_sapiens/homo_sapiens_snoBIRD_s1_cs5.tsv',
        ttherm = 'results/predictions/snoBIRD/final/tetrahymena_thermophila/tetrahymena_thermophila_s5.tsv',
        plasmo = 'results/predictions/snoBIRD/final/plasmodium_falciparum/plasmodium_falciparum_s5_cs5.tsv',
        pombe = 'results/predictions/snoBIRD/final/schizosaccharomyces_pombe/snoBIRD_complete_predictions.tsv'
    params:
        link_danio = 'https://zenodo.org/records/15276930/files/danio_rerio_s5_cs100.tsv',
        link_droso = 'https://zenodo.org/records/15276930/files/drosophila_melanogaster_s5_cs10.tsv',
        link_gallus = 'https://zenodo.org/records/15276930/files/gallus_gallus_s5_cs10.tsv',
        link_macaca = 'https://zenodo.org/records/15276930/files/macaca_mulatta_s5_cs10.tsv',
        link_human = 'https://zenodo.org/records/15276930/files/homo_sapiens_snoBIRD_s1_cs5.tsv',
        link_ttherm = 'https://zenodo.org/records/15276930/files/tetrahymena_thermophila_s5.tsv',
        link_plasmo = 'https://zenodo.org/records/15276930/files/plasmodium_falciparum_s5_cs5.tsv',
        link_pombe = 'https://zenodo.org/records/15276930/files/snoBIRD_complete_predictions.tsv',
    shell:
        "wget -O {output.danio} {params.link_danio} && "
        "wget -O {output.droso} {params.link_droso} && "
        "wget -O {output.gallus} {params.link_gallus} && "
        "wget -O {output.macaca} {params.link_macaca} && "
        "wget -O {output.human} {params.link_human} && "
        "wget -O {output.ttherm} {params.link_ttherm} && "
        "wget -O {output.plasmo} {params.link_plasmo} && "
        "wget -O {output.pombe} {params.link_pombe}"
