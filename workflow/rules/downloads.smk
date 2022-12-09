rule download_blat:
    """ Download BLAT to find coordinates 
        of a given sequence in a genome 
        (in fasta). Must enter sudo password at 
        a given point in the download."""
    output:
        tmp_file = "data/blat_test.txt"
    params:
        libpng = config['download']['libpng'],
        blat = config['download']['blat']
    shell:
        "cd data/references/ && wget {params.libpng} && "
        "tar -xvzf libpng-1.6.2.tar.gz && cd libpng-1.6.2/ && "
        "./configure --prefix=`pwd` && make && "
        "sudo make install && LIBPNGDIR=`pwd` && cd ../ && "
        "wget {params.blat} && unzip blatSrc35.zip && cd blatSrc/ && "
        "cp $LIBPNGDIR/png.h lib/ && cp $LIBPNGDIR/pngconf.h lib/ && "
        "cp $LIBPNGDIR/pnglibconf.h lib/ && echo $MACHTYPE && "
        "MACHTYPE=x86_64 && export MACHTYPE && mkdir -p ~/bin/$MACHTYPE && "
        "make && echo 'export MACHTYPE=x86_64' >> ~/.bashrc && "
        "echo 'export PATH=$PATH:~/bin/$MACHTYPE/blat' >> ~/.bashrc && "
        "source ~/.bashrc && cd ../../../ && touch {output.tmp_file}"

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
              "oryza_sativa"], remove=True)
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
    """ Download the reference genome (fasta file) of S. pombe and 
        S. cerevisiae from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/{species}_genome.fa'
    wildcard_constraints:
        species=join_list(config['species_tgirt'], ["saccharomyces_cerevisiae",
              "schizosaccharomyces_pombe"])
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
    """ Download the reference genome (fasta file) of S. cerevisiae 
        from ENSEMBL ftp servers."""
    output:
        gtf = 'data/references/gtf/saccharomyces_cerevisiae.gtf'
    params:
        link = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-55/gtf/saccharomyces_cerevisiae/*5.gtf.gz"
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