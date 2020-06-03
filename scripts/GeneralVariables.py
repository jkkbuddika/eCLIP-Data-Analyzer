class GeneralVariables:

    ## General
    seq_method = 'single'

    ## Name of the Biomart rRNA list in add_mat directory
    rRNA_list = 'Dro_rRNA.txt'

    ## Link to the Biomart genome file of interest
    genome = 'ftp://ftp.ensembl.org/pub/release-100/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna_sm.toplevel.fa.gz'

    ## Link to the annotation file of interest (i.e., Drosophila annotation from ensembl FTP site)
    feature = 'ftp://ftp.ensembl.org/pub/release-100/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.100.gtf.gz'

    ## Strandedness of the experiment: '0', '1' or '2'
    stranded = '0'

    ## Include a list of features to be quantified
    diff_features = ['gene', 'CDS', 'five_prime_utr', 'three_prime_utr', 'exon']
