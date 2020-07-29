import os
import UserDefinedOptions

class CommonOptions:
    gv = UserDefinedOptions.UserDefinedOptions()

    ## General
    home_dir = os.path.dirname(os.getcwd()) + '/'
    raw_sequences_dir = home_dir + 'raw_sequences/'
    add_mat = home_dir + 'add_mat/'
    Threads = '8'
    extensions = ['.fastq', '.sam', '.csv', '.txt', '.bam', '.bw', '.bed']
    summary_dir = home_dir + 'summary_files/'
    mode = ['single', 'paired']
    seq_method = mode[int(gv.runMode)]

    ## Potential Adapter Sequences
    L19 = ''
    X1A = '-A ATATAGGNNNNNAGATCGGAAGAGCGTCGTGTAG'
    X1B = '-A AATAGCANNNNNAGATCGGAAGAGCGTCGTGTAG'
    X2A = '-A AAGTATANNNNNAGATCGGAAGAGCGTCGTGTAG'
    X2B = '-A AGAAGATNNNNNAGATCGGAAGAGCGTCGTGTAG'
    A01 = '-A ATTGCTTAGATCGGAAGAGCGTCGTGTAG'
    B06 = '-A ACAAGCCAGATCGGAAGAGCGTCGTGTAG'
    r2_adapter = ['L19', 'X1A', 'X1B', 'X2A', 'X2B', 'A01', 'B06']
    r2_adapter_seq = [L19, X1A, X1B, X2A, X2B, A01, B06]

    ## TagDust Variables
    tagdust_singu = add_mat + 'gostripes.simg'
    tagdust_out = home_dir + 'tagdust_out/'
    rRNA_path = add_mat + gv.species + '_rRNA.txt'

    ## Cutadapt Variables
    cutadapt_dir = home_dir + 'cutadapt/'
    dupCollapse_dir = home_dir + 'dupCollapsed/'

    ## UMI Tools Variables
    umi_options = ['NNNNN', 'NNNNNNNNNN']
    umi_seq = umi_options[int(gv.umi)]
    umi_extract = home_dir + 'umi_extract/'
    umi_dedup = home_dir + 'umi_dedup/'

    ## Variables for FastQC
    fastqc_raw = 'fastqc_raw'

    ## Support Genomes
    dm_genome = 'ftp://ftp.ensembl.org/pub/release-100/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna_sm.toplevel.fa.gz'
    ce_genome = 'ftp://ftp.ensembl.org/pub/release-100/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz'
    hs_genome = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz'
    mm_genome = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz'
    dr_genome = 'ftp://ftp.ensembl.org/pub/release-100/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.toplevel.fa.gz'
    if gv.species == 'custom': custom_genome = input('Ensembl url to the genome to download: ')

    ## Support Annotations
    dm_annotation = 'ftp://ftp.ensembl.org/pub/release-100/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.100.gtf.gz'
    ce_annotation = 'ftp://ftp.ensembl.org/pub/release-100/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.100.gtf.gz'
    hs_annotation = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz'
    mm_annotation = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz'
    dr_annotation = 'ftp://ftp.ensembl.org/pub/release-100/gtf/danio_rerio/Danio_rerio.GRCz11.100.gtf.gz'
    if gv.species == 'custom': custom_annotation = input('Ensembl url to the annotation to download: ')

    ## Genome Sequence and Annotation Details
    genome_file = os.path.basename(eval(gv.species + '_genome'))
    genome_path = os.path.dirname(eval(gv.species + '_genome')) + '/'
    genome_dir_name = 'genome'
    genome_dir = home_dir + 'genome/'
    genome_fa = genome_dir + os.path.splitext(genome_file)[0]
    feature_file = os.path.basename(eval(gv.species + '_annotation'))
    feature_path = os.path.dirname(eval(gv.species + '_annotation')) + '/'
    feature_dir_name = 'genome_feature'
    feature_dir = home_dir + 'genome_feature/'
    genes_gtf = feature_dir + os.path.splitext(feature_file)[0]

    ## STAR Alignment
    ref_genome = home_dir + 'star_genome/'
    star_aligned = home_dir + 'star_aligned/'

    ## Sam Tools Sorting
    sam_sorted = home_dir + 'sam_sorted/'
    sam_retrieved_r2 = home_dir + 'sam_retrieved_r2/'

    ## BedTools and Crosslink Data
    crosslink_data = home_dir + 'crosslink_data/'

    ## DeepTools BigWig Files
    bigwig_files = home_dir + 'bigwig_files/'

    ## FeatureCounts
    fc_output = home_dir + 'feature_counts'
    diff_features = ['gene', 'exon']

    ## PureCLIP Peak Caller
    pure_clip_pc = home_dir + 'pure_clip_pc/'
