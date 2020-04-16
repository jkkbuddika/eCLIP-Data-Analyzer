import os
import GeneralVariables

class CommonVariables:
    gv = GeneralVariables.GeneralVariables()

    ## General
    home_dir = os.path.dirname(os.getcwd()) + '/'
    raw_sequences_dir = home_dir + 'raw_sequences/'
    add_mat = home_dir + 'add_mat/'
    Threads = '8'
    extensions = ['.fastq', '.sam', '.csv', '.txt', '.bam', '.bw', '.bed']
    summary_dir = home_dir + 'summary_files/'

    ## Potential adaptor sequences
    L19 = ''
    X1A = '-a ATATAGGNNNNNAGATCGGAAGAGCGTCGTGTAG'
    X1B = '-a AATAGCANNNNNAGATCGGAAGAGCGTCGTGTAG'
    X2A = '-a AAGTATANNNNNAGATCGGAAGAGCGTCGTGTAG'
    X2B = '-a AGAAGATNNNNNAGATCGGAAGAGCGTCGTGTAG'
    A01 = '-a ATTGCTTAGATCGGAAGAGCGTCGTGTAG'
    B06 = '-a ACAAGCCAGATCGGAAGAGCGTCGTGTAG'
    r2_adapter = ['L19', 'X1A', 'X1B', 'X2A', 'X2B', 'A01', 'B06']
    r2_adapter_seq = [L19, X1A, X1B, X2A, X2B, A01, B06]

    ## TagDust Variables
    tagdust_singu = add_mat + 'gostripes.simg'
    tagdust_out = home_dir + 'tagdust_out/'
    rRNA_path = add_mat + gv.rRNA_list

    ## Cutadapt Variables
    cutadapt_dir = home_dir + 'cutadapt/'
    dupCollapse_dir = home_dir + 'dupCollapsed/'


    ## UMI Tools Variables
    umi_extract = home_dir + 'umi_extract/'
    umi_dedup = home_dir + 'umi_dedup/'

    ## Variables for FastQC
    fastqc_raw = 'fastqc_raw'

    ## Genome sequence and annotation
    genome_file = os.path.basename(gv.genome)
    genome_path = os.path.dirname(gv.genome) + '/'
    genome_dir_name = 'genome'
    genome_dir = home_dir + 'genome/'
    genome_fa = genome_dir + os.path.splitext(genome_file)[0]
    feature_file = os.path.basename(gv.feature)
    feature_path = os.path.dirname(gv.feature) + '/'
    feature_dir_name = 'genome_feature'
    feature_dir = home_dir + 'genome_feature/'
    genes_gtf = feature_dir + os.path.splitext(feature_file)[0]

    ## STAR Alignment
    ref_genome = home_dir + 'star_genome/'
    star_aligned = home_dir + 'star_aligned/'

    ## Sam Tools Sorting
    sam_sorted = home_dir + 'sam_sorted/'
    sam_retrieved_r2 = home_dir + 'sam_retrieved_r2/'

    ## DeepTools BigWig Files
    bigwig_files = home_dir + 'bigwig_files/'

    ## FeatureCounts
    fc_output = home_dir + 'feature_counts'

    ## PureCLIP Peak Caller
    pure_clip_pc = home_dir + 'pure_clip_pc/'