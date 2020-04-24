import os
import glob
import subprocess as sp
import ColorTextWriter

class StarAligner:

    def __init__(self, home_dir, input_dir, threads, genome_dir, extensions, genes_gtf, seq_method):
        self.home_dir = home_dir
        self.input_dir = input_dir
        self.threads = threads
        self.genome_dir = genome_dir
        self.extensions = extensions
        self.genes_gtf = genes_gtf
        self.seq_method = seq_method

    def aligner(self):

        outdir = os.path.join(self.home_dir, 'star_aligned')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        #### Sequence Alignment using STAR
        se_reads = sorted(glob.glob(self.input_dir + '*tagdustout.fq'))
        r1_reads = sorted(glob.glob(self.input_dir + '*tagdustout_READ1.fq'))
        r2_reads = sorted(glob.glob(self.input_dir + '*tagdustout_READ2.fq'))

        ctw = ColorTextWriter.ColorTextWriter()

        print('\n' + ctw.CBEIGE + ctw.CBOLD + 'Running Star Aligner ...' + ctw.CEND)

        if self.seq_method == 'single':
            for i in se_reads:

                output_file = outdir + '/' + os.path.basename(i).split('tagdustout.fq')[0] + 'aligned' + self.extensions[4]

                if i.endswith('_tagdustout.fq'):
                    print('\n' + ctw.CBEIGE + ctw.CBOLD + 'Mapping: ' + ctw.CBLUE + os.path.basename(i) + ctw.CBEIGE + ctw.CBOLD + ' ...' + ctw.CEND + '\n')

                    command = [
                        'STAR --runThreadN', self.threads,
                        '--genomeDir', self.genome_dir,
                        '--outFileNamePrefix',
                        output_file.split(self.extensions[4])[0], '>', output_file,
                        '--outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999',
                        '--outReadsUnmapped Fastx --outSJfilterReads Unique --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1',
                        '--outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate',
                        '--readFilesIn', i
                    ]

                    command = ' '.join(command)
                    sp.check_call(command, shell=True)

        elif self.seq_method == 'paired':
            for (i, j) in zip(r1_reads, r2_reads):

                output_file = outdir + '/' + os.path.basename(i).split('tagdustout_READ1.fq')[0] + 'aligned' + self.extensions[4]

                if i.endswith('_tagdustout_READ1.fq') and j.endswith('_tagdustout_READ2.fq'):
                    print('\n' + ctw.CBEIGE + ctw.CBOLD + 'Mapping: ' + ctw.CBLUE + os.path.basename(i) + ' and ' + os.path.basename(j) + ctw.CBEIGE + ctw.CBOLD + ' ...' + ctw.CEND + '\n')

                    command = [
                        'STAR --runThreadN', self.threads,
                        '--genomeDir', self.genome_dir,
                        '--outFileNamePrefix',
                        output_file.split(self.extensions[4])[0], '>', output_file,
                        '--outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999',
                        '--outReadsUnmapped Fastx --outSJfilterReads Unique --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1',
                        '--outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate',
                        '--readFilesIn', i, j
                    ]

                    command = ' '.join(command)
                    sp.check_call(command, shell=True)

        print('\n' + ctw.CRED + ctw.CBOLD + 'Star Alignment Completed!!!' + ctw.CEND + '\n')

        #### Mapping quality control using Qualimap
        bam_files = sorted(glob.glob(outdir + '/' + '*.bam'))

        for i in bam_files:
            command = [
                'qualimap rnaseq',
                '-bam', i,
                '-gtf', self.genes_gtf
            ]

            command = ' '.join(command)
            sp.check_call(command, shell=True)

            command = [
                'qualimap bamqc',
                '-bam', i,
                '-gff', self.genes_gtf,
                '-sd -c'
            ]

            command = ' '.join(command)
            sp.check_call(command, shell=True)

        print('\n' + ctw.CRED + ctw.CBOLD + 'Alignment Quality Assessment Completed!!!' + ctw.CEND + '\n')