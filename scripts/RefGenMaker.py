import os
import subprocess as sp
import ColorTextWriter

class RefGenMaker:

    def __init__(self, home_dir, threads, genome_fasta, genes_gtf):
        self.home_dir = home_dir
        self.threads = threads
        self.genome_fasta = genome_fasta
        self.genes_gtf = genes_gtf

    def refgen(self):

        outdir = os.path.join(self.home_dir, 'star_genome')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        ctw = ColorTextWriter.ColorTextWriter()

        print(ctw.CBEIGE + ctw.CBOLD + 'Creating the STAR Reference Genome ...' + ctw.CEND + '\n')

        param = [
            'STAR --runThreadN', self.threads,
            '--runMode genomeGenerate --genomeSAindexNbases 12',
            '--genomeDir', outdir + '/',
            '--genomeFastaFiles', self.genome_fasta,
            '--sjdbGTFfile', self.genes_gtf
        ]

        command = ' '.join(param)
        sp.check_call(command, shell=True)

        print('\n' + ctw.CBEIGE + ctw.CBOLD + 'Reference Genome Created!!!' + ctw.CEND)