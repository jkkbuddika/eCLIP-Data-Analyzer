import os
import glob
import subprocess as sp
import ColorTextWriter

class BigWigFileMaker():

    def __init__(self, home_dir, input_dir, extensions):
        self.home_dir = home_dir
        self.input_dir = input_dir
        self.extensions = extensions

    def bigwig(self):

        outdir = os.path.join(self.home_dir, 'bedgraphs')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        bam_list = sorted(glob.glob(self.input_dir + '*.bam'))

        ctw = ColorTextWriter.ColorTextWriter()

        for i in bam_list:

            print('\n' + ctw.CRED + 'Generating BigWig file for ' + ctw.CBLUE + os.path.basename(i) + ctw.CRED + ' ...' + ctw.CEND + '\n')

            output_file = outdir + '/' + os.path.basename(i).split('.bam')[0] + self.extensions[5]

            command = [
                'samtools index', i
            ]

            command = ' '.join(command)
            sp.check_call(command, shell=True)

            command = [
                'bamCoverage --normalizeUsing CPM --binSize 1',
                '-b', i,
                '--outFileFormat bigwig',
                '-o', output_file
                ]

            command = ' '.join(command)
            sp.check_call(command, shell=True)