import os
import glob
import subprocess as sp
import ColorTextWriter

class FastQCRunner:

    def __init__(self, home_dir, output_dir, input_dir):
        self.home_dir = home_dir
        self.output_dir = output_dir
        self.input_dir = input_dir

    def fastqc(self):

        outdir = os.path.join(self.home_dir, self.output_dir)
        if not os.path.isdir(outdir): os.mkdir(outdir)

        file_list = sorted(glob.glob(self.input_dir + '*R2.fastq'))

        ctw = ColorTextWriter.ColorTextWriter()

        print(ctw.CBEIGE + ctw.CBOLD + 'Running FastQC ...' + ctw.CEND)

        for i in file_list:

            print('\n')

            command = [
                'fastqc', i,
                '-o', outdir
            ]

            command = ' '.join(command)
            sp.check_call(command, shell=True)

        print('\n' + ctw.CBEIGE + ctw.CBOLD + 'Running FastQC done!!!' + ctw.CEND + '\n')
