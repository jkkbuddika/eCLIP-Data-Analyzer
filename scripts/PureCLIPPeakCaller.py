import os
import glob
import subprocess as sp
import ColorTextWriter

class PureCLIPPeakCaller():

    def __init__(self, home_dir, input_dir, genome_fa, extensions):
        self.home_dir = home_dir
        self.input_dir = input_dir
        self.genome_fa = genome_fa
        self.extensions = extensions

    def peak_caller(self):

        outdir = os.path.join(self.home_dir, 'pure_clip_pc')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        file_list = sorted(glob.glob(self.input_dir + '*_UV*.bam'))
        input_bam = sorted(glob.glob(self.input_dir + '*_In*.bam'))

        ctw = ColorTextWriter.ColorTextWriter()

        for i in file_list:
            print(ctw.CRED + 'Peak Calling: ' + ctw.CBLUE + os.path.basename(i) + ctw.CRED + ' ...' + ctw.CEND + '\n')

            output_file = outdir + '/' + os.path.basename(i).split('.bam')[0] + '_peakCall' + self.extensions[6]

            param = [
                'pureclip',
                '-i', i, '-bai', i + '.bai',
                '-g', self.genome_fa, '-nt 10',
                '-o', output_file,
                '-ibam', input_bam[0], '-ibai', input_bam[0] + '.bai'
            ]

            command = ' '.join(param)
            sp.check_call(command, shell=True)

        print(ctw.CBEIGE + ctw.CBOLD + 'Peak Calling Completed!!!' + ctw.CEND)