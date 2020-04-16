import os
import glob
import subprocess as sp
import ColorTextWriter

class MultiQCRunner:

    def __init__(self, home_dir):
        self.home_dir = home_dir

    def multiqc(self):

        outdir = os.path.join(self.home_dir, 'MultiQC_Summary')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        dir_list = sorted(glob.glob(self.home_dir + '*/'))

        ctw = ColorTextWriter.ColorTextWriter()

        print(ctw.CBEIGE + ctw.CBOLD + 'Running MultiQC ...' + ctw.CEND + '\n')

        for i in dir_list:
            print('\n' + ctw.CBEIGE + ctw.CBOLD + 'Generating MultiQC Reports for ' + i + ' ...' + ctw.CEND + '\n')

            command = [
                'multiqc',
                i, '-o', outdir,
                '-n', i.split(self.home_dir)[1].split('/')[0] + '_MultiQC_Report.html'
            ]

            command = ' '.join(command)
            sp.check_call(command, shell=True)

        print('\n' + ctw.CBEIGE + ctw.CBOLD + 'Running MultiQC done!!!' + ctw.CEND + '\n')