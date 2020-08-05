import os
import glob
import subprocess as sp
import ColorTextWriter

class Tagduster:

    def __init__(self, home_dir, threads, tagdust_sing, input_dir, rrna, extensions, seq_method):
        self.home_dir = home_dir
        self.threads = threads
        self.tagdust_sing = tagdust_sing
        self.input_dir = input_dir
        self.rrna_list = rrna
        self.extensions = extensions
        self.seq_method = seq_method

    def tagdust(self):

        outdir = os.path.join(self.home_dir, 'tagdust_out')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        r1_reads = sorted(glob.glob(self.input_dir + '*R1_trimmed_UMI.fastq'))
        r2_reads = sorted(glob.glob(self.input_dir + '*R2_trimmed_UMI.fastq'))

        ctw = ColorTextWriter.ColorTextWriter()

        print('\n' + ctw.CBEIGE + ctw.CBOLD + 'Running TagDust ...' + ctw.CEND)

        if self.seq_method == 'single':
            for i in r2_reads:
                print('\n' + ctw.CBEIGE + ctw.CBOLD + 'Tagdusting: ' + ctw.CBLUE + os.path.basename(i) + ctw.CBEIGE + ctw.CBOLD + ' ...' + ctw.CEND + '\n')

                output_file = outdir + '/' + os.path.basename(i).split('_R2_trimmed_UMI.fastq')[0] + '_tagdustout'

                command = [
                    'tagdust -t', self.threads,
                    '-1 R:N', '-o', output_file,
                    '-ref', self.rrna_list,'-fe 0', i
                ]

                command = ' '.join(command)
                sp.check_call(command, shell=True)

        elif self.seq_method == 'paired':
            for (i, j) in zip(r1_reads, r2_reads):
                print('\n' + ctw.CBEIGE + ctw.CBOLD + 'Tagdusting: ' + ctw.CBLUE + os.path.basename(i) + ' and ' + os.path.basename(j) + ctw.CBEIGE + ctw.CBOLD + ' ...' + ctw.CEND + '\n')

                output_file = outdir + '/' + os.path.basename(i).split('_R1_trimmed_UMI.fastq')[0] + '_tagdustout'

                command = [
                    'tagdust -t', self.threads,
                    '-1 R:N', '-o', output_file,
                    '-ref', self.rrna_list, '-fe 0', i, j
                ]

                command = ' '.join(command)
                sp.check_call(command, shell=True)

        print('\n' + ctw.CRED + ctw.CBOLD + 'Running TagDust Completed!!!' + ctw.CEND + '\n')
