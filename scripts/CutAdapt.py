import os
import glob
import subprocess as sp
import ColorTextWriter

class CutAdapt:

    def __init__(self, home_dir, input_dir, adapter, r2_adapter_seq, extensions, seq_method):
        self.home_dir = home_dir
        self.input_dir = input_dir
        self.adapter = adapter
        self.r2_adapter_seq = r2_adapter_seq
        self.extensions = extensions
        self.seq_method = seq_method

    def cutadapt(self):

        outdir_1 = os.path.join(self.home_dir, 'cutadapt')
        if not os.path.isdir(outdir_1): os.mkdir(outdir_1)

        ctw = ColorTextWriter.ColorTextWriter()

        print(ctw.CBEIGE + ctw.CBOLD + 'Running CutAdapt ...' + ctw.CEND + '\n')

        r1_reads = sorted(glob.glob(self.input_dir + '*R1.fastq'))
        r2_reads = sorted(glob.glob(self.input_dir + '*R2.fastq'))

        x = 0
        while x < len(self.adapter):

            for (i, j) in zip(r1_reads, r2_reads):
                if i.endswith(self.adapter[x] + '_R1.fastq') and j.endswith(self.adapter[x] + '_R2.fastq'):
                    print('\n' + ctw.CBEIGE + ctw.CBOLD + 'CutAdapting: ' + ctw.CBLUE + os.path.basename(i) + ' and ' + os.path.basename(j) + ctw.CBEIGE + ctw.CBOLD + ' ...' + ctw.CEND + '\n')

                    trim_file_R1 = outdir_1 + '/' + os.path.basename(i).split(os.path.splitext(i)[1])[0] + '_trimmed' + self.extensions[0]
                    trim_file_R2 = outdir_1 + '/' + os.path.basename(j).split(os.path.splitext(i)[1])[0] + '_trimmed' + self.extensions[0]

                    command = [
                        'cutadapt -f fastq --match-read-wildcards --times 2 -e 0.1 -O 1 --quality-cutoff 5 -m 20',
                        self.r2_adapter_seq[x],
                        '-o', trim_file_R1, '-p', trim_file_R2, i, j,
                        '>', trim_file_R1.split('_R1_trimmed')[0] + '_trim.matrics' + self.extensions[3]
                    ]

                    command = ' '.join(command)
                    sp.check_call(command, shell=True)

            x = x + 1

        print('\n' + ctw.CRED + ctw.CBOLD + 'CutAdapt Trimming Completed!!!' + ctw.CEND + '\n')
