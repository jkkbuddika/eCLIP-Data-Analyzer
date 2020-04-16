import os
import glob
import subprocess as sp
import ColorTextWriter

class UMITools:

    def __init__(self, home_dir, input_dir, extensions, seq_method):
        self.home_dir = home_dir
        self.input_dir = input_dir
        self.extensions = extensions
        self.seq_method = seq_method

    def extract_UMI(self):

        outdir = os.path.join(self.home_dir, 'umi_extract')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        ctw = ColorTextWriter.ColorTextWriter()

        print(ctw.CBEIGE + ctw.CBOLD + 'Extracting UMIs ...' + ctw.CEND + '\n')

        r1_reads = sorted(glob.glob(self.input_dir + '*R1_unique_trimmed.fastq'))
        r2_reads = sorted(glob.glob(self.input_dir + '*R2_unique_trimmed.fastq'))

        if self.seq_method == 'single':
            for i in r2_reads:
                print(ctw.CBEIGE + ctw.CBOLD + 'UMI Extraction: ' + ctw.CBLUE + os.path.basename(i) + ctw.CBEIGE + ctw.CBOLD + ' ...' + ctw.CEND + '\n')

                output_file = outdir + '/' + os.path.basename(i).split('.fastq')[0] + '_UMI' + self.extensions[0]

                param = [
                    'umi_tools extract --extract-method=string',
                    '--stdin=' + i, '--bc-pattern=NNNNNNNNNN',
                    '--stdout=' + output_file,
                    '-L', output_file.split('_R2')[0] + '_UMI_extract.log'
                ]

                command = ' '.join(param)
                sp.check_call(command, shell=True)

        elif self.seq_method == 'paired':
            for (i, j) in zip(r1_reads, r2_reads):
                print('\n' + ctw.CBEIGE + ctw.CBOLD + 'UMI Extraction: ' + ctw.CBLUE + os.path.basename(i) + ' and ' + os.path.basename(j) + ctw.CBEIGE + ctw.CBOLD + ' ...' + ctw.CEND + '\n')

                output_file_R1 = outdir + '/' + os.path.basename(i).split('.fastq')[0] + '_UMI' + self.extensions[0]
                output_file_R2 = outdir + '/' + os.path.basename(j).split('.fastq')[0] + '_UMI' + self.extensions[0]

                command = [
                    'umi_tools extract --extract-method=string',
                    '--stdin=' + j, '--stdout=' + output_file_R2,
                    '--bc-pattern=NNNNNNNNNN', '--read2-in=' + i, '--read2-out=' + output_file_R1,
                    '-L', output_file_R1.split('_R1')[0] + '_UMI_extract.log'
                ]

                command = ' '.join(command)
                sp.check_call(command, shell=True)

        print('\n' + ctw.CRED + ctw.CBOLD + 'UMI Extraction Completed!!!' + ctw.CEND + '\n')

    def dedup(self):

        outdir = os.path.join(self.home_dir, 'umi_dedup')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        bam_list = sorted(glob.glob(self.input_dir + '*.bam'))

        ctw = ColorTextWriter.ColorTextWriter()

        for i in bam_list:

            print(ctw.CRED + 'Removing duplicates: ' + ctw.CBLUE + os.path.basename(i) + ctw.CRED + ' ...' + ctw.CEND + '\n')

            output_file = outdir + '/' + os.path.basename(i).split('.bam')[0] + '_dupRm' + self.extensions[4]

            command = [
                'umi_tools dedup',
                '-I', i
            ]

            if self.seq_method == 'paired': command.extend(['--paired'])

            command.extend([
                '-S', output_file,
                '--output-stats=' + output_file.split('.bam')[0]
            ])

            command = ' '.join(command)
            sp.check_call(command, shell=True)

        print(ctw.CBEIGE + ctw.CBOLD + 'Deduplicate Removal Completed!!!' + ctw.CEND + '\n')
