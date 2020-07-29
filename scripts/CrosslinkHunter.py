import os
import glob
import pandas as pd
import subprocess as sp
import ColorTextWriter

class CrosslinkHunter:

    def __init__(self, home_dir, input_dir, genome_fasta, threads, extensions):
        self.home_dir = home_dir
        self.input_dir = input_dir
        self.genome_fasta = genome_fasta
        self.threads = threads
        self.extensions = extensions

    def crosslink(self):

        outdir = os.path.join(self.home_dir, 'crosslink_data')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        bam_list = sorted(glob.glob(self.input_dir + '*UV*.bam'))

        ctw = ColorTextWriter.ColorTextWriter()

        #### Generate a chromosome length file using SAMTools
        sp.check_call(' '.join(['samtools faidx', self.genome_fasta]), shell=True)
        fai_file = self.genome_fasta + '.fai'

        for i in bam_list:

            print('\n' + ctw.CRED + 'Detection of cross-linked nucleotides: ' + ctw.CBLUE + os.path.basename(i) + ctw.CRED + ' ...' + ctw.CEND + '\n')

            output_file = outdir + '/' + os.path.basename(i).split('.bam')[0] + '_shifted.bed'

            #### Convert bam to bed and shift intervals 1bp upstream to identify the cross-linked nucleotide
            bamtoshiftedbed = [
                'bamToBed -i', i, '-bed12', '|',
                'bedtools shift -m 1 -p -1 -g', fai_file, '>', output_file
            ]

            bamtoshiftedbed = ' '.join(bamtoshiftedbed)
            sp.check_call(bamtoshiftedbed, shell=True)

            #### Convert shifted bed to bam
            shiftedbedtobam = [
                'bedToBam -i', output_file,
                '-g', fai_file, '-bed12', '>', output_file.split('.bed')[0] + '.bam'
            ]

            shiftedbedtobam = ' '.join(shiftedbedtobam)
            sp.check_call(shiftedbedtobam,shell=True)

            command = [
                'samtools sort -@', self.threads, '-T', outdir + '/', output_file.split('.bed')[0] + '.bam',
                '-o', output_file.split('.bed')[0] + '.bam'
            ]

            command = ' '.join(command)
            sp.check_call(command, shell=True)

            #### Make coverage tracks for plus and minus strands
            strand = ['+', '-']
            bedgraph_extension = ['_plus.bedgraph', '_minus.bedgraph']
            bigwig_extension = ['_plus.bw', '_minus.bw']

            x = 0

            for j in strand:

                command = [
                    'bedtools genomecov -bg -strand', j, '-5 -i', output_file, '-g', fai_file, #'-scale', str(scale_fac),
                    '>', output_file.split('.bed')[0] + bedgraph_extension[x]
                ]

                command = ' '.join(command)
                sp.check_call(command, shell=True)

                command = [
                    'LC_COLLATE=C sort -k1,1 -k2,2n', output_file.split('.bed')[0] + bedgraph_extension[x],
                    '>', output_file.split('.bed')[0] + '_sorted' + bedgraph_extension[x]
                ]

                command = ' '.join(command)
                sp.check_call(command, shell=True)

                command = [
                    'bedGraphToBigWig', output_file.split('.bed')[0] + '_sorted' + bedgraph_extension[x], fai_file,
                    output_file.split('.bed')[0] + bigwig_extension[x]
                ]

                command = ' '.join(command)
                sp.check_call(command, shell=True)

                x = x + 1

        print(ctw.CBEIGE + ctw.CBOLD + 'Identification of cross-linked nucleotides is completed!!!' + ctw.CEND + '\n')

    def cldata_extractor(self):

        outdir = os.path.join(self.home_dir, 'crosslink_data')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        outdir_summary = os.path.join(self.home_dir, 'summary_files')
        if not os.path.isdir(outdir_summary): os.mkdir(outdir_summary)

        plus = sorted(glob.glob(outdir + '/' + '*UV*_plus.bedgraph'))
        minus = sorted(glob.glob(outdir + '/' + '*UV*_minus.bedgraph'))

        ctw = ColorTextWriter.ColorTextWriter()

        cl_summary = pd.DataFrame({'Summary': ['Total cross-link events', 'Number of stacked cross-link events', 'Number of nucleotides with stacked cross-link events']})

        count = 0
        frames = []
        count_list = []
        column_headers = []

        for (i, j) in zip(plus, minus):

            print(ctw.CRED + 'Cross-link data extraction from ' + ctw.CBLUE + os.path.basename(i) + ctw.CRED + ' and ' + ctw.CBLUE + os.path.basename(j) + ctw.CRED + ' ...' + ctw.CEND + '\n')

            column_headers.append(os.path.basename(i).split('_shifted')[0])

            command = [
                'cat', i, j, '|',
                "awk 'BEGIN{ totalcount=0 }{ totalcount += (($3 - $2) * $4) }END{ print totalcount }'"
            ]

            command = ' '.join(command)
            results = sp.check_output(command, shell=True, universal_newlines=True)
            count_list.append(int(results))

            command = [
                'cat', i, j, '|',
                "awk 'BEGIN{ totalstackedcount=0 }{ if($4 > 1) totalstackedcount += (( $3 - $2) * $4) }END{ print totalstackedcount }'"
            ]

            command = ' '.join(command)

            results = sp.check_output(command, shell=True, universal_newlines=True)
            count_list.append(int(results))

            command = [
                'cat', i, j, '|',
                "awk 'BEGIN{ totalstackedpos=0 }{ if($4 > 1) totalstackedpos += ($3 - $2) }END{ print totalstackedpos }'"
            ]

            command = ' '.join(command)

            results = sp.check_output(command, shell=True, universal_newlines=True)
            count_list.append(int(results))

            count_details = pd.DataFrame({i: count_list})

            if count == 0:
                frames = [cl_summary, count_details]
                count = 1
                count_list = []

            else:
                frames.append(count_details)
                count_list = []

        df = pd.concat(frames, axis=1)
        summary_column_header = 'Summary'
        final_column_headers = [summary_column_header] + column_headers

        df.columns = [final_column_headers]
        df.to_csv(outdir_summary + '/' + 'Cross-link_Data_Summary.csv', index=False)
        return df

    def ins_del_finder(self):

        outdir_summary = os.path.join(self.home_dir, 'summary_files')
        if not os.path.isdir(outdir_summary): os.mkdir(outdir_summary)

        bam_list = sorted(glob.glob(self.input_dir + '*UV*.bam'))

        ctw = ColorTextWriter.ColorTextWriter()

        summary = pd.DataFrame({'Summary': ['Mapped with Insertions', 'Mapped with Deletions']})

        count = 0
        frames = []
        count_list = []
        column_headers = []

        for i in bam_list:

            print(ctw.CRED + 'Detection of Insertions/Deletions: ' + ctw.CBLUE + os.path.basename(i) + ctw.CRED + ' ...' + ctw.CEND + '\n')

            column_headers.append(os.path.basename(i).split(self.extensions[4])[0])

            #### Convert from bam to sam format

            output_file = i.split(self.extensions[4])[0] + self.extensions[1]

            command = [
                'samtools view', i, '-o', output_file
            ]

            command = ' '.join(command)
            sp.check_call(command, shell=True)

            #### Total number of reads mapped with Insertions

            command = [
                'cut -f 6', output_file, '| grep I | wc -l'
            ]

            command = ' '.join(command)
            results = sp.check_output(command, shell=True, universal_newlines=True)
            count_list.append(int(results))

            #### Total number of reads mapped with Deletions

            command = [
                'cut -f 6', output_file, '| grep D | wc -l'
            ]

            command = ' '.join(command)
            results = sp.check_output(command, shell=True, universal_newlines=True)
            count_list.append(int(results))

            count_details = pd.DataFrame({i: count_list})

            if count == 0:
                frames = [summary, count_details]
                count = 1
                count_list = []

            else:
                frames.append(count_details)
                count_list = []

        df = pd.concat(frames, axis=1)
        summary_column_header = 'Summary'
        final_column_headers = [summary_column_header] + column_headers

        df.columns = [final_column_headers]
        df.to_csv(outdir_summary + '/' + 'Insertion_Deletion_Data_Summary.csv', index=False)
        return df