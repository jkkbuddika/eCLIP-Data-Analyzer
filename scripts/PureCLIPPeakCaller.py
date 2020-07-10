import os
import glob
import subprocess as sp
import ColorTextWriter
import BigWigFileMaker

class PureCLIPPeakCaller():

    def __init__(self, home_dir, input_dir, genome_fa, extensions):
        self.home_dir = home_dir
        self.input_dir = input_dir
        self.genome_fa = genome_fa
        self.extensions = extensions

    def peak_caller(self):

        outdir = os.path.join(self.home_dir, 'pure_clip_pc')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        ctw = ColorTextWriter.ColorTextWriter()

        uv1_list = sorted(glob.glob(self.input_dir + '*UV1*.bam'))
        uv2_list = sorted(glob.glob(self.input_dir + '*UV2*.bam'))
        input_list = sorted(glob.glob(self.input_dir + '*IN*.bam'))

        x = 0

        for (i, j) in zip(uv1_list, uv2_list):

            print(ctw.CRED + 'Merging replicates and indexing:' + ctw.CBLUE + os.path.basename(i) + ctw.CRED + ' and ' + ctw.CBLUE + os.path.basename(j) + ctw.CRED + ' ...' + ctw.CEND + '\n')

            merged_bam = outdir + '/' + os.path.basename(i).split('_UV')[0] + '_UV_aligned_sorted_dupRm_merged.bam'

            command = [
                'samtools merge', merged_bam, i, j
            ]

            command = ' '.join(command)
            sp.check_call(command, shell=True)
            sp.check_call(' '.join(['samtools index', merged_bam]), shell=True)

            #### Peak Calling with Input Normalization

            sp.check_call(' '.join(['samtools index', input_list[x]]), shell=True)

            outfile_prefix = outdir + '/' + os.path.basename(i).split(self.extensions[4])[0] + '_pureclip_crosslink'

            command = [
                'pureclip',
                '-i', merged_bam, '-bai', merged_bam + '.bai',
                '-g', self.genome_fa, '-ld -nt 8',
                '-o', outfile_prefix + '_INnorm_sites.bed',
                '-or', outfile_prefix + '_INnorm_regions.bed',
                '-ibam', input_list[x], '-ibai', input_list[x] + '.bai'
            ]

            command = ' '.join(command)
            sp.check_call(command, shell=True)

            command = [
                'cat', outfile_prefix + '_INnorm_sites.bed', '|',
                'cut -f 1,2,3,4,5,6 >', outfile_prefix + '_INnorm_sites_short.bed'
            ]

            command = ' '.join(command)
            sp.check_call(command, shell=True)

            x = x + 1

        #### Generate wiggle files for merged bam files
        bw = BigWigFileMaker.BigWigFileMaker(self.home_dir, outdir + '/', self.extensions)
        bw.bigwig()

        print(ctw.CBEIGE + ctw.CBOLD + 'Peak Calling Completed!!!' + ctw.CEND)
