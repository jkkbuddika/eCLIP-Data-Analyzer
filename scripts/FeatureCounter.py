import os
import glob
import pandas as pd
import subprocess as sp
import ColorTextWriter

class FeatureCounter:

    def __init__(self, home_dir, input_dir, gfeature, stranded, feature_dir, feature_file, extensions, seq_method):
        self.home_dir = home_dir
        self.input_dir = input_dir
        self.gfeature = gfeature
        self.stranded = stranded
        self.feature_dir = feature_dir
        self.feature_file = feature_file
        self.extensions = extensions
        self.seq_method = seq_method

    def feature(self):

        file_list = sorted(glob.glob(self.input_dir + '*.bam'))

        ctw = ColorTextWriter.ColorTextWriter()

        print(ctw.CBEIGE + ctw.CBOLD + 'Feature counting started ...' + ctw.CEND + '\n')

        i = 0
        while i < len(self.gfeature):

            outdir = os.path.join(self.home_dir, 'featureCounts' + '_' + self.gfeature[i])
            if not os.path.isdir(outdir): os.mkdir(outdir)

            print('\n' + 'Quantifying ' + self.gfeature[i] + 's ...' + '\n')

            command = [
                'featureCounts -t', self.gfeature[i],
                '-F GTF -g gene_name -O -M -s', self.stranded
            ]

            if self.seq_method == 'paired': command.extend(['-p -B -C'])

            command.extend([
                '-a', self.feature_dir + self.feature_file.split('.gz')[0],
                '-o', outdir + '/' + self.gfeature[i] + self.extensions[3],
                ' '.join([self.input_dir + '{0}'.format(j.split(self.input_dir)[1]) for j in file_list])
            ])

            command = ' '.join(command)
            sp.check_call(command, shell=True)

            ### Manipulating the FeatureCounts output
            data = pd.read_csv(outdir + '/' + self.gfeature[i] + self.extensions[3], sep='\t', header=0, index_col=0, skiprows=1)

            data = data.drop(data.iloc[:, 0:5], axis=1)

            data.to_csv(outdir + '/' + self.gfeature[i] + '_DESeq2_Input' + self.extensions[3], sep='\t')

            i = i + 1

        if i == len(self.gfeature):
            print('\n' + ctw.CBEIGE + ctw.CBOLD + 'Feature counting done!!!' + ctw.CEND + '\n')

