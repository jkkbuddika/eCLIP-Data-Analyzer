import os
import tarfile
import subprocess as sp
import ColorTextWriter

class TarExtractor:

    def __init__(self, zip_file, output_dir):
        self.zip_file = zip_file
        self.output_dir = output_dir

    def tar_extractor(self):

        ctw = ColorTextWriter.ColorTextWriter()

        print('\n' + ctw.CRED + 'Unzipping ' + ctw.CBLUE + os.path.basename(self.zip_file) + ctw.CRED +' ...' + ctw.CEND + '\n')

        self.zip_formats = ['.tar', '.tar.gz']

        if self.zip_file.endswith(tuple(self.zip_formats)):
            tar = tarfile.open(self.zip_file)
            tar.extractall(self.output_dir)
            tar.close()
        else:
            sp.call(['gunzip', self.zip_file])