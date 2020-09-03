import os
import subprocess as sp
import TarExtractor
import ColorTextWriter

class WebDownloader:

    def __init__(self, home_dir, output_dir, input_path, download_file):
        self.home_dir = home_dir
        self.output_dir = output_dir
        self.input_path = input_path
        self.download_file = download_file

    def download(self):

        outdir = os.path.join(self.home_dir, self.output_dir)
        if not os.path.isdir(outdir): os.mkdir(outdir)

        ctw = ColorTextWriter.ColorTextWriter()

        if len(os.listdir(outdir + '/')) == 0:

            print(ctw.CRED + 'Downloading ' + ctw.CBLUE + self.download_file + ctw.CRED + ' ...' + ctw.CEND + '\n')

            sp.call('wget ' + self.input_path + self.download_file, shell=True)

            print(ctw.CBLUE + self.download_file + ctw.CRED + ' Downloaded!!!' + ctw.CEND)

            dir_path = os.path.dirname(os.path.realpath(__file__))
            te = TarExtractor.TarExtractor(dir_path + '/' + self.download_file, self.output_dir)
            te.tar_extractor()

            if self.download_file.endswith('.tar.gz'):
                sp.call('rm -r ' + self.download_file, shell=True)

            else:
                new_file = self.download_file.split('.gz')[0]
                sp.call(['mv', new_file, outdir + '/' + new_file])

        else:
            print(ctw.CRED + "Destination directory contain files. Ignore download!!!" + ctw.CEND + '\n')