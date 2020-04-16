import os
import glob
import pandas as pd

lines = ['total input reads', 'successfully extracted', 'extracted', 'problems with architecture']

class TDSummaryProcessor:

    def __init__(self, home_dir, input_dir):
        self.home_dir = home_dir
        self.input_dir = input_dir

    def data_mining(self, summary_line, line_dis):
        self.summary_line = summary_line
        self.line_dis = line_dis
        data = self.summary_line.split(']')[1].split(self.line_dis)[0]
        data = data.strip()
        return data

    def td_summary(self):

        outdir = os.path.join(self.home_dir, 'summary_files')
        if not os.path.isdir(outdir): os.mkdir(outdir)

        file_list = sorted(glob.glob(self.input_dir + '*.txt'))

        count = 0
        frames = []
        column_headers = []

        for i in file_list:

            column_headers.append(os.path.basename(i).split('_tagdustout')[0])

            with open(i) as f:
                content = f.readlines()
                id = 0
                summary = ''
                summary_list = []
                td_summary = pd.DataFrame({'Summary': [lines[0], lines[1], lines[2], lines[3]]})

                for line in content:
                    if lines[id] in line:
                        data = self.data_mining(line, lines[id])
                        summary_list.append(data)
                        summary += data + '\n'
                        id = id + 1

                    if id == 4:
                        break

                td_details = pd.DataFrame({i: summary_list})

                if count == 0:
                    frames = [td_summary, td_details]
                    count = 1

                else:
                    frames.append(td_details)

        df = pd.concat(frames, axis=1)
        summary_column_header = 'Summary'
        final_column_headers = [summary_column_header] + column_headers

        df.columns = [final_column_headers]
        df.to_csv(outdir + '/' + 'TagDust_summary.csv', index=False)
        return df
