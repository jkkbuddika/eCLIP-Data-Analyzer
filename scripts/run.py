import UserDefinedOptions
import CommonOptions
import FastQCRunner
import Tagduster
import TDSummaryProcessor
import CutAdapt
import UMITools
import WebDownloader
import RefGenMaker
import GenomeAligner
import SamTools
import CrosslinkHunter
import BigWigFileMaker
import FeatureCounter
import MultiQCRunner
import PureCLIPPeakCaller
import ColorTextWriter

#########################
# Executing the Program #
#########################

ctw = ColorTextWriter.ColorTextWriter()
print('\n' + ctw.CRED + ctw.CBOLD + 'Initiating eCLIP data analyzer ...' + ctw.CEND + '\n')

trimming = input('Have the adaptor sequences are trimmed (Answer Yes or No)? ')

gv = UserDefinedOptions.UserDefinedOptions()
cv = CommonOptions.CommonOptions()

qc_1 = FastQCRunner.FastQCRunner(cv.home_dir, cv.fastqc_raw, cv.raw_sequences_dir)
qc_1.fastqc()

if trimming == 'Yes':
    ex_umi = UMITools.UMITools(cv.home_dir, cv.raw_sequences_dir, cv.extensions, cv.umi_seq, cv.seq_method)
    ex_umi.extract_UMI()

elif trimming == 'No':
    ca = CutAdapt.CutAdapt(cv.home_dir, cv.raw_sequences_dir, cv.r2_adapter, cv.r2_adapter_seq, cv.extensions, cv.seq_method)
    ca.cutadapt()

    ex_umi = UMITools.UMITools(cv.home_dir, cv.cutadapt_dir, cv.extensions, cv.umi_seq, cv.seq_method)
    ex_umi.extract_UMI()

td = Tagduster.Tagduster(cv.home_dir, cv.Threads, cv.tagdust_singu, cv.umi_extract, cv.rRNA_path, cv.extensions, cv.seq_method)
td.tagdust()

tdsp = TDSummaryProcessor.TDSummaryProcessor(cv.home_dir, cv.tagdust_out)
tdsp.td_summary()

wd = WebDownloader.WebDownloader(cv.home_dir, cv.genome_dir_name, cv.genome_path, cv.genome_file)
wd.download()

wd = WebDownloader.WebDownloader(cv.home_dir, cv.feature_dir_name, cv.feature_path, cv.feature_file)
wd.download()

rg = RefGenMaker.RefGenMaker(cv.home_dir, cv.Threads, cv.genome_fa, gv.species, cv.genes_gtf)
rg.refgen()

sa = GenomeAligner.GenomeAligner(cv.home_dir, cv.tagdust_out, cv.Threads, cv.ref_genome, cv.extensions, cv.genes_gtf, cv.seq_method)
sa.aligner()

ss = SamTools.SamTools(cv.home_dir, cv.star_aligned, cv.Threads, cv.extensions, cv.seq_method)
ss.sam_filtering()

dedup_umi = UMITools.UMITools(cv.home_dir, cv.sam_sorted, cv.extensions, cv.umi_seq, cv.seq_method)
dedup_umi.dedup()

cl = CrosslinkHunter.CrosslinkHunter(cv.home_dir, cv.umi_dedup, cv.genome_fa, cv.Threads, cv.extensions)
cl.crosslink()
cl.cldata_extractor()
cl.ins_del_finder()

bw = BigWigFileMaker.BigWigFileMaker(cv.home_dir, cv.umi_dedup, cv.extensions)
bw.bigwig()

bw = BigWigFileMaker.BigWigFileMaker(cv.home_dir, cv.crosslink_data, cv.extensions)
bw.bigwig()

fc = FeatureCounter.FeatureCounter(cv.home_dir, cv.umi_dedup, cv.diff_features, cv.feature_dir, cv.feature_file, cv.extensions, cv.seq_method)
fc.feature()

mqc = MultiQCRunner.MultiQCRunner(cv.home_dir)
mqc.multiqc()

if cv.seq_method == 'paired':
    sr2r = SamTools.SamTools(cv.home_dir, cv.umi_dedup, cv.Threads, cv.extensions, cv.seq_method)
    sr2r.sam_retrieval()

    pc = PureCLIPPeakCaller.PureCLIPPeakCaller(cv.home_dir, cv.sam_retrieved_r2, cv.genome_fa, gv.pcMode, cv. extensions)
    pc.peak_caller()

elif cv.seq_method == 'single':
    pc = PureCLIPPeakCaller.PureCLIPPeakCaller(cv.home_dir, cv.umi_dedup, cv.genome_fa, gv.pcMode, cv.extensions)
    pc.peak_caller()

bw = BigWigFileMaker.BigWigFileMaker(cv.home_dir, cv.pure_clip_pc, cv.extensions)
bw.bigwig()

ctw = ColorTextWriter.ColorTextWriter()
print('\n' + ctw.CGREEN + ctw.CBOLD + ctw.CBLINK + ctw.CURL + 'Data analysis successfully completed !!!' + ctw.CBLACK + ctw.CEND + '\n')
