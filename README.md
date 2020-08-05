# eCLIP Data Analyzer
This repository contains a series of python based modules to automate the analysis of sequencing data generated from eCLIP (enhanced UV crosslinking and immunoprecipitation) experiments. I ***highly recommend*** reading through this step-by-step manual *carefully* before you start analyzing your data. If you require a protocol for eCLIP library preparation please look at the Yeo Lab eCLIP library preparation protocol given in their original [paper](https://www.nature.com/articles/nmeth.3810).

## Requirements
The eCLIP data analyzer requires following tools to be installed for data analysis (see [Step 2](https://github.com/jkkbuddika/eCLIP-Data-Analyzer/blob/master/USERGUIDE.md#step-2-setup-the-miniconda-environment) of the User Guide).

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) : Quality assessment
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) : Adaptor trimming
- [UMI-tools](https://github.com/CGATOxford/UMI-tools) : Extraction of UMIs and deduplication
- [TagDust2](http://tagdust.sourceforge.net/) : Removal of rRNA reads
- [STAR](https://github.com/alexdobin/STAR) : Mapping reads to a given genome
- [QualiMap](http://qualimap.bioinfo.cipf.es/) : Mapping quality assessment
- [SAMtools](https://github.com/samtools/samtools) : Sorting and indexing of mapped reads
- [deepTools](https://github.com/deeptools/deepTools/) : Generate bigwig files for IGV visualization
- [Subread](http://subread.sourceforge.net/) : Count features
- [MultiQC](https://github.com/ewels/MultiQC) : Summarize logs
- [BedTools](https://github.com/arq5x/bedtools2) : Extraction of cross-linked nucleotides
- [PureCLIP](https://github.com/skrakau/PureCLIP) : eCLIP-peak calling

We thank developers of these valueble tools!

## Analysis process
All analyzed data will be saved to subdirectories inside the home directory where you have deposited the *scripts* directory. The pipeline first use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to assess the quality of raw input files. Conventional 3'-eCLIP adaptors (*i.e.*, based on the adaptor identity embedded in the file name) are removed from input sequence files using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/). Next, the pipeline utilize [UMI-tools](https://github.com/CGATOxford/UMI-tools) to extract 5' 10 nucleotide UMIs of R2 reads. These extracted UMIs are used for deduplication later. Subsequently, [TagDust2](http://tagdust.sourceforge.net/) is being used to remove rRNA contaminants from individual datasets. A datamining module in the pipeline will summarize TagDust2 rRNA removal logs and deposit mined data onto a file named **TagDust_summary.csv** which will be saved onto a directory named *summary_files*. Subsequently, user defined genome and annotation files are downloaded and a reference genome is generated using [STAR](https://github.com/alexdobin/STAR). Then, rRNA-depleted sequences are mapped to the given genome using [STAR](https://github.com/alexdobin/STAR) genome aligner and only uniquely mapped reads are kept. The pipeline use [QualiMap](http://qualimap.bioinfo.cipf.es/) to assess the quality of sequence alignment. The analysis scheme use [SAMtools](https://github.com/samtools/samtools) to coordinate sort and index alignment output files. Sorted bam files along with indices are then imported to [UMI-tools](https://github.com/CGATOxford/UMI-tools) to remove potential PCR duplicates. Then deduplicated bam files are used to generate (1) bigwig files for [IGV](https://software.broadinstitute.org/software/igv/) visualization using [deepTools](https://github.com/deeptools/deepTools/), (2) count tables for user defined features using subread package [featureCounts](http://subread.sourceforge.net/), (3) data related to crosslinking events and (4) peak calling output files using [PureCLIP](https://github.com/skrakau/PureCLIP). Furthermore, the pipeline integrates [MultiQC](https://github.com/ewels/MultiQC) to generate summary files in an interactive manner. As mentioned above,output files of all these steps will be saved to subdirectories in your home directory.

Cross-link sites can be recovered using R2 reads of eCLIP data. Therefore, the pipeline allows analysis in ***two different*** modes: (1) ***Single-end mode*** where only *R2 reads* of the pair are used for analysis (we frequently use this mode for our analyzes) and (2) ***Paired-end mode*** where both R1/R2 reads are used at earlier steps of analysis, but extract only *mapped R2 reads* following deduplication.

Now that you know the general outline of the analysis process, go through the step-by-step guide given [here](https://github.com/jkkbuddika/eCLIP-Data-Analyzer/blob/master/USERGUIDE.md) to analyze your eCLIP data.

If you use eCLIP Data Analyzer, please cite: https://doi.org/10.1101/2020.06.27.175174
