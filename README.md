# eCLIP Data Analyzer
This repository contains a series of python based modules to automate sequencing data generated from eCLIP (enhanced UV crosslinking and immunoprecipitation) experiments. I ***highly recommend*** reading through this step-by-step manual *carefully* before you start analyzing your data. If you require a protocol for eCLIP library preparation please look at the Yeo Lab eCLIP library preparation protocol given in their original [paper](https://www.nature.com/articles/nmeth.3810).

## Requirements
The eCLIP data analyzer requires following tools to be installed for data analysis.

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) : Quality assessment
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) : Adaptor trimming
- [BBMap](https://github.com/BioInfoTools/BBMap) : To collapse exact PCR duplicates
- [UMI-tools](https://github.com/CGATOxford/UMI-tools) : Manipulation of UMIs (Unique Molecular Identifiers)
- [TagDust2](http://tagdust.sourceforge.net/) : Remove rRNA reads
- [STAR](https://github.com/alexdobin/STAR) : Mapping reads to a given genome
- [QualiMap](http://qualimap.bioinfo.cipf.es/) : Mapping quality assessment
- [SAMtools](https://github.com/samtools/samtools) : Sorting and indexing of mapped reads
- [deepTools](https://github.com/deeptools/deepTools/) : Generate bigwig files for IGV visualization
- [Subread](http://subread.sourceforge.net/) : Count features
- [MultiQC](https://github.com/ewels/MultiQC) : Summarize logs
- [PureCLIP](https://github.com/skrakau/PureCLIP) : CLIP-peak calling

We thank developers of these valueble tools!

## Analysis process
All analyzed data will be saved onto the home directory where you deposited the *scripts* directory. The pipeline first use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to assess the quality of raw input files. Conventional 3'-eCLIP adaptors (*i.e.*, based on the adaptor identity embedded in the file name) are removed from input sequence files using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and exact PCR duplicates are collapsed using clumpify.sh available through [BBMap](https://github.com/BioInfoTools/BBMap). Next, the pipeline utilize [UMI-tools](https://github.com/CGATOxford/UMI-tools) to extract 5' 10 nucleotide UMIs of R2 reads. These extracted UMIs are used for deduplication later. Subsequently, [TagDust2](http://tagdust.sourceforge.net/) is being used to remove rRNA contaminants from individual datasets. A datamining module in the pipeline will summarize TagDust2 rRNA removal logs and deposit mined data onto a file named **TagDust_summary.csv** which will be saved onto a directory named *summary_files*. Subsequently, user defined genome and annotation files are downloaded and a reference genome is generated using [STAR](https://github.com/alexdobin/STAR). Next, rRNA-depleted sequences are mapped to the given genome using [STAR](https://github.com/alexdobin/STAR) genome aligner. The pipeline use [QualiMap](http://qualimap.bioinfo.cipf.es/) to assess the quality of sequence alignment. The analysis scheme use [SAMtools](https://github.com/samtools/samtools) to coordinate sort, remove unmapped reads and index alignment output files. Sorted bam files along with index are then imported to [UMI-tools](https://github.com/CGATOxford/UMI-tools) to remove potential PCR duplicates. Then deduplicated bam files are used to generate (1) bigwig files for [IGV](https://software.broadinstitute.org/software/igv/) visualization using [deepTools](https://github.com/deeptools/deepTools/) and (2) count tables for user defined features using subread package [featureCounts](http://subread.sourceforge.net/). Furthermore, the pipeline integrates [MultiQC](https://github.com/ewels/MultiQC) to generate summary files in an interactive manner.

Cross-link sites can be recovered using R2 reads of eCLIP data. The pipeline allows analysis in ***two different*** modes: (1) Single-end mode where only *R2 reads* of the pair are used for analysis and (2) Paired-end mode where both R1/R2 reads are used at earlier steps of analysis, but extract only *mapped R2 reads* following deduplication. Subsequently, the pipeline use [PureCLIP](https://github.com/skrakau/PureCLIP) for peak calling. Output files of all these steps will be saved to your home directory. 

Now that you know the general outline of the analysis process, go through the step-by-step guide given [here](https://github.com/jkkbuddika/eCLIP-Data-Analyzer/blob/master/USERGUIDE.md) to analyze your eCLIP data.
