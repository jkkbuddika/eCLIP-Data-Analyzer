# User guide
Please read through this step-by-step guide to setup and begin analysis of your data (**very IMPORTANT!!!**). This is a simple ***FOUR*** step process:

### Step 1: Clone the repository
To clone the current repository on to your home directory using terminal, first navigate to the ***home directory*** (i.e., where you want analyzed data to be deposited), paste and enter the following command:
```
git clone https://github.com/jkkbuddika/eCLIP-Data-Analyzer.git
ls
```
> eCLIP-Data-Analyzer       

Once the cloning is completed:
```
mv eCLIP-Data-Analyzer/*/ ./
rm -rf eCLIP-Data-Analyzer
ls
```
> add_mat   
> scripts             
> environment               

### Step 2: Setup the miniconda environment
There are two ways to set up the conda environment: (1) Using the *environment.yml* file in the **environment** directory or (2) manually creating a conda environment and installing all required packages. The advantage of using the first approach is, it gives the conda environment I used when I was writing this python pipeline. However, the second approach let you install the latest versions of required packages. Note that based on the version of a similar tool, the output results can be varied. Let's go through how to both of these.

#### Install miniconda
Before creating the environment, [install](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda) miniconda and add conda to your PATH variable (see this [post](https://developers.google.com/earth-engine/python_install-conda) to learn more). Then update conda by running ```conda update conda```.

#### Setting up the miniconda environment with the *environment.yml* file
Simply copy, paste and run the following command on your terminal window.

```
cd environment
conda env create -f environment.yml
cd ..
```
Running this command will create a conda environment named ***dataanalyzer*** on your local computer.

#### Setting up the miniconda environment by manually installing required packages and dependencies
In this scenario, to setup the conda environment (i.e., dataanalyzer), run following terminal commands.

```
conda create -n dataanalyzer -c conda-forge -c bioconda python=3.7
conda install -n dataanalyzer -c conda-forge -c bioconda fastqc cutadapt star qualimap samtools deeptools subread multiqc pandas umi_tools pureclip bamtools ucsc-bedgraphtobigwig ucsc-bigWigMerge
```
#### Activate and deactivate the miniconda environment
To activate the enironment:
```
source activate dataanalyzer
```
To deactivate the environment run the following terminal command or simply close the terminal window:
```
source deactivate
```

For more details on managing conda enviroments [click here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#).

#### Install additional packages
In addition to above packages, eCLIP Data Analyzer require sucessful installation of [TagDust2](https://github.com/TimoLassmann/tagdust). Download TagDust2 from [here](https://sourceforge.net/projects/tagdust/) and install as described [here](http://tagdust.sourceforge.net/#install).

### Step 3: Input data preparation
The pipeline uses input files in .fastq format for analysis. To upload input data, navigate first to the home directory and create a directory *raw_sequences*.
```
mkdir raw_sequences
ls
```
> add_mat  
> environment       
> raw_sequences   
> scripts                       

Then upload input sequences to the *raw_sequences* directory. Naming of files is ***very important*** and follow the recommended naming scheme. Name of an input fastq file must follow the following order: ***'_sample_adaptor_R1.fastq'*** or/and ***'_sample_adaptor_R2.fastq'***    

> Note that ***'sample'*** above supports *three* options: (1) **IN** for ***input*** datasets (which usually is a single dataset), (2) **UV1**, **UV2**, ..., **UVn** for any number of ***UV samples***, and (3) **nonUV** for ***non-UV controls*** (which usually is a one dataset).                   

> Note that the ***'adaptor'*** above denotes the 3'-adaptor used in that particular library. Options are limited to conventional eCLIP adaptors: **L19**, **X1A**, **X1B**, **X2A**, **X2B**, **A01** or **B06**. Make sure to use the adaptor identities given here.        

> Here is a pair of acceptable input file names: ***pum2_IN_L19_R1.fastq***, ***pum2_IN_L19_R2.fastq***         

Remember, this naming scheme is vitally important!!! Pay extra attention to this.

You can use simple bash commands like below to quickly automate the renaming for you. Let's assume that the name of the pair of reads you have is *pum2_In_S12_R1_001.fastq* and *pum2_In_S12_R2_001.fastq*. Following two bash commands should convert those file names to the correct format: ***pum2_IN_L19_R1.fastq*** and ***pum2_IN_L19_R2.fastq***

```
for i in `ls *In*R1*`; do
newname="${i/%In_S12_R1_001.fastq/IN_L19_R1.fastq}"
mv -- "$i" "$newname"; 
done

for i in `ls *In*R2*`; do
newname="${i/%In_S12_R2_001.fastq/IN_L19_R2.fastq}"
mv -- "$i" "$newname"; 
done
```

### Step 4: Executing the pipeline
All executables of the pipeline are written onto *run.py* module. To start data analysis, activate the conda environment above, navigate to the scripts directory and execute *run.py* using python.
```
source activate dataanalyzer
cd scripts
python run.py
```
This should intiate running the analysis pipeline. Immediately, a couple questions will pop-up that you have to answer.
- **Enter the species code (Options: hs, mm, dm, ce, dr or custom):** Answer based on the species, **hs**: human, **mm**: mouse, **dm**: fruit fly, **ce**: *C. elegans*, **dr**: zebra fish or **custom**: any other model organism
- **Enter the Run Mode (Options: 0 for single-end or 1 for paired-end):** Pipeline supports analysis of both single-end or paired-end data. Answer **0** if single-end. Answer **1** if paired-end.
- **Enter the UMI Preference (Options: 0 for 5N or 1 for 10N):** Answer based on the used UMI-sequence length, **0**: NNNNN (5N), **1**: NNNNNNNNNN (10N).
- **Enter the PureCLIP Run Mode (Options: 0 for short-defined or 1 for long-undefined):** The pipeline allows running PureCLIP in two modes, **0**: short defined binding regions, **1**: larger binding regions.

If the species is **custom**, you have to:
- **FTP link to the genome to download:** Enter the link to the genome FASTA to download. For instance, if the custom species is yeast here is the Ensembl url to download the genome.
> ftp://ftp.ensembl.org/pub/release-100/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz               

- **FTP link to the annotation to download:** Enter the link to the corresponding GTF to download. For instance, if the custom species is yeast here is the Ensembl url to download the GTF.
> ftp://ftp.ensembl.org/pub/release-100/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.100.gtf.gz        

- Make sure you have transferred the custom rRNA sequence file to the *add_mat* directory. Note that the *add_mat* directory by default contain rRNA sequences for **hs, mm, dm, ce and dr**. Download rRNA sequences of the custom species from Ensembl [BioMart](https://www.ensembl.org/biomart/martview/), name the file **custom_rRNA.txt** and transfer into the *add_mat* directory.

You are all set!!! Let it run. Depending on the size of each file and the number of datasets, run time can vary so much!

## Retrieve additional information
1. If the directory of interest have a series of *.fastq* files, you can use the following bash command to get read counts saved into a *.txt* file in the same directory. As an example let's save read counts of the *raw_sequences* directory.
```
cd raw_sequences

for i in `ls *.fastq`; do
c=`cat $i | wc -l`
c=$((c/4))
echo $i $c
done > raw_readCounts.txt
```
> Executing the above bash command will save a file named *raw_readCounts.txt* in the *raw_sequences* directory with file name and number of reads in each file.

2. If the directory of interest have a series of *.bam* files, you can use the following bash command that uses [SAMtools](https://github.com/samtools/samtools). As an example let's save mapped read counts of the *star_aligned* directory.
```
cd star_aligned

for i in `ls *.bam`; do
echo ${i} $(samtools view -c $i)
done > bam_readCounts_aligned.txt
```
> Executing the above bash command will save a file named *bam_readCounts_aligned.txt* in the *star_aligned* directory with bam file names and number of reads that are mapped to the reference genome.
