# User guide
Please read through this step-by-step guide to setup and begin analysis of your data (**very IMPORTANT!!!**). This is a simple ***six*** step process:

### Step 1: Clone the repository
To clone the current repository on to your home directory using terminal, first navigate to the ***home directory*** (i.e., where you want analyzed data to be deposited), paste and enter the following command:
```
git clone https://github.com/jkkbuddika/eCLIP-Data-Analyzer.git
ls
```
> eCLIP-Data-Analyzer       

Once the cloning is completed:
```
mv eCLIP_Data_Analyzer/*/ ./
rm -rf eCLIP_Data_Analyzer
ls
```
> add_mat   
> scripts             
> environment               

### Step 2: Setup the miniconda environment
There are two ways to set up the conda environment: (1) Using the *environment.yml* file in the **environment** directory or (2) manually creating a conda environment and installing all required packages. The advantage of using the first approach is, it gives the conda environment I used when I was writing this python pipeline. However, the second approach let you install the latest versions of required packages. Note that based on the version of a similar tool, the output results can be varied. Let's go through how to both of these.

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
conda install -n dataanalyzer -c conda-forge -c bioconda fastqc cutadapt star qualimap samtools deeptools subread multiqc pandas singularity umi_tools pureclip bamtools ucsc-bedgraphtobigwig ucsc-bigWigMerge
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

### Step 3: Download additional materials
#### Download the singularity image from GoSTRIPES workflow to use TagDust2 for rRNA removal
To do this:
```
cd add_mat
singularity pull --name gostripes.simg shub://BrendelGroup/GoSTRIPES
cd ..
```
Note that the *add_mat* directory by default contain rRNA sequences for *D. melanogaster*, *C. elegans* and *S. cerevisiae* downloaded from Ensembl [BioMart](https://www.ensembl.org/biomart/martview/b1eec568acae1f43251215e8bd8f26fd).
```
cd add_mat
ls
cd ..
```
> Celegans_rRNA.txt	  
> Dro_rRNA.txt		    
> gostripes.simg        
> Yeast_rRNA.txt    

Here is the link to [GoSTRIPES](https://github.com/BrendelGroup/GoSTRIPES) workflow hosted by the [Brendel Group](http://brendelgroup.org/). 

### Step 4: Analysis mode selection and defining additional experiment specific variables
To specify experiment specific variables, open and update "GeneralVariables.py" module using emacs text editor.
```
cd scripts
emacs GeneralVariables.py
```
- **seq_method** : Run-mode. Options are 'single' (single-end data analysis) or 'paired' (paired-end data analysis).
- **rRNA_list** : Name of the rRNA sequence list in *add_mat* directory. Ex: 'Dro_rRNA.txt'
- **genome** : Biomart link to the genome of interest. Ex: Link to the Drosophila genome is 'ftp://ftp.ensembl.org/pub/release-99/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna_sm.toplevel.fa.gz'
- **feature** : Biomart link to the genome annotation of interest. Ex: Link to the Drosophila genome annotation is 'ftp://ftp.ensembl.org/pub/release-99/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.99.gtf.gz'
- **stranded** : Strandedness of the library preparation. Options are **0** (unstranded), **1** (stranded) or **2** (reversely stranded)
- **diff_features** : Features to be counted using featureCounts. To check supported features, download and open the annotation file using "less" command.

Once necessary changes are being made:
```
Ctrl+x+s then Ctrl+x+c ## To save and quit emacs
cd ..
ls
```
> add_mat  
> scripts                     
> environment

### Step 5: Input data preparation
The pipeline uses input files in .fastq format for analysis. To upload input data, navigate first to the home directory and create a directory *raw_sequences*.
```
mkdir raw_sequences
ls
```
> add_mat  
> raw_sequences   
> scripts                       
> environment

Then upload input sequences to the *raw_sequences* directory. Naming of files is ***very important*** and follow the recommended naming scheme. Name of an input fastq file must follow the following order:
***'_sample_adaptor_R1.fastq'*** or/and ***'_sample_adaptor_R2.fastq'***
> Note that ***'sample'*** above supports *three* options: (1) **IN** for ***input*** datasets (which usually is a single dataset), (2) **UV1**, **UV2**, ..., **UVn** for any number of ***UV samples***, and (3) **nonUV** for ***non-UV controls*** (which usually is a one dataset).                   

> Note that the ***'adaptor'*** above denotes the 3'-adaptor used in that particular library. Options are limited to conventional eCLIP adaptors: **L19**, **X1A**, **X1B**, **X2A**, **X2B**, **A01** or **B06**. Make sure to use the adaptor identities given here.        

> Here is a pair of acceptable input file names: ***pum2_IN_L19_R1.fastq***, ***pum2_IN_L19_R2.fastq***         

Remember, this naming scheme is vitally important!!! Pay extra attention to this.

You can use simple bash commands like below to quickly automate the renaming for you. Let's assume that the name of the pair of reads you have is *pum2_In_S12_R1_001.fastq* and *pum2_In_S12_R2_001.fastq*. Following two bash commands should convert those file names to the correct format:

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

### Step 6: Executing the pipeline
All executables of the pipeline are written onto *run.py* module. To start analyzing data activate the conda environment above, navigate to the scripts directory and execute *run.py* using python.
```
source activate dataanalyzer
cd scripts
python run.py
```

This will start running the pipeline! Time to get some rest. Watch a movie or a couple of episodes of a TV show! You got plenty of time!!!

After the run, I am sure you are eager to track the number of reads in each step. Read the following section which simplifies looking into this. Bash commands are helpful in this matter.

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

