# Isoform psaudo-alignment and quantification with Kallisto

## <ins> Overview </ins>

This document serves as documentation and pseudo-tutorial for processing RNAseq fastq files to obtain gene expression measurements broken down by annotated isoforms of genes. 

Gene isoforms are forms of mature mRNAs that contain specifc exon combinations. For example, Gene A might have 3 known exons, and two established isoforms - one that contains exons 1 and 3 and one that contains exons 2 and 3. 

Some of the contents in this document is specific to the CAMH high-performance computing cluster (SCC), though most of the information should be generally applicable. Contents - code, settings to use, etc. are currently being tested/confirmed, but a basic version of the workflow has been worked out.

### *References*

Kallisto homepage - https://pachterlab.github.io/kallisto/

Kallisto publication - https://www.nature.com/articles/nbt.3519

Kallisto manual - https://pachterlab.github.io/kallisto/manual

Kallisto information for starting (main reference for the code in our current in-development workflow) - https://pachterlab.github.io/kallisto/starting

## <ins> Input </ins>

### *RNAseq reads (necessary)*

* These are stored essentially in text files, one or two for each sample. 
* The most common format is .fastq, but could be in another format, such as .txt. 
* These files are often compressed, with common compression extensions being ".gz" and ".tar". Fastq.gz files do not need to be decompressed for use in Kallisto. 
* RNAseq experiments can be single or double stranded - meaning 1 or 2 files exist per sample, respectively. Double stranded files may end with "_1" and "_2" after a specimen name and before file extensions. 

Example filenames of RNAseq reads:

Samplename1.txt for a sample1

Samplename2_1.fastq.gz and Samplename2_2.fastq.gz for a sample2

### *Reference transcriptome (necessary)*

* These are similar to genome reference assemblies (sizeable, strictly formatted text files), but contain sequences of mRNA only.
* We have used .fna files with Kallisto without issue - such as "GCF_000001635.23_GRCm38.p3_rna.fna.gz" from NCBI.
* Any reference (Refseq, Ensembl, etc.) can be used, but it would be good to consider downstream analyses when picking which to use. For example, if you want to combine or compare multiple RNAseq datasets, it may be good to use the same reference.
* Similar to above, make sure to note/keep track of which reference is going to be used for reproducibility and future use cases.

### *Metadata (preferred/as needed)*

Most important for using publicly available datasets. Helps you understand what the RNAseq reads refer to/which samples the reads come from, and also can help identify whether the reads come from a single-stranded or double-stranded RNAseq experiment.

## <ins> Running Kallisto </ins>

### *First time - setting up the transcriptome reference*

If working with Kallisto for the first time on a given computer, follow installation directions in https://pachterlab.github.io/kallisto/download .

If working on the CAMH SCC, Kallisto is already installed with Pyhon 3.8.5 and can simply be loaded in a linux terminal (such as via MobaXterm):

```bash {cmd}

module load lang/Python/3.8.5-Anaconda3-2021.03

```

Especially when working on the SCC, remember to note down the version of Kallisto and Python you use for future reproducibility and reporting.

For the first time aligning and quantifying reads with Kallisto, we need to set up the reference for use.

```bash {cmd}

kallisto index \ #running Kallisto for transcriptome reference index processing
         -i "/path_to_where_you_want_to_store_processed_reference/some_reference_name.idx" \ #where we want to store the processed transcriptome as an "idx" file [1]
         "/path_to_reference_file_like_a_fasta/" \ #the mRNA reference file
  
```

Note: note/write down what the exact command was for processing the mRNA reference for later use; there is no specific output of this command that helps track it down later

### *Basic run for one sample*

After setting up the transcriptome reference, we can align RNAseq read files to it and get our quantified isoform reads

```bash {cmd}

kallisto quant \ #running Kallisto for pseudo-alignment and quantification
         -i "/path_to_stored_processed_reference/some_reference_name.idx" \ #same index path/file as above [1]
         -o "/path_to_where_you_want_to_store_Kallisto_output/" \ #where we want to store the output of Kallisto
         -b 100 \ #number of times/samples to bootstrap, 100 is the example used in the "getting started" page of the Kallisto website
         "~/Data/Samplename2_1.fastq.gz" \ #first of a pair of fastq files for a sample
         "~/Data/Samplename2_2.fastq.gz" #second of a pair of fastq files for a sample

```

#### Other options of interest

In addition to the basic parameters above, we can also try:

* "--single" to specify single-end reads/fastq files
* "-l [integer]" length of fragments (necessary for processing single-end reads) 
* "-s [integer]" standard deviation of fragment lengths (necessary for processing single-end reads)

#### Outputs

Outputs of the Kallisto run are saved in the "/path_to_where_you_want_to_store_Kallisto_output/" specified with the run command. Two files in particular are useful to note:

* "abundance.tsv" contains the result proper, it gives an estimate for every isoform in the reference used for pseudo-alignment/quantification for the sample being processed
* "run_info.json" contains details of the specific run, including the bash command used, which you can refer back later to if needed

### *Processing multiple files*

* See https://github.com/sonnyc247/PSQ_Pipeline/blob/master/Code/Alignment/Isoform-Quantifiers/Kallisto/kallistobatch.sh for an example script of how we currently process multiple samples with Kallisto.
* The strategy/script is based on https://combine-lab.github.io/salmon/getting_started/, where isoform quantifications were performed for multiple files in a loop with Salmon. 
* In essence, put all the fastq files in a directory, and loop through all the files to process them with Kallisto.
* We will switch away from this in the future as it currently takes ~20-30 minutes to process one single-cell RNAseq sample.
 
* Afterwards, we can use something like https://github.com/sonnyc247/PSQ_Pipeline/blob/master/Code/Alignment/Isoform-Quantifiers/Kallisto/combine_csv_isoform.R to go through all output files and combine results into one count-matrix.
