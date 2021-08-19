# Isoform psaudo-alignment and quantification with Salmon

## <ins> Overview </ins>

This document serves as documentation and pseudo-tutorial for processing RNAseq fastq files to obtain gene expression measurements broken down by annotated isoforms of genes. 

Gene isoforms are forms of mature mRNAs that contain specifc exon combinations. For example, Gene A might have 3 known exons, and two established isoforms - one that contains exons 1 and 3 and one that contains exons 2 and 3. 

Some of the contents in this document is specific to the CAMH high-performance computing cluster (SCC), though most of the information should be generally applicable. Contents - code, settings to use, etc. are currently being tested/confirmed, but a basic version of the workflow has been worked out.

### *References*

Salmon homepage - https://combine-lab.github.io/salmon/

Salmon publication - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/

Salmon manual/documentation - https://salmon.readthedocs.io/en/latest/salmon.html

Salmon information for starting (main reference for the code in our current in-development workflow) - https://combine-lab.github.io/salmon/getting_started/

## <ins> Input </ins>

### *RNAseq reads (necessary)*

* These are stored essentially in text files, one or two for each sample. 
* The most common format is .fastq, but could be in another format, such as .txt. 
* These files are often compressed, with common compression extensions being ".gz" and ".tar". Fastq.gz files do not need to be decompressed for use in Salmon. 
* RNAseq experiments can be single or double stranded - meaning 1 or 2 files exist per sample, respectively. Double stranded files may end with "_1" and "_2" after a specimen name and before file extensions. 

Example filenames of RNAseq reads:

Samplename1.txt for a sample1

Samplename2_1.fastq.gz and Samplename2_2.fastq.gz for a sample2

### *Reference transcriptome (necessary)*

* These are similar to genome reference assemblies (sizeable, strictly formatted text files), but contain sequences of mRNA only.
* We have used .fna files with Salmon without issue - such as "GCF_000001635.23_GRCm38.p3_rna.fna.gz" from NCBI.
* Any reference (Refseq, Ensembl, etc.) can be used, but it would be good to consider downstream analyses when picking which to use. For example, if you want to combine or compare multiple RNAseq datasets, it may be good to use the same reference.
* Similar to above, make sure to note/keep track of which reference is going to be used for reproducibility and future use cases.

### *Metadata (preferred/as needed)*

Most important for using publicly available datasets. Helps you understand what the RNAseq reads refer to/which samples the reads come from, and also can help identify whether the reads come from a single-stranded or double-stranded RNAseq experiment.

## <ins> Running Salmon </ins>

### *First time - setting up the transcriptome reference*

If working with Salmon for the first time on a given computer, follow installation directions in https://combine-lab.github.io/salmon/getting_started/#obtaining-salmon .

If working on the CAMH SCC, Salmon is already installed with Pyhon 3.8.5 and can simply be loaded in a linux terminal (such as via MobaXterm):

```bash {cmd}

module load lang/Python/3.8.5-Anaconda3-2021.03

```

Especially when working on the SCC, remember to note down the version of Salmon and Python you use for future reproducibility and reporting.

For the first time aligning and quantifying reads with Salmon, we need to set up the reference for use.

```bash {cmd}

salmon index \ #running Salmon for reference processing
       -t "/path_to_reference_file_like_a_fasta/" \ #the mRNA reference file, can be compressed
       -i "/path_to_where_you_want_to_store_processed_reference/" \ #where we want to store the processed transcriptome reference [1]
  
```

Note: note/write down what the exact command was for processing the mRNA reference for later use; there is no specific output of this command that helps track it down later, though the output directory has other files with details of the command's run, such as Salmon version in "versionInfo.json".

### *Basic run for one sample*

After setting up the transcriptome reference, we can align RNAseq read files to it and get our quantified isoform reads

```bash {cmd}

salmon quant \ #running Salmon for pseudo-alignment and quantification
       -i "/path_to_stored_processed_reference/" \ #same processed transcriptome path as above [1]
       -l A \ #automatically infer library type
       -1 "~/Data/Samplename2_1.fastq.gz" \ #first of a pair of fastq files for a sample
       -2 "~/Data/Samplename2_2.fastq.gz" \ #second of a pair of fastq files for a sample 
       -p 8 \ #number of threads to use
       --validateMappings \ #use more sensitive/accurate mapping algorithm
       -o "/path_to_where_you_want_to_store_Salmon_output/" #where we want to store the output of Salmon

```

#### Outputs

Outputs of the Salmon run are saved in the "/path_to_where_you_want_to_store_Salmon_output/" specified with the run command. 

* "quant.sf" contains the result proper, it gives an estimate for every isoform in the reference used for pseudo-alignment/quantification for the sample being processed
* "cmd_info.json" contains details of the specific run, but does not contain the exact command entered for the run
* "salmon_quant.log" in the "logs" folder contains printed output/feedback from during the run

### *Processing multiple files*

* See https://github.com/sonnyc247/PSQ_Pipeline/blob/master/Code/Alignment/Isoform-Quantifiers/Salmon/salmonbatch.sh for an example script of how we currently process multiple samples with Salmon.
* The strategy/script is based on https://combine-lab.github.io/salmon/getting_started/ . 
* In essence, put all the fastq files in a directory, and loop through all the files to process them with Salmon.
* This runs quite quickly. With the way we run it, Salmon currently takes ~1 minute or less to process one single-cell RNAseq sample. However, there may be reason to use bootstrapping with Salmon in the future, which would increase runtime and warrant higher priority of parallelization of the runs instead.
 
* Afterwards, we can use something like https://github.com/sonnyc247/PSQ_Pipeline/blob/master/Code/Alignment/Isoform-Quantifiers/Salmon/combine_csv_isoform.R to go through all output files and combine results into one count-matrix.
