# Tripathy lab RNAseq alignment and QC pipeline

This document serves as documentation and pseudo-tutorial for RNAseq (including bulk, sc/snRNAseq, and patch-seq) data processing and quality control methods in the Tripathy lab at the Krembil Centre for Neuroinformatics in the Centre for Addiction and Mental Health (CAMH). Some of the contents in this document is specific to the CAMH high-performance computing cluster (SCC), though most of the information should be generally applicable.

## Alignment

### Overview

This is the process of aligning raw/unmapped RNA fragments ("reads") to a reference genome. We use STAR for this.

### References

STAR repo - https://github.com/alexdobin/STAR
STAR manual - https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
STAR article - https://pubmed.ncbi.nlm.nih.gov/23104886/

### Input

#### RNAseq reads (necessary)

* These are stored essentially in text files, one for each sample. 
* The most common format is .fastq, but could be in another format, such as .txt. We have been able to run STAR with .txt files previously without any change to code/command.
* These files are often compressed, with common compression extensions being ".gz" and ".tar". Compressed files can be read with STAR directly by parameter specification in the command (see below). 
* RNAseq experiments can be single or double stranded. Double stranded files may end with "_1" and "_2" after a specimen name and before extensions. Inputting a pair of fastq files into star can also be specified in STAR parameters (see below).

Example filenames of RNAseq reads:

Samplename1.txt for a sample1
Samplename2_1.fastq.gz and Samplename2_2.fastq.gz for a sample2

#### Reference genome (necessary)

* These are genome reference assemblies (may be considered as large, strictly formatted text files too).
* Often in .fasta format. Could also be .fa or .fna (others possible as well). We have been able to use .fa and .fna without any change to STAR commands.
* Any genome reference can be used, but it would be good to consider downstream analyses when picking the genome to use. For example, if you want to combine or compare multiple RNAseq datasets, it may be good to use the same genome reference + release version.
* Similar to above, make sure to note/keep track of which genome reference + release is going to be used for reproducibility and future use cases.
* Most of our experience to date (Feb 2021) has been to use Ensembl references. There have been no issues.

#### Reference genome annotation (preferred/recommended)

* These are files (if available) that annotate their corresponding genome references, such as indicating which part of a genome is an exonic region for a given gene. These are also large, strictly formatted text files.
* These files are most important for quantifying and assigning reads to genes, but can also help with the alignment process.
* Common file extensions are .gtf and .gff. We have been able to use both without issue and without changes to the command in STAR.

#### Metadata (preferred/as needed)

Most important for using publicly available datasets. Helps you understand what the RNAseq reads refer to/which samples the reads come from, and also helps identify whether the reads come from a single-stranded or double-stranded RNAseq experiment.

### Running STAR

#### First time

If working with STAR the first time on a given computer, follow installation directions in https://github.com/alexdobin/STAR .

If working with CAMH SCC, STAR is already installed and can simply be loaded in linux terminal (such as via MobaXterm):

```bash {cmd}
module load STAR
```

For the first time aligning reads to any given genome, we need to set up the genome for STAR to use.

```bash {cmd}
STAR --runMode genomeGenerate \ #tells STAR that we are running a command for setting up the genome for later use
     --genomeDir "/path_to_where_you_want_to_store_processed_genome/" \ #this folder will be referenced when aligning reads later [1]
     --genomeFastaFiles "/path_to_genome_file_likely_a_fasta/" \ #the genome file
     --sjdbGTFfile "/path_to_reference_file_likely_a_gtf/" \ #the genome annotation file    
```

In addition to the basic parameters above, we have also used previously:
* "--runThreadN" to specify number of threads
* "--sjdbGTFtagExonParentTranscript Parent" for accomodating an NCBI/refseq genome

#### Basic alignment

After setting up the genome, we can align RNAseq read files to it

```bash {cmd}
# For single-stranded files

STAR --genomeDir "/path_to_where_you_stored_processed_genome/" \ #same genomeDir as above [1]
     --readFilesIn "/path_to_RNAseq_reads_file_likely_a_fastq/" \ #the RNAseq reads
     --outSAMtype BAM SortedByCoordinate \ #we want sorted BAM files as output 
```

In addition to the basic parameters above, we have also used previously:
* "--runThreadN" to specify number of threads
* "--outFilterMultimapNmax 1" to only output uniquely mapped reads

## Acknowledgements

Special thanks to Justin Chee (https://github.com/cheejus2) and Jordan Sicherman (https://github.com/jsicherman) for help/consultation in developing this pipeline.

This work was supported in part by funding provided by Brain Canada, in partnership with Health Canada, for the Canadian Open Neuroscience Platform initiative.

