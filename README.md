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
* RNAseq experiments can be single or double stranded. Double stranded files may end with "_1" and "_2" after a specimen name and before extensions. A single STAR command with two read/fastq files as input is how you specify to the program that you have paired reads.

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

Note: if you ever need to track down what your exact command was for generating the genome reference for star, you can find it in "genomeParameters.txt" in the genomeDir

#### Basic alignment

After setting up the genome, we can align RNAseq read files to it

```bash {cmd}
# For general/basic STAR alignment

STAR --genomeDir "/path_to_where_you_stored_processed_genome/" \ #same genomeDir as above [1]
     --sjdbGTFfile "/path_to_reference_file_likely_a_gtf/" \ #reference gtf file again, most important if you want STAR to quantify gene reads
     --readFilesIn "/path_to_ONE_RNAseq_reads_file_likely_a_fastq/" #the RNAseq reads (use only one file per command if your files are single-stranded, only an appropriate pair of files if your data come as paired reads)

# For example:

STAR --genomeDir "~/Genomic_references/Ensembl/Human/Release_103/USE_THIS_genomeDir/" \ 
     --sjdbGTFfile "~/Genomic_references/Ensembl/Human/Release_103/Raw/Homo_sapiens.GRCh38.103.gtf" \ #reference gtf file, was also used to generate "/USE_THIS_genomeDir/"
     --readFilesIn "~/Data/Samplename2_1.fastq.gz" "~/Data/Samplename2_2.fastq.gz" #NOT two samples, this is a pair of fastq files for ONE sample

```

In addition to the basic parameters above, we have also used previously:

Important
=========
* "--readFilesCommand zcat" to read in compressed files like .gz
* "--quantMode ["TranscriptomeSAM" and (separated by one space)/or "GeneCounts"] for a transcript-coordinate bam used in RSEM quantification and/or STAR to output its own gene counts, respectively
* "--outSAMtype ["BAM SortedByCoordinate" or "None"]" for sorted genomic-coordinates of reads as a bam file or to not output a sam or bam (but still output a TranscrptomeSAM specified above, for example), respectively
* "--outFilterMultimapNmax 1" to only output uniquely mapped reads
* "--runRNGseed [number]" to fix RNG for primary assignment of multimapped reads, which involves RNG (may be important for reproducibility)

Other/convenience
=================

* "--runThreadN [number]" to specify number of threads
* "--outFileNamePrefix ["/path_of_desired_output/output_prefix"]" to specify output dir AND add a prefix to the output files if desired
* "--outReadsUnmapped Fastx" to output unmapped reads into another file (that can be processed/analyzed further)

Putting the above together, we can use the following basic command to run a STAR alignment destined for RSEM quantification afterwards:

```bash {cmd}
# For STAR alignment into RSEM

STAR --genomeDir "/path_to_where_you_stored_processed_genome/" \ #same genomeDir as above [1]
     --sjdbGTFfile "/path_to_reference_file_likely_a_gtf/" \ #reference gtf file again, most important if you want STAR to quantify gene reads
     --readFilesIn "/path_to_ONE_RNAseq_reads_file_likely_a_fastq/" \ #the RNAseq reads (one sample; use only one file or pair of files for paired-reads)
     --quantMode TranscriptomeSAM \ #generate a transcript-coordinate bam for RSEM
     --outSAMtype None # we only need the BAM above; to save space/unless otherwise desired, we can choose not to output a genome-coordinate sam or bam
     
```


## Acknowledgements

Special thanks to Justin Chee (https://github.com/cheejus2) and Jordan Sicherman (https://github.com/jsicherman) for help/consultation in developing this pipeline.

This work was supported in part by funding provided by Brain Canada, in partnership with Health Canada, for the Canadian Open Neuroscience Platform initiative.

