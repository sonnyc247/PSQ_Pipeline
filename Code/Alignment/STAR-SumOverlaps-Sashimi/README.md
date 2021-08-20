# RNAseq alignment for Sashimi plots (such as with Gviz) and quantification (such as with GenomicAlignments)  

This document serves as documentation and pseudo-tutorial for RNAseq (including bulk, sc/snRNAseq, and patch-seq) data processing in the Tripathy lab at the Krembil Centre for Neuroinformatics in the Centre for Addiction and Mental Health (CAMH). 

This particular pipeline/workflow version generates coordinate-sorted BAM files that are then used in Gviz, GenomicAlignments, and other tools.

Some of the contents in this document is specific to the CAMH high-performance computing cluster (SCC), though most of the information should be generally applicable.

## <ins> Alignment </ins>

### *Overview*

This is the process of aligning raw/unmapped RNA sequence fragments ("reads") to a reference genome. We use STAR for this.

### *References*

STAR repo - https://github.com/alexdobin/STAR

STAR manual - https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

STAR article - https://pubmed.ncbi.nlm.nih.gov/23104886/

### *Input*

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
* Most of our experience to date (Aug 2021) has been to use Ensembl references, which have worked with no issues.
* We have tried Refseq references from NCBI, with some success only for the STAR step (and not for downstream commands with RSEM).

#### Reference genome annotation (preferred/recommended)

* These are files (if available) that annotate their corresponding genome references, such as indicating which part of a genome is an exonic region for a given gene. These are also large, strictly formatted text files.
* These files are most important for quantifying and assigning reads to genes, but can also help with the alignment process.
* Common file extensions are .gtf and .gff. We have been able to use both without issue and without changes to the command in STAR.

#### Metadata (preferred/as needed)

Most important for using publicly available datasets. Helps you understand what the RNAseq reads refer to/which samples the reads come from, and also helps identify whether the reads come from a single-stranded or double-stranded RNAseq experiment.

### *Running STAR*

#### First time

If working with STAR for the first time on a given computer, follow installation directions in https://github.com/alexdobin/STAR .

If working on the CAMH SCC, STAR is already installed and can simply be loaded in linux terminal (such as via MobaXterm):

```bash {cmd}

module load STAR

```

Especially when working on the SCC, remember to note down the version of STAR you use for future reproducibility and reporting.

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

#### Basic alignment to generate coordinate-sorted BAM files

The general command for aligning RNAseq read files to a genome reference and generating a coordinate-sorted BAM is:

```bash {cmd}

# For STAR alignment to generate coordinate-sorted BAM files

STAR --genomeDir "/path_to_where_you_stored_processed_genome/" \ #same genomeDir as above [1]
     --sjdbGTFfile "/path_to_reference_file_likely_a_gtf/" \ #reference gtf file again, most important if you want STAR to quantify gene reads
     --readFilesIn "/path_to_ONE_RNAseq_reads_file_likely_a_fastq/" \ #the RNAseq reads (one sample; use only one file or pair of files for paired-reads)
     --quantMode GeneCounts \ #likely default option, for quantifying gene counts downstream
     --outSAMtype BAM SortedByCoordinate #generate coordinate-sorted BAM for RSEM
     
```

In addition to the basic parameters above, we may also use other options:

<ins> Important </ins>

* "--readFilesCommand zcat" to read in compressed files like .gz
* "--outFilterMultimapNmax 1" to only output uniquely mapped reads
* "--runRNGseed [integer]" to fix RNG for primary assignment of multimapped reads, which involves RNG (may be important for reproducibility)

<ins> Other/convenience </ins>

* "--runThreadN [integer]" to specify number of threads
* "--outFileNamePrefix ["/path_of_desired_output/output_prefix"]" to specify output dir AND add a prefix to the output files if desired
* "--outReadsUnmapped Fastx" to output unmapped reads into another file (that can be processed/analyzed further)

#### Outputs

There are a few files of interest that STAR outputs into "/path_of_desired_output/":

* "[output_prefix]Aligned.sortedByCoord.out.bam" is our main output of interest for downstream processing and analyses - such as in Gviz and GenomicAlignments
* "[output_prefix]Log.out" for a record of the input command and terminal feedback messages/printouts
* "[output_prefix]Log.final.out" for some summartive quality metrics

## <ins> Gene quantification using coordinate-sorted BAM files and GenomicAlignments </ins>

### *Overview*

In this step, we quantify the number of reads for each gene in each sample, now that the reads have been aligned to the (annotated) genome. This workflow uses GenomicAlignments (R package)-based scripts from the Allen Institute for Brain Science (AIBS). 

### *References*

These scripts are based on 5 R packages:

* GenomicFeatures - https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003118 ; https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
* Rsamtools - https://bioconductor.org/packages/release/bioc/html/Rsamtools.html
* GenomicAlignments - https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html
* edgeR - https://academic.oup.com/bioinformatics/article/26/1/139/182458 ; https://bioconductor.org/packages/release/bioc/html/edgeR.html
* rtracklayer - https://academic.oup.com/bioinformatics/article/25/14/1841/225816?login=true ; https://bioconductor.org/packages/release/bioc/html/rtracklayer.html

### *Input*

#### Aligned RNAseq reads (necessary)

These are the coordinate-sorted BAM files generated above

#### Reference genome annotation (necessary)

* Same kind of reference genome annotation as used for STAR - a .gtf file, for example
* For a given sample, please use the exact same genome annotation reference for STAR alignment and gene quantification

### *Running quantification scripts*

#### First time - reference processing

If working for the first time with a given genome annotation on a given computer, run code from the following script:

https://github.com/sonnyc247/PSQ_Pipeline/blob/master/Code/Alignment/STAR-SumOverlaps-Sashimi/AIBS_reference_building.R

Use import.gff() on line 11 to import the reference annotation (gtf) file.

Then run all the rest of the code in the script before saving a reference Rdata file for future use with the next step of gene quantification.

#### Quantification of BAM files

We have modified a script from the AIBS for getting gene counts from BAM files (https://github.com/sonnyc247/PSQ_Pipeline/blob/master/Code/Alignment/STAR-SumOverlaps-Sashimi/adapted_AIBS_sample_bam_to_count_csv.R)

Modify the script on line 24, and load the appropriate reference Rdata file, such as one generated in the previous step

Afterwards, the script can be run simply with

```bash {cmd}

module load R #load R on the SCC
~/adapted_AIBS_sample_bam_to_count_csv.R "/path_to_directory_of_coord_sorted_bam_files/"

```

This will generate gene counts from BAM files in a folder, which can then be combined with a script like https://github.com/sonnyc247/PSQ_Pipeline/blob/master/Code/Alignment/STAR-SumOverlaps-Sashimi/combine_csv.R that operates on the output directory from "adapted_AIBS_sample_bam_to_count_csv.R", into one count matrix.

## <ins> Sashimi plots with coordinate-sorted BAM files </ins>

### *Overview*

In this step, we generate sashimiplots/read pile-ups for samples aligned to the (annotated) genome. 

### *References*

These scripts are based on X R packages:

* GenomicFeatures - https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003118 ; https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html
* Rsamtools - https://bioconductor.org/packages/release/bioc/html/Rsamtools.html
* GenomicAlignments - https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html
* edgeR - https://academic.oup.com/bioinformatics/article/26/1/139/182458 ; https://bioconductor.org/packages/release/bioc/html/edgeR.html
* rtracklayer - https://academic.oup.com/bioinformatics/article/25/14/1841/225816?login=true ; https://bioconductor.org/packages/release/bioc/html/rtracklayer.html

### *Input*

#### Aligned RNAseq reads (necessary)

These are the coordinate-sorted BAM files generated above

#### Reference genome annotation (necessary)

* Same kind of reference genome annotation as used for STAR - a .gtf file, for example
* For a given sample, please use the exact same genome annotation reference for STAR alignment and gene quantification

### *Running quantification scripts*
