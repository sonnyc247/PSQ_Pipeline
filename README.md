# Tripathy lab RNAseq data processing and QC pipeline

This repository holds documentation and pseudo-tutorials for RNAseq (including bulk, sc/snRNAseq, and patch-seq) data processing and quality control methods in the Tripathy lab at the Krembil Centre for Neuroinformatics in the Centre for Addiction and Mental Health (CAMH). Some of the contents in this document is specific to the CAMH high-performance computing cluster (SCC), though most of the information should be generally applicable.

## <ins> Data processing </ins>

### *Overview*

Data processing of RNAseq data generally involves transforming the data from random fragments of RNA sequences into human-readable information. Most commonly, we perform alignment of RNAseq reads to an annotated reference genome followed by quantification to produce count matrices. We can also produce plots of aligned RNAseq reads to better understand mRNA species and composition in samples.

Currently, we have four pipeline-tutorials at various stages of development:

1. STAR-RSEM for alignment and quantification (also recommended for beginners to RNAseq data processing): https://github.com/sonnyc247/PSQ_Pipeline/tree/master/Code/Alignment/STAR-RSEM
2. STAR-generation of coordinated-sorted BAM files for gene quantification and read-pile-up plots: https://github.com/sonnyc247/PSQ_Pipeline/tree/master/Code/Alignment/STAR-SumOverlaps-Sashimi
3. Kallisto for isoform-level pseudo-alignment and quantification: https://github.com/sonnyc247/PSQ_Pipeline/tree/master/Code/Alignment/Isoform-Quantifiers/Kallisto
4. Salmon for isoform-level pseudo-alignment and quantification: https://github.com/sonnyc247/PSQ_Pipeline/tree/master/Code/Alignment/Isoform-Quantifiers/Salmon

## <ins> Quality assessment and control </ins>

### *Overview*

It is difficult to know form processed data whether the RNAseq data from a particular experiment or sample is of high or low quality. Information from before, during, and after RNAseq data processing can be used in combination to assess the quality of RNAseq data from a particular sample.

Currently, this work and pipeline is in development. Some example code and tasks that we currently use to assess RNAseq quality are stored in https://github.com/sonnyc247/PSQ_Pipeline/tree/master/Code/Quality_Assessment.

## <ins> Acknowledgements </ins>

Special thanks to Justin Chee (https://github.com/cheejus2), Jordan Sicherman (https://github.com/jsicherman), and Derek Howard (https://github.com/derekhoward) for help and consultation in developing this pipeline.

This work was supported in part by funding provided by Brain Canada, in partnership with Health Canada, for the Canadian Open Neuroscience Platform initiative.

