#!/usr/bin/env Rscript
### contains a lot of code adapted from the Allen Institute for Brain Science, courtesy of Olivia Fong and Zizhen Yao
# script that produces individual, sample-level, count matrices from directory of bamfiles
# this is for MOUSE data
# arguments should be: 
# 1. directory of bamfiles 

# reading in and checking arguments from the command line
args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("One argument must be supplied (directory of bamfiles)", call.=FALSE)
} else if (length(args) > 1) {
  stop("One argument must be supplied (directory of bamfiles)", call.=FALSE)
} else if (!dir.exists(args[1])) {
  stop("Invalid directory", call.False)
}

library(Rsamtools)
library(GenomicAlignments)
library(edgeR)
library(rtracklayer)

load("/external/rprshnas01/kcni/ychen/References/Ens.GRCm38.98.Gtf.reference.count.Rdata") #loading sequence reference
mito <- "NC_005089.1" #set mito... maybe setting species? (unsure, from AIBS)
params<-ScanBamParam(tag = "NH", what=c("qname","flag")) #set parameters to retreive from bam files
dir.create(paste(args[1], "/InExCounts/", sep = "")) #create directory for output

# list of bam files in given directory
folder_files <- list.files(path = args[1], pattern = "\\.bam$")

for (i in folder_files) {
  
  fn = paste0(args[1], "/", i)
  
  #Error handling if BAM has no reads ##use readGAlignments for working with single stranded data; originally AIBS provided us paired-end data script using readGAlignmentPairs
  reads =  tryCatch(readGAlignments(fn,param=params),
                    error=function(error_message){message("There are no genome reads mapped in BAM") 
                      reads=c() }) 
  
  #No ercc counts file is generated if there are no reads. ##add first() to the next set of code, around "reads", within values() for dealing with paired-end data; you can do tests/checks comparing values(reads) and values(first(reads)) from single-stranded and paired data, respectively, to confirm (in terms of format/structure, etc.; if not specific data)
  if(!is.null(reads)){
    
    total = length(unique(values(reads)$qname))
    reads.unique = reads[values(reads)$flag < 1024]
    
    #counting introns first, then exons
    junc=junctions(reads.unique)
    no_junc = elementNROWS(junc)==0
    
    
    intron.olap=overlapsAny(reads.unique,intron.list, minoverlap=3, ignore.strand=TRUE)
    intron.count.tmp =summarizeOverlaps(intron.list,reads.unique[no_junc & intron.olap],"IntersectionNotEmpty",ignore.strand=TRUE)
    intron.count=matrix(0,length(exon.list),1)
    
    
    exon.count =assay(summarizeOverlaps(exon.list,reads.unique[!no_junc | !intron.olap],"IntersectionNotEmpty",ignore.strand=TRUE))
    rownames(intron.count)=rownames(exon.count)
    intron.count[rownames(intron.count.tmp),1]=assay(intron.count.tmp)
    colnames(exon.count)="exon_count"
    
    
    #calculate fpkm
    if (sum(exon.count!=0)){
      fpkms<-round(rpkm(exon.count,exonic.gene.sizes), digits=2)
    }else{
      fpkms=exon.count
    }
    colnames(fpkms)="fpkm"
    
    
    if (!is.null(mito)){ #Mouse or Human
      mito.exon.list =  exon.list[sum(seqnames(exon.list)==mito)] #mouse: NC_005089.1 human:NC_012920.1
      mito.reads = subsetByOverlaps(reads, mito.exon.list, ignore.strand=TRUE)
      mito.counts = length(unique(values(mito.reads)$qname))
      
      tRNA.reads = subsetByOverlaps(reads, tRNA.exon.list, ignore.strand=TRUE)
      tRNA.counts = length(unique(values(tRNA.reads)$qname))
      
      rRNA_rmsk.reads = subsetByOverlaps(reads, rRNA.exon.list, ignore.strand=TRUE)
      rRNA_rmsk.counts = length(unique(values(rRNA_rmsk.reads)$qname))
      
      
      ncRNA.exon.list =  exon.list[which(anno$type_from_gtf=="ncRNA")]
      ncRNA.reads = subsetByOverlaps(reads, ncRNA.exon.list, ignore.strand=TRUE)
      ncRNA.counts = length(unique(values(ncRNA.reads)$qname))
      
      rRNA_gtf.exon.list =  exon.list[which(anno$type_from_gtf=="rRNA")]
      rRNA.reads = subsetByOverlaps(reads, rRNA_gtf.exon.list, ignore.strand=TRUE)
      rRNA.counts = length(unique(values(rRNA.reads)$qname))
      
    }else{ #Gorilla (gtf is ensembl)
      
      mito.counts=0
      tRNA.counts=0
      ncRNA.exon.list =  exon.list[which(anno$type_from_gtf=="3prime_overlapping_ncRNA")]
      ncRNA.reads = subsetByOverlaps(reads, ncRNA.exon.list, ignore.strand=TRUE)
      ncRNA.counts = length(unique(values(ncRNA.reads)$qname))
      
      rRNA_rmsk.reads = subsetByOverlaps(reads, rRNA.exon.list, ignore.strand=TRUE)
      rRNA_rmsk.counts = length(unique(values(rRNA_rmsk.reads)$qname))
      
      rRNA_gtf.exon.list =  exon.list[which(anno$type_from_gtf=="rRNA")]
      rRNA.reads = subsetByOverlaps(reads, rRNA_gtf.exon.list, ignore.strand=TRUE)
      rRNA.counts = length(unique(values(rRNA.reads)$qname))
      
    }
    
  }else{
    exon.count=rep(0,length(exon.list))
    intron.count=rep(0,length(exon.list))
    fpkms=rep(0,length(exon.list))
    mito.counts=0
    tRNA.counts=0
    rRNA.counts=0
    ncRNA.counts=0
    rRNA_rmsk.counts=0
    reads=NULL
    reads.unique=NULL
  }
  
  
  df<- data.frame(gene_id=anno$gene_id,gene_symbol=names(exon.list),typeofgene=anno$type_from_gtf,mito=anno$mt,length=exonic.gene.sizes,intron_count=intron.count,exon_count=exon.count,fpkm=fpkms)
  
  # writing sample-level dataframe (counts of exons, introns, fpkm, etc.) as csv file
  write.csv(df, file = paste(args[1], "/InExCounts/", tools::file_path_sans_ext(i), "_count.csv", sep = ""))
  
}
