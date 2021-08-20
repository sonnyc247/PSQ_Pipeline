### working from/adapting code from AIBS 

#.libPaths(c("//allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Olivia/3.5/"))

library(Rsamtools)
library(GenomicAlignments)
library(edgeR)
library(rtracklayer)

args<- commandArgs(TRUE) 
suppressMessages(library(Rsamtools)) 
suppressMessages(library(GenomicAlignments)) 
suppressMessages(library(edgeR))
suppressMessages(library(rtracklayer)) 

message(library(Rsamtools))
message(library(GenomicAlignments))
message(library(edgeR))
message(library(rtracklayer))

fn="/external/rprshnas01/kcni/ychen/RandomTestDir/AIBS_InExPipeTest/AIBS_STAR/Aligned.sortedByCoord.out.bam"
fn="/external/rprshnas01/netdata_kcni/stlab/RamseyMielnik/STAR_results/coord_bams/343_1Aligned.sortedByCoord.out.bam"
fn = "/external/rprshnas01/netdata_kcni/stlab/kristina_patchseq_processed/STAR_results/sortedbams/100_sequenceAligned.sortedByCoord.out.bam"
expc="testcase"
species="mouse"
path="/external/rprshnas01/kcni/ychen/git/Vgatcollab/"

#Load counting reference
load("/external/rprshnas01/kcni/ychen/git/Vgatcollab/reference.count.Rdata")

mito <- "NC_005089.1"

# params<-ScanBamParam(tag = "NH", what=c("qname","flag")) ##modified to the below

total.range <- GRanges("2", IRanges(25291177, 25319253))

params<-ScanBamParam(tag = "NH", what=c("qname","flag"))
params<-ScanBamParam(what=scanBamWhat())
params<-ScanBamParam(which=total.range, tag = "NH", what=c("qname","flag"))

readGAlignments(fn,param=params)

#Error handling if BAM has no reads
reads =  tryCatch(readGAlignments(fn,param=params),
                  error=function(error_message){message("There are no genome reads mapped in BAM") 
                    reads=c() })

reads = readGAlignmentPairs(fn,param=params)

#No ercc counts file is generated if there are no reads.
if(!is.null(reads)){
  
  total = length(unique(values(first(reads))$qname))
  reads.unique = reads[values(first(reads))$flag < 1024]
  
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
    mito.reads =subsetByOverlaps(reads, mito.exon.list, ignore.strand=TRUE)
    mito.counts = length(unique(values(first(mito.reads))$qname))
    
    tRNA.reads = subsetByOverlaps(reads, tRNA.exon.list, ignore.strand=TRUE)
    tRNA.counts = length(unique(values(first(tRNA.reads))$qname))
    
    rRNA_rmsk.reads = subsetByOverlaps(reads, rRNA.exon.list, ignore.strand=TRUE)
    rRNA_rmsk.counts = length(unique(values(first(rRNA_rmsk.reads))$qname))
    
    
    ncRNA.exon.list =  exon.list[which(anno$type_from_gtf=="ncRNA")]
    ncRNA.reads = subsetByOverlaps(reads, ncRNA.exon.list, ignore.strand=TRUE)
    ncRNA.counts = length(unique(values(first(ncRNA.reads))$qname))
    
    rRNA_gtf.exon.list =  exon.list[which(anno$type_from_gtf=="rRNA")]
    rRNA.reads = subsetByOverlaps(reads, rRNA_gtf.exon.list, ignore.strand=TRUE)
    rRNA.counts = length(unique(values(first(rRNA.reads))$qname))
    
  }else{ #Gorilla (gtf is ensembl)
    
    mito.counts=0
    tRNA.counts=0
    ncRNA.exon.list =  exon.list[which(anno$type_from_gtf=="3prime_overlapping_ncRNA")]
    ncRNA.reads = subsetByOverlaps(reads, ncRNA.exon.list, ignore.strand=TRUE)
    ncRNA.counts = length(unique(values(first(ncRNA.reads))$qname))
    
    rRNA_rmsk.reads = subsetByOverlaps(reads, rRNA.exon.list, ignore.strand=TRUE)
    rRNA_rmsk.counts = length(unique(values(first(rRNA_rmsk.reads))$qname))
    
    rRNA_gtf.exon.list =  exon.list[which(anno$type_from_gtf=="rRNA")]
    rRNA.reads = subsetByOverlaps(reads, rRNA_gtf.exon.list, ignore.strand=TRUE)
    rRNA.counts = length(unique(values(first(rRNA.reads))$qname))
    
    
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

## OUR TRY at "No ercc counts file is generated if there are no reads".
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
    mito.reads =subsetByOverlaps(reads, mito.exon.list, ignore.strand=TRUE)
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

genecounttestdf <- PCRtestdf[PCRtestdf$exon_count==0,]
sum(genecounttestdf$intron_count > 0)

sum((df$exon_count + df$intron_count) > 0 )

write.csv(df,paste0(path,"/",expc,"_results.csv"),quote=F,row.names=F)



#output QC
qc=NULL
qc<-c(paste("reads_aligned_to_exons",sum(exon.count,na.rm=T)),
      paste("reads_aligned_to_introns",sum(intron.count,na.rm=T)),
      paste("reads_aligned_to_rrna_rmsk", rRNA_rmsk.counts),
      paste("reads_aligned_to_rrna",rRNA.counts),
      paste("reads_aligned_to_trna", tRNA.counts),
      paste("reads_aligned_to_ncrna",ncRNA.counts),
      paste("reads_aligned_to_mt_exons",mito.counts),
      paste("reads_aligned_to_intergenic",length(reads.unique)-(sum(exon.count,na.rm=T)+sum(intron.count,na.rm=T))),       #changed to unique reads
      paste("reads_aligned_unique", length(reads.unique))
)

write.table(qc,paste0(path,"/",expc,"_qctmp.txt"),sep="\t",quote=F,row.names=F,col.names=F)
