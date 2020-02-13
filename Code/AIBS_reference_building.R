### code adapted from the Allen Institute for Brain Science, courtesy of Olivia Fong and Zizhen Yao
#this is for building a reference object for quantifying introns and exons in of ____-seq reads

library(GenomicFeatures)
library(rtracklayer)

# setwd("//allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Ensembl_Release93")
# GTFs
# reference <- import.gff("Mus_musculus.GRCm38.93.gff3")

reference <- import.gff("/external/rprshnas01/kcni/ychen/References/Mus_musculus.GRCm38.98.gtf")
keeplevels=levels(seqnames(reference))[1:22]


reference=keepSeqlevels(reference, keeplevels, pruning.mode="coarse")

mchr=read.table("mouse_chr.txt",as.is=T,check.names=F) ##skipped
rownames(mchr)=mchr$V1 ##skipped
seqlevels(reference)=mchr[seqlevels(reference),"V7"] ##skipped


#get exons from gtf
reference.exons = reference[reference$type %in% "exon"]

#split exons according to transcript, get start/end, seqname and strand info, to get 3' UTR 
exon.list = split(reference.exons, reference.exons$transcript_id)
ts.start=sapply(start(exon.list), min)
ts.end=sapply(end(exon.list), max)
ts.sn = sapply(seqnames(exon.list), function(x)runValue(x)[1])
ts.strand = sapply(strand(exon.list), function(x)runValue(x)[1])
reference.transcript = GRanges(IRanges(start=ts.start, end = ts.end), seqnames=ts.sn,strand=ts.strand)
reference.transcript$transcript_id = names(ts.start)
pos.transcript = reference.transcript[strand(reference.transcript)=="+"]
pos.transcript$tag = paste(pos.transcript$transcript_id, end(pos.transcript))
neg.transcript = reference.transcript[strand(reference.transcript)=="-"]
neg.transcript$tag = paste(neg.transcript$transcript_id, start(neg.transcript))
reference.exons$end.tag = paste(reference.exons$transcript_id, end(reference.exons))
reference.exons$start.tag = paste(reference.exons$transcript_id, start(reference.exons))
pos.3utr = reference.exons$end.tag %in% pos.transcript$tag
neg.3utr = reference.exons$start.tag %in% neg.transcript$tag


###Extend 3'UTR by 100bp
flank = 100
end(reference.exons[pos.3utr]) = end(reference.exons[pos.3utr]) + flank
start(reference.exons[neg.3utr]) = start(reference.exons[neg.3utr]) - flank



#get genes
#reference.mRNA = reference[reference$type %in% "mRNA"] #changed for gorilla
#gene.start = tapply(start(reference.mRNA), reference.mRNA$gene_name, min)  ##maybe check later
#gene.end = tapply(end(reference.mRNA), reference.mRNA$gene_name, max) ##maybe check later

reference.gene = reference[reference$type %in% "gene"]



#get unique exons
unique.exons = unique(reference.exons[,c("gene_id", "gene_name","transcript_id")])
exons=unique.exons


#get unique gene symbols
gene.name = setNames(reference$gene_name, reference$gene_id)

tmp = !duplicated(gene.name)
gene.name=gene.name[tmp]

#split by gene symbol
exon.list = split(reference.exons, reference.exons$gene_name)
exon.list=GenomicRanges::reduce(exon.list)
exon.list=as(exon.list,"GRangesList") #Add this line to make it compatible to LIMS, until LIMS is updated to new GenomicRanges.  New GenomicRanges package make CompressedGRangesList by default which LIMS library does not recongize ##performed the first time

save(exon.list, file="reference.exons.rda")


#Get introns
tmp1=reference.gene
strand(tmp1) = "*"
tmp2=exons
strand(tmp2)="*"
###find all introns
introns = setdiff(tmp1,tmp2)


###Assign introns to genes
tmp=as.matrix(findOverlaps(introns, tmp1, type="within", select="all"))

intron.list=split(introns[tmp[,1]],tmp1$gene_name[tmp[,2]])
intron.list=as(intron.list,"GRangesList")
save(intron.list, file="reference.introns.rda")
exonic.gene.sizes <- unlist(lapply(exon.list,function(x){sum(width(x))}))
save(exonic.gene.sizes,file="exonic.gene.sizes")

idx=match(names(exon.list),gene.name)
ind=match(names(exon.list),reference.gene$gene_name)


anno=data.frame(gene_name=names(exon.list),gene_id=names(gene.name)[idx],typeofgene=reference.gene$gene_biotype[ind],type_from_gtf=reference.gene$gene_biotype[ind],mt=FALSE,stringsAsFactors=F)
gene_id=anno$gene_id


mtgenes=unique(reference$gene_name[which(seqnames(reference)=="NC_005089.1")])
mtidx=match(mtgenes,names(exon.list))
anno$mt[mtidx]=TRUE




rRNA.exon.list=reference[reference$gene_biotype %in% "rRNA"]
tRNA.exon.list=reference[reference$gene_biotype %in% "tRNA"]
save(exon.list,intron.list,exonic.gene.sizes,rRNA.exon.list,tRNA.exon.list,anno,gene_id,file="reference.count.Rdata")




