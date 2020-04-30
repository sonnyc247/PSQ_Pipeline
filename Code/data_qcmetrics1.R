### this script is not yet generalized, taken from another project, continuing from Moffit_comparison.R in that project

### use metadata spreadsheet as a starting point

# do the below, if needed/skipping moffit/previous steps, etc.
library(readxl)
Phenodata <- read_excel("/external/rprshnas01/netdata_kcni/stlab/kristina_patchseq/Phenodata.xlsx")
Phenodata[21,1] <- "ps71" #fixing typo from collaborator (this is specific to Kristina's data)

Pheno_result_metrics <- as.data.frame(Phenodata) #initialize matrix for qc metrics results

### getting/importing/checking intron+exon data count matrices, working with PCR-duplicate data

combined_exons <- read.csv("/external/rprshnas01/netdata_kcni/stlab/kristina_patchseq_processed/STAR_results/coord_bams/InExCounts/gs_combined_exons.csv", row.names=1, stringsAsFactors=FALSE)
rownames(combined_exons) <- combined_exons$Gene_ID #since we want matrices, we want to move gene names to row names and remove the gene name column
combined_exons <- combined_exons[, 2:37]
colnames(combined_exons) <- paste0("ps",(gsub("[^0-9]", "", colnames(combined_exons)))) #fixing sample names

combined_introns <- read.csv("/external/rprshnas01/netdata_kcni/stlab/kristina_patchseq_processed/STAR_results/coord_bams/InExCounts/gs_combined_introns.csv", row.names=1, stringsAsFactors=FALSE)
rownames(combined_introns) <- combined_introns$Gene_ID #since we want matrices, we want to move gene names to row names and remove the gene name column
combined_introns <- combined_introns[, 2:37]
colnames(combined_introns) <- paste0("ps",(gsub("[^0-9]", "", colnames(combined_introns)))) #fixing sample names

InEx_Matrix <- combined_exons + combined_introns #generating combined matrix

### getting number of reads for exons, introns, and total reads

qc_metric_holder <- as.data.frame(colSums(combined_exons)) #sum of reads per sample (column)
Pheno_result_metrics <- merge(Pheno_result_metrics, qc_metric_holder, by.x = "sample", by.y = "row.names") #merge new metric to metric df
colnames(Pheno_result_metrics)[5] <- "Exonic_reads" #label the metric

qc_metric_holder <- as.data.frame(colSums(combined_introns)) #sum of reads per sample (column)
Pheno_result_metrics <- merge(Pheno_result_metrics, qc_metric_holder, by.x = "sample", by.y = "row.names") #merge new metric to metric df
colnames(Pheno_result_metrics)[6] <- "Intronic_reads" #label the metric

Pheno_result_metrics$Genic_reads <- Pheno_result_metrics$Exonic_reads + Pheno_result_metrics$Intronic_reads #getting total reads form adding introns and exons

### getting number of genes detected by Exons alone, Introns alone, and all detected genes (exons alone, introns alone, as well as genes detected with both introns and exons)

qc_metric_holder <- as.data.frame(colSums(InEx_Matrix != 0) - colSums(combined_introns != 0)) #getting genes that are only detected via exons (number of total detected genes minus number of genes with at least 1 intron)
Pheno_result_metrics <- merge(Pheno_result_metrics, qc_metric_holder, by.x = "sample", by.y = "row.names") #merge new metric to metric df
colnames(Pheno_result_metrics)[8] <- "Exon_only_genes" #label the metric

qc_metric_holder <- as.data.frame(colSums(InEx_Matrix != 0) - colSums(combined_exons != 0)) #getting genes that are only detected via introns (number of total detected genes minus number of genes with at least 1 exon)
Pheno_result_metrics <- merge(Pheno_result_metrics, qc_metric_holder, by.x = "sample", by.y = "row.names") #merge new metric to metric df
colnames(Pheno_result_metrics)[9] <- "Intron_only_genes" #label the metric

qc_metric_holder <- as.data.frame(abs(colSums(InEx_Matrix != 0))) #getting all detected genes (can be detected by exons only, introns only, or both)
Pheno_result_metrics <- merge(Pheno_result_metrics, qc_metric_holder, by.x = "sample", by.y = "row.names") #merge new metric to metric df
colnames(Pheno_result_metrics)[10] <- "All_detected_genes" #label the metric

### getting information from log files, based on https://stackoverflow.com/questions/51611127/grep-summary-statistics-from-log-files ; initially/originally used, not too important
# specifically, we are getting numbers of input/sequenced reads, multimapped reads, and unmapped reads

log_dir <- '/external/rprshnas01/netdata_kcni/stlab/kristina_patchseq_processed/STAR_results/star_logs/Init_QCMetrics/' #set directory of star qc log files
log_name_list = list.files(path = log_dir, recursive = F, full.names = F, pattern = "Log.final.out")  #get list of files

library(tidyverse) #maybe unneeded

qc_metric_holder <- data.frame(Sample=character(), Sequenced_reads=double(), Uniquely_aligned_reads=double(), Multi_mapped_reads=double(), Unmapped_reads=double(), stringsAsFactors=FALSE) #reinitialize the df into a blank holder df

for (i in 1:length(log_name_list)){
  
  sample_name = paste0("ps",(gsub("[^0-9]", "", log_name_list[i]))) #setting/cleaning sample name
  df <- read.delim(paste0(log_dir,log_name_list[i]), header=FALSE, comment.char="#", as.is = TRUE) #load log file
  qc_metric_holder[i,1] = sample_name #adding sample name to holder df
  qc_metric_holder[i,2] = as.numeric(df[5,2]) #adding sequenced/input reads to holder df
  qc_metric_holder[i,3] = as.numeric(df[8,2]) #adding uniquely aligned reads to holder df
  qc_metric_holder[i,4] = as.numeric(df[25,2]) #adding multi-mapped reads to holder df; this code works on our STAR code that only allows 1 alignment per read; will need to modify if some multimapped reads are allowed
  qc_metric_holder[i,5] = as.numeric(df[28,2]) + as.numeric(df[30,2]) + as.numeric(df[32,2]) #getting the sum of all types of unmapped reads (too many mismatches, too short, etc.) and adding the sums to the holder df
  print(c(sample_name, as.numeric(df[5,2]), as.numeric(df[8,2]), as.numeric(df[25,2]), qc_metric_holder[i,5])) #print out the values for a manual/visual sanity check; not necessary

}

Pheno_result_metrics <- merge(Pheno_result_metrics, qc_metric_holder, by.x = "sample", by.y = "Sample") #add new metrics to the metric df

### getting intergenic reads

Pheno_result_metrics$Intergenic_reads <- Pheno_result_metrics$Uniquely_aligned_reads - Pheno_result_metrics$Genic_reads #calculating intergenic reads and adding the metric

### getting mt genes/reads

#mtassess_df <- merge(InEx_Matrix, Gene_anno, by.x = 'row.names', by.y = 'Gene') #if needed, old/manual method of converting ensembl ids to mgi symbols; not needed if using allen code, which can output unique gene symbols instead of ensembl ids
mtassess_df <- InEx_Matrix[grep("mt-", row.names(InEx_Matrix)),] #isolating mitochondrial genes
mtassess_df <- mtassess_df[2:38,] #removing "Bhmt-ps1" row
qc_metric_holder <- as.data.frame(colSums(mtassess_df)) #getting sum of mt reads into qc df
Pheno_result_metrics <- merge(Pheno_result_metrics, qc_metric_holder, by.x = "sample", by.y = "row.names") #merge new metric to metric df
colnames(Pheno_result_metrics)[16] <- "mt_gene_reads" #label the metric

### getting some percentages (rather than just absolute read amounts that we were using previously)

Pheno_result_metrics$Unique_align_pct <- round(Pheno_result_metrics$Uniquely_aligned_reads/Pheno_result_metrics$Sequenced_reads*100, digits = 2) #percentage of uniquely aligned reads vs all input/sequenced reads
Pheno_result_metrics$Unmapped_pct <- round(Pheno_result_metrics$Unmapped_reads/Pheno_result_metrics$Sequenced_reads*100, digits = 2) #percentage of unmapped reads vs all input/sequenced reads
Pheno_result_metrics$mt_gene_pct <- round((Pheno_result_metrics$mt_gene_reads/Pheno_result_metrics$Genic_reads)*100,2) #percentage of mitochondrial reads vs all genic (intron + exon) reads

### plot the metrics

Pheno_result_metrics$phenotype <- as.factor(Pheno_result_metrics$phenotype) #change phenotype to factor
Pheno_result_metrics$phenotype <- factor(Pheno_result_metrics$phenotype, levels(Pheno_result_metrics$phenotype)[c(1,4,2,3)]) #reorder factors to the order that we want to view them in

library(ggplot2)

ggplot(Pheno_result_metrics, aes(x = phenotype, y = mt_gene_pct)) + geom_boxplot() + geom_jitter(width = 0.1) #template/example code for plotting any metric as jitter plot grouped by phenotype; in this case, percentage of mitochondrial reads

ggplot(Pheno_result_metrics, aes(x = Sequenced_reads, y = All_detected_genes)) + geom_point() + geom_smooth(method = "lm") #template/example code for plotting two metrics against each other; in this case, number of input/sequenced reads vs the number of all detected genes
ggplot(Pheno_result_metrics, aes(x = Unmapped_pct, y = All_detected_genes)) + geom_point() + geom_smooth(method = "lm") #template/example code for plotting two metrics against each other; in this case, percentage of unmapped reads vs the number of all detected genes

# summary plots
library(tidyr)
gathered_Pheno_result_metrics <- Pheno_result_metrics[,c(1,4:6,15,13,14)] #subsetting data to the metrics/reads we that we want to plot; how sequenced reads are subdivided, etc.
gathered_Pheno_result_metrics <- gather(gathered_Pheno_result_metrics, "Read_type", "Read_amount", c(3:7)) #format the data so we can do facet wrap 
gathered_Pheno_result_metrics$Read_type <- gsub("_reads", "", gathered_Pheno_result_metrics$Read_type) #clean the "Read_type" to remove "_reads" from the values
gathered_Pheno_result_metrics$Read_type <- ordered(gathered_Pheno_result_metrics$Read_type, levels = c("Exonic", "Intronic", "Intergenic", "Multi_mapped", "Unmapped")) #reorder the factors so that they're displayed in the order that we want

# plot for absolute number of reads
ggplot() +
  geom_bar(data=gathered_Pheno_result_metrics, aes(y = Read_amount, x = sample, fill = Read_type), stat="identity",
           position='stack') +
  theme_bw() + 
  facet_wrap( ~ phenotype, scales ="free_x") +
  scale_fill_manual(values = c("red", "blue", "yellow", "purple", "grey")) +
  ylab("Number of reads") + xlab("Sample ID") + labs(fill = "Type of read")

# plot for pct/proportion of reads
ggplot() +
  geom_bar(data=gathered_Pheno_result_metrics, aes(y = Read_amount, x = sample, fill = Read_type), stat="identity",
           position='fill') +
  theme_bw() + 
  facet_wrap( ~ phenotype, scales ="free") +
  scale_fill_manual(values = c("red", "blue", "yellow", "purple", "grey")) +
  scale_y_continuous(labels = scales::percent_format()) +
  ylab("Read-type percentage") + xlab("Sample ID") + labs(fill = "Type of read")

### export the csv

write.csv(Pheno_result_metrics, file = "/external/rprshnas01/kcni/ychen/git/Vgatcollab/Output_results/Pheno_result_metrics.csv")
