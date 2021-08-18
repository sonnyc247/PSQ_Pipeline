#!/usr/bin/env Rscript

# script that reads in a directory of count matrices (csv format), divides them
# by their experimental group (based on metadata in xlsx file) and then plots the data accordingly

# Arguments should be the metadata file followed by the dir containing the csv files

# to run script in command line type:
# "Rscript plot_csv.R arg1 arg2"
# where arg1 is the xlsx metadata file and arg2 is the directory containing csv files

library(xlsx)
library(ggplot2)
library(data.table)

# reading in and checking arguments from the command line
args = commandArgs(trailingOnly = TRUE)

if (length(args) > 2 | length(args) < 2) {
  stop("2 arguments must be supplied:
       1. xlsx file containing metadata
       2. path to directory containing csv files", call.=FALSE)
} else if (!dir.exists(csv_dir)) {
  stop("Invalid csv directory", call.False)
}

# hard code directories so you can go through plots line by line 
pheno <- "/external/rprshnas01/netdata_kcni/stlab/kristina_patchseq/Phenodata.xlsx"

csv_dir <- "/external/rprshnas01/netdata_kcni/stlab/kristina_patchseq_processed/STAR_results/pcrless_bams/InExCounts/"

# loading in metadata
meta <- read.xlsx(file = pheno, sheetIndex = 1)

# sample names so they can be compared to csv files
meta$sample <- gsub("[^0-9]", "", meta$sample)

# separating metadata by group
flr <- meta[grep("^fluorescent$", meta$phenotype), ]
nonflr <- meta[grep("not-fluorescent", meta$phenotype), ]
ctrlwp <- meta[grep("nc-w-pr", meta$phenotype), ]
ctrlwop <- meta[grep("nc-wo-pr", meta$phenotype), ]

# list of csv files in given directory
csv_files <- list.files(path = csv_dir, pattern = "\\.csv$")

# ignore csv produced by script combine_csv.R
csv_files <- grep(csv_files, pattern='^combined', inv=T, value=T)

# initializing dataframes to hold plot data for each group
flr_stack <- data.frame(sample = character(), condition = character(), value = integer())
nonflr_stack <- data.frame(sample = character(), condition = character(), value = integer())
ctrlwp_stack <- data.frame(sample = character(), condition = character(), value = integer())
ctrlwop_stack <- data.frame(sample = character(), condition = character(), value = integer())

# reading in combined count matrices
combined_exons = read.csv(paste0(csv_dir,'combined_exons.csv'))
combined_introns = read.csv(paste0(csv_dir,'combined_introns.csv'))
combined_genes = read.csv(paste0(csv_dir,'combined_genes.csv'))

# counting number of genes detected (counts > 0)
gene_count_exon = colSums(combined_exons[,3:38] > 0, na.rm = T)
gene_count_intron = colSums(combined_introns[,3:38] > 0, na.rm = T)
gene_count_exons_introns = colSums(combined_genes[,3:38] > 0, na.rm = T)

# separating gene counts based on location
gene_count_exon_only <- colSums((combined_exons[,3:38] > 0 & 
                                           combined_introns[,3:38] == 0),
                                        na.rm = T)
gene_count_intron_only <- colSums((combined_exons[,3:38] == 0 & 
                                           combined_introns[,3:38] > 0),
                                        na.rm = T)
gene_count_exons_introns_only <- colSums((combined_exons[,3:38] > 0 & 
                                           combined_introns[,3:38] > 0),
                                        na.rm = T)

# reformatting df to contain location and sample group info
exon_only_df <- data.frame("gene_counts" = gene_count_exon_only, 
                           "location" = rep("exons only", length(gene_count_exon_only)))
setDT(exon_only_df, keep.rownames = "sample")
intron_only_df <- data.frame("gene_counts" = gene_count_intron_only, 
                           "location" = rep("introns only", length(gene_count_intron_only)))
setDT(intron_only_df, keep.rownames = "sample")
exon_intron_df <- data.frame("gene_counts" = gene_count_exons_introns_only, 
                           "location" = rep("exons + introns", length(gene_count_exons_introns_only)))
setDT(exon_intron_df, keep.rownames = "sample")

# combined df so that stacked barplots can be generated
gene_count_stack_df <- do.call("rbind", list(exon_only_df, intron_only_df, exon_intron_df))
gene_count_stack_df$sample <- gsub("[^0-9]", "", gene_count_stack_df$sample)

gene_count_stack_df$group <- "temp"

# adding sample group info
gene_count_stack_df[gene_count_stack_df$sample %in% flr$sample]$group <- "flr"
gene_count_stack_df[gene_count_stack_df$sample %in% nonflr$sample]$group <- "nonflr"
gene_count_stack_df[gene_count_stack_df$sample %in% ctrlwp$sample]$group <- "ctrlwp"
gene_count_stack_df[gene_count_stack_df$sample %in% ctrlwop$sample]$group <- "ctrlwop"

# plotting stacked gene counts
jpeg(paste(csv_dir, "/", "gene_counts_stacked.jpg", sep = ""))

ggplot(gene_count_stack_df, aes(fill=location, y=gene_counts, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  xlab("Sample") + ylab("Gene Count") + ggtitle("Allen Institute PCRless Gene Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap( ~ group, scales = "free_x") + 
  guides(fill=guide_legend(title="Read location"))

dev.off()

# plotting percent stacked gene counts
jpeg(paste(csv_dir, "/", "gene_counts_percent_stacked.jpg", sep = ""))

ggplot(gene_count_stack_df, aes(fill=location, y=gene_counts, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  xlab("Sample") + ylab("Gene Count") + ggtitle("Allen Institute PCRless Gene Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap( ~ group, scales = "free_x") + 
  guides(fill=guide_legend(title="Read location"))

dev.off()

# gene counts exons vs exons + introns
gene_count_df <- data.frame("gene_count_exon_introns_only" = gene_count_exon_introns_only,
                            "gene_count_exon_only" = gene_count_exon_only,
                            "gene_count_intron_only" = gene_count_intron_only)

gene_count_df %>% ggplot(aes(x = gene_count_exon_introns_only, y = gene_count_exon_only)) +
  geom_point()

read_count_exon_only = colSums(combined_exons[,3:38], na.rm = T)
read_count_intron_only = colSums(combined_introns[,3:38], na.rm = T)
read_count_exons_introns = colSums(combined_genes[,3:38], na.rm = T)

meta_sample_df = data.frame('gene_count_exon_only' = gene_count_exon_only,
           'gene_count_intron_only' = gene_count_intron_only,
           'gene_count_exons_introns' = gene_count_exons_introns,
           'read_count_exon_only' = read_count_exon_only,
           'read_count_intron_only' = read_count_intron_only,
           'read_count_exons_introns' = read_count_exons_introns
           )
meta_sample_df$sample = gsub("[^0-9]", "", rownames(meta_sample_df))
rownames(meta_sample_df) = meta_sample_df$sample %>% make.names()

meta_sample_df %>% ggplot(aes(x = gene_count_exons_introns, y = gene_count_exon_only)) +
  geom_point()

# $exons, $introns for Rsubread, $exon_count, $intron_count for Allen Institute
for (i in csv_files) {
  
  if (gsub("[^0-9]", "", i) %in% flr$sample) {
    
    count_m <- read.csv(file = paste(csv_dir, "/", i, sep = ""))
    
    temp <- data.frame(sample = rep(gsub("[^0-9]", "", i), 2),
                       condition = c("exons", "introns"),
                       value = c(sum(count_m$exon_count), sum(count_m$intron_count)))
    
    flr_stack <- rbind(flr_stack, temp)
    
  } else if (gsub("[^0-9]", "", i) %in% nonflr$sample) {

    count_m <- read.csv(file = paste(csv_dir, "/", i, sep = ""))
    
    temp <- data.frame(sample = rep(gsub("[^0-9]", "", i), 2),
                       condition = c("exons", "introns"),
                       value = c(sum(count_m$exon_count), sum(count_m$intron_count)))
    
    nonflr_stack <- rbind(nonflr_stack, temp)
    
  } else if (gsub("[^0-9]", "", i) %in% ctrlwp$sample) {

    count_m <- read.csv(file = paste(csv_dir, "/", i, sep = ""))
    
    temp <- data.frame(sample = rep(gsub("[^0-9]", "", i), 2),
                       condition = c("exons", "introns"),
                       value = c(sum(count_m$exon_count), sum(count_m$intron_count)))
    
    ctrlwp_stack <- rbind(ctrlwp_stack, temp)
    
  } else if (gsub("[^0-9]", "", i) %in% ctrlwop$sample) {

    count_m <- read.csv(file = paste(csv_dir, "/", i, sep = ""))
    
    temp <- data.frame(sample = rep(gsub("[^0-9]", "", i), 2),
                       condition = c("exons", "introns"),
                       value = c(sum(count_m$exon_count), sum(count_m$intron_count)))
    
    ctrlwop_stack <- rbind(ctrlwop_stack, temp)
    
  }
}

# adding group column to dfs
flr_stack$group <- rep("flr", nrow(flr_stack))
nonflr_stack$group <- rep("nonflr", nrow(nonflr_stack))
ctrlwp_stack$group <- rep("ctrlwp", nrow(ctrlwp_stack))
ctrlwop_stack$group <- rep("ctrlwop", nrow(ctrlwop_stack))

group_stack <- do.call("rbind", list(flr_stack, nonflr_stack, ctrlwp_stack, ctrlwop_stack))

jpeg(paste(csv_dir, "/", "group_stacked.jpg", sep = ""))

ggplot(group_stack, aes(fill = condition, y = value, x = sample)) +
  geom_bar(position = "stack", stat = "identity") +
  xlab("Sample") + ylab("Counts") + ggtitle("Allen Institute PCRless") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap( ~ group, scales = "free_x") +
  guides(fill=guide_legend(title="Read location"))

dev.off()

jpeg(paste(csv_dir, "/", "group_percent_stacked.jpg", sep = ""))

ggplot(group_stack, aes(fill = condition, y = value, x = sample)) +
  geom_bar(position = "fill", stat = "identity") +
  xlab("Sample") + ylab("Ratio") + ggtitle("Allen Institute PCRless") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap( ~ group, scales = "free_x") + 
  guides(fill=guide_legend(title="Read location"))

dev.off()

# meta_sample_df
ggplot(merge(meta_sample_df, group_stack %>% dplyr::filter(condition == 'exons'), by = 'sample'), aes(x = sample, y = gene_count_exons_introns)) +
  geom_bar(stat = "identity") + facet_wrap(~group, scales = "free_x") + ylab("Unique genes")

# preliminary plotting examples

# # plotting flr samples
# jpeg(paste(csv_dir, "/", "flr_stacked.jpg", sep = ""))
# 
# ggplot(flr_stack, aes(fill = condition, y = value, x = sample)) + 
#   geom_bar(position = "stack", stat = "identity") + 
#   xlab("Sample") + ylab("Counts") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# dev.off()
# 
# jpeg(paste(csv_dir, "/", "flr_percent_stacked.jpg", sep = ""))
# 
# ggplot(flr_stack, aes(fill = condition, y = value, x = sample)) + 
#   geom_bar(position = "fill", stat = "identity") + 
#   xlab("Sample") + ylab("percentage") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# dev.off()
# 
# # plotting nonflr samples
# jpeg(paste(csv_dir, "/", "nonflr_stacked.jpg", sep = ""))
# 
# ggplot(nonflr_stack, aes(fill = condition, y = value, x = sample)) + 
#   geom_bar(position = "stack", stat = "identity") + 
#   xlab("Sample") + ylab("Counts") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# dev.off()
# 
# jpeg(paste(csv_dir, "/", "nonflr_percent_stacked.jpg", sep = ""))
# 
# ggplot(nonflr_stack, aes(fill = condition, y = value, x = sample)) + 
#   geom_bar(position = "fill", stat = "identity") + 
#   xlab("Sample") + ylab("percentage") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# dev.off()
# 
# # plotting ctrlwp samples
# jpeg(paste(csv_dir, "/", "ctrlwp_stacked.jpg", sep = ""))
# 
# ggplot(ctrlwp_stack, aes(fill = condition, y = value, x = sample)) + 
#   geom_bar(position = "stack", stat = "identity") + 
#   xlab("Sample") + ylab("Counts") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# dev.off()
# 
# jpeg(paste(csv_dir, "/", "ctrlwp_percent_stacked.jpg", sep = ""))
# 
# ggplot(ctrlwp_stack, aes(fill = condition, y = value, x = sample)) + 
#   geom_bar(position = "fill", stat = "identity") + 
#   xlab("Sample") + ylab("percentage") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# dev.off()
# 
# # plotting ctrlwop samples
# jpeg(paste(csv_dir, "/", "ctrlwop_stacked.jpg", sep = ""))
# 
# ggplot(ctrlwop_stack, aes(fill = condition, y = value, x = sample)) + 
#   geom_bar(position = "stack", stat = "identity") + 
#   xlab("Sample") + ylab("Counts") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# dev.off()
# 
# jpeg(paste(csv_dir, "/", "ctrlwop_percent_stacked.jpg", sep = ""))
# 
# ggplot(ctrlwop_stack, aes(fill = condition, y = value, x = sample)) + 
#   geom_bar(position = "fill", stat = "identity") + 
#   xlab("Sample") + ylab("percentage") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# dev.off()
