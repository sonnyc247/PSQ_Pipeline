#!/usr/bin/env Rscript

# script that reads in a directory of count matrices (csv format), divides them
# by their experimental group (based on metadata in xlsx file) and then plots the data accordingly

# Arguments should be the metadata file followed by the dir containing the csv files

# to run script in command line type:
# "Rscript plot_csv.R arg1 arg2"
# where arg1 is the xlsx metadata file and arg2 is the directory containing csv files

library(xlsx)
library(ggplot2)

# reading in and checking arguments from the command line
args = commandArgs(trailingOnly = TRUE)

if (length(args) > 2 | length(args) < 2) {
  stop("2 arguments must be supplied: 
       1. xlsx file containing metadata 
       2. path to directory containing csv files", call.=FALSE)
} else if (!dir.exists(csv_dir)) {
  stop("Invalid csv directory", call.False)
}

pheno <- "/external/rprshnas01/netdata_kcni/stlab/kristina_patchseq/Phenodata.xlsx"

csv_dir <- "/external/rprshnas01/netdata_kcni/stlab/kristina_patchseq_processed/STAR_results/sortedbams/CSV/"

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

for (i in csv_files) {
  
  if (gsub("[^0-9]", "", i) %in% flr$sample) {
    
    count_matrix <- read.csv(file = paste(csv_dir, "/", i, sep = ""))
    
    temp <- data.frame(sample = rep(gsub("[^0-9]", "", i), 2),
                       condition = c("exons", "introns"),
                       value = c(sum(count_matrix$exons), sum(count_matrix$introns)))
    
    flr_stack <- rbind(flr_stack, temp)
    
  } else if (gsub("[^0-9]", "", i) %in% nonflr$sample) {
    
    count_matrix <- read.csv(file = paste(csv_dir, "/", i, sep = ""))
    
    temp <- data.frame(sample = rep(gsub("[^0-9]", "", i), 2),
                       condition = c("exons", "introns"),
                       value = c(sum(count_matrix$exons), sum(count_matrix$introns)))
    
    nonflr_stack <- rbind(nonflr_stack, temp)
    
  } else if (gsub("[^0-9]", "", i) %in% ctrlwp$sample) {
    
    count_matrix <- read.csv(file = paste(csv_dir, "/", i, sep = ""))
    
    temp <- data.frame(sample = rep(gsub("[^0-9]", "", i), 2),
                       condition = c("exons", "introns"),
                       value = c(sum(count_matrix$exons), sum(count_matrix$introns)))
    
    ctrlwp_stack <- rbind(ctrlwp_stack, temp)
    
  } else if (gsub("[^0-9]", "", i) %in% ctrlwop$sample) {
    
    count_matrix <- read.csv(file = paste(csv_dir, "/", i, sep = ""))
    
    temp <- data.frame(sample = rep(gsub("[^0-9]", "", i), 2),
                       condition = c("exons", "introns"),
                       value = c(sum(count_matrix$exons), sum(count_matrix$introns)))
    
    ctrlwop_stack <- rbind(ctrlwop_stack, temp)
    
  }
}

# adding group column to dfs
flr_stack$group <- rep("flr", nrow(flr_stack))
nonflr_stack$group <- rep("nonflr", nrow(nonflr_stack))
ctrlwp_stack$group <- rep("ctrlwp", nrow(ctrlwp_stack))
ctrlwop_stack$group <- rep("ctrlwop", nrow(ctrlwop_stack))

group_stack <- do.call("rbind", list(flr_stack, nonflr_stack, ctrlwp_stack, ctrlwop_stack))
print(group_stack)

jpeg(paste(csv_dir, "/", "group_stacked.jpg", sep = ""))

ggplot(group_stack, aes(fill = condition, y = value, x = sample)) +
  geom_bar(position = "stack", stat = "identity") +
  xlab("Sample") + ylab("Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap( ~ group)

dev.off()

jpeg(paste(csv_dir, "/", "group_percent_stacked.jpg", sep = ""))

ggplot(group_stack, aes(fill = condition, y = value, x = sample)) +
  geom_bar(position = "fill", stat = "identity") +
  xlab("Sample") + ylab("Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap( ~ group)

dev.off()

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
