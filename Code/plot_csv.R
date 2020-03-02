#!/usr/bin/env Rscript

# script that plots a directory of csv files produced by "count_matrix.R" or "count_matrix_dir.R"

# Argument should be the dir containing the csv files (individual count matrices)

# to run script in command line type:
# "Rscript plot_csv.R arg1"
# where arg1 is the directory containing csv files

library(ggplot2)

# reading in and checking arguments from the command line
args = commandArgs(trailingOnly = TRUE)

if (length(args) > 1 | length(args) < 1) {
  stop("1 argument must be supplied: 1. path to directory containing csv files", call.=FALSE)
} else if (!dir.exists(args[1])) {
  stop("Invalid csv directory", call.False)
}

# list of csv files in given directory
csv_files <- list.files(path = args[1], pattern = "\\.csv$")

# ignore csv produced by script combine_csv.R
csv_files <- grep(csv_files, pattern='combined_counts.csv', inv=T, value=T)

# initialize dataframe containing stacked barplot data
stack_df <- data.frame(sample = character(), condition = character(), value = integer())

for (i in csv_files) {
  
  count_matrix <- read.csv(file = paste(args[1], "/", i, sep = ""))
  
  temp <- data.frame(sample = rep(gsub("[^0-9]", "", i), 2),
                     condition = c("exons", "introns"),
                     value = c(sum(count_matrix$exons), sum(count_matrix$introns)))
  
  stack_df <- rbind(stack_df, temp)
}

jpeg(paste(args[1], "/", "stacked.jpg", sep = ""))

ggplot(stack_df, aes(fill = condition, y = value, x = sample)) + 
  geom_bar(position = "stack", stat = "identity") + 
  xlab("Sample") + ylab("Counts") + ggtitle("Stacked Counts")
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

jpeg(paste(args[1], "/", "percent_stacked.jpg", sep = ""))

ggplot(stack_df, aes(fill = condition, y = value, x = sample)) + 
  geom_bar(position = "fill", stat = "identity") + 
  xlab("Sample") + ylab("percentage") + ggtitle("Percentage Counts")
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
