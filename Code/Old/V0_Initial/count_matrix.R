#!/usr/bin/env Rscript

# script that produces count matrix from SINGLE bam file
# arguments should be: 
# 1. path to bamfile 2. path to gtf file 3. output directory for csv files

# to run script on command line type in:
# "Rscript count_matrix.R arg1 arg2 arg3"
# e.g.
# "Rscript count_matrix.R datafiles/mutant_female_Grin1.bam datafiles/Mus_musculus.GRCm38.98.gtf datafiles/csv"

library("Rsubread")
library("dplyr")

# function to create dataframe with exon, gene and intron data
combine_counts <- function(e_count, g_count) {
  
  # merge the exon and gene data
  temp <- merge(e_count, g_count, by = "row.names")
  
  # dataframe to store count data 
  df <- data.frame(Gene_ID = character(),
                   exons = integer(), 
                   genes = integer(), 
                   introns = integer())
  
  for (i in 1:length(temp[,1])) {
    if (temp[i, 2] != 0 && temp[i, 3] != 0) {
      if (temp[i, 3] - temp[i, 2] >= 0) {
        df <- rbind(df, data.frame("Gene_ID" = temp[i, 1], "exons" = temp[i, 2], 
                                   "genes" = temp[i, 3], "introns" = temp[i, 3] - temp[i, 2]))
      } else {
        df <- rbind(df, data.frame("Gene_ID" = temp[i, 1], "exons" = temp[i, 2], 
                                   "genes" = temp[i, 3], "introns" = 0))
      }
    }
  }
  
  return(df)
}

# reading in and checking arguments from the command line
args = commandArgs(trailingOnly = TRUE)

if (length(args) > 3 | length(args) < 3) {
  stop("3 arguments must be supplied: 
       1. path to bamfile 2. path to gtf file 3. output directory for csv files", call.=FALSE)
} else if (!dir.exists(args[3])) {
  stop("Invalid csv directory", call.False)
}

# variable holding exon counts
temp_e <- featureCounts(files = args[1], 
                        annot.ext = args[2],
                        isGTFAnnotationFile = TRUE, GTF.featureType="exon", GTF.attrType="gene_id", 
                        isPairedEnd = TRUE)
# variable holding gene counts
temp_g <- featureCounts(files = args[1], 
                        annot.ext = args[2],
                        isGTFAnnotationFile = TRUE, GTF.featureType="gene", GTF.attrType="gene_id", 
                        isPairedEnd = TRUE)

# temp dataframe holding exon, gene and intron counts
temp_df <- combine_counts(temp_e$counts, temp_g$counts)

# write csv file containing count matrix
write.csv(temp_df, file = paste(args[3], "/", 
                                tools::file_path_sans_ext(basename(args[1])), 
                                "_count.csv", 
                                sep = ""))
