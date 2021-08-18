
csv_dir <- "/external/rprshnas01/netdata_kcni/stlab/File_transfer/Mini_Test_Data/test_output/RSEM_results/Gene_counts/"

# list of files in given directory
csv_files <- list.files(path = csv_dir, pattern = "\\.results$")

# ignore csv produced by script combine_csv.R
#csv_files <- grep(csv_files, pattern='^combined', inv=T, value=T)

# for every file
for (i in 1:length(csv_files)) {

  curr_file <- csv_files[i] # get current file 
  sample_matrix <- read.delim(file = paste(csv_dir, curr_file, sep = "")) # intake current file table
  sample_name <- gsub("\\..*","", curr_file) # get sample name (approximate)
  colnames(sample_matrix)[5] <- sample_name # set sample name over expected counts column
  sample_matrix <- sample_matrix[,c("gene_id", sample_name)] # take columns of interest
  
  if (i == 1) {
    
  count_matrix <- sample_matrix
    
  } else {
    
  count_matrix <- merge(count_matrix, sample_matrix, by = "gene_id")
    
  }
  
}

row.names(count_matrix) <- count_matrix$gene_id # set row names as gene names
count_matrix <- count_matrix[-1] # remove gene name column to make object pure matrix

# writing combined count matrices
write.csv(count_matrix, file = "/external/rprshnas01/netdata_kcni/stlab/File_transfer/Mini_Test_Data/test_output/RSEM_results/gene_count_matrix.csv")
