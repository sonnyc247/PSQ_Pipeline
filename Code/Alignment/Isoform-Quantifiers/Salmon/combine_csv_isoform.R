#### salmon ####

csv_dir <- "/external/rprshnas01/netdata_kcni/stlab/File_transfer/RefSeqTest/Pilot_PVALBvsL23_inVisp/Pilot_quant/salmon_batch_results/"

# list of dir in given directory
csv_files <- list.dirs(path = csv_dir, recursive = F)

# ignore csv produced by script combine_csv.R
#csv_files <- grep(csv_files, pattern='^combined', inv=T, value=T)

# for every file
for (i in 1:length(csv_files)) {
  
  curr_file <- csv_files[i] # get current file 
  sample_matrix <- read.delim(file = paste(curr_file, "/quant.sf", sep = "")) # intake current file table
  sample_name <- basename(curr_file) # get sample name (approximate)
  sample_name <- substr(sample_name, 1, nchar(sample_name)-6)
  colnames(sample_matrix)[5] <- sample_name # set sample name over expected counts column
  sample_matrix <- sample_matrix[,c("Name", sample_name)] # take columns of interest
  
  if (i == 1) {
    
    count_matrix <- sample_matrix
    
  } else {
    
    count_matrix <- merge(count_matrix, sample_matrix, by = "Name")
    
  }
  
}

row.names(count_matrix) <- count_matrix$Name # set row names as gene names
count_matrix <- count_matrix[-1] # remove gene name column to make object pure matrix

# writing combined count matrices
write.csv(count_matrix, file = "salmon_count_mtx.csv")

