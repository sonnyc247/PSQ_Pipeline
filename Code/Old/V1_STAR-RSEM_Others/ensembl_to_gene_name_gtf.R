
mouse_ensembl_gtf_anno <- read.csv("/external/rprshnas01/kcni/ychen/git/PSQ_Pipeline/Data/mouse_ensembl_gtf_anno.csv", row.names=1) # load mouse ensembl gtf anno

count_matrix <- read.csv("/external/rprshnas01/netdata_kcni/stlab/File_transfer/Mini_Test_Data/test_output/RSEM_results/gene_count_matrix.csv") # load count matrix to convert

converted_count_matrix <- merge(mouse_ensembl_gtf_anno[,c("gene_name", "gene_id")], count_matrix, by.x = "gene_id", by.y = "X") # merge dataframes

row.names(converted_count_matrix) <- converted_count_matrix$gene_name # assign gene names to row names

converted_count_matrix <- converted_count_matrix[,3:ncol(converted_count_matrix)] # remove row names and ids 

write.csv(count_matrix, file = "/external/rprshnas01/netdata_kcni/stlab/File_transfer/Mini_Test_Data/test_output/RSEM_results/gene_name_count_mtx_.csv") # save converted matrix as separate file
