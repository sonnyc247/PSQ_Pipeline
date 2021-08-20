### for use generating coverage plots; use this script after data_load script and Ensembl_to_Mgi, also need coord bams from STAR
### check out "fromjustin.R" for where some of the code below comes from (mainly for genomic alignments library work)

###get totals for nomalizing counts (by library size minus mitochondrial counts)

#get Mgi names of genes
Mgi_count_results <- merge(x = count_results, y = Gene_anno, by.x = "row.names", by.y = "Gene", all.x = TRUE)

#isolate mithondrial genes (that start with "mt-")
Mt_count_results <- Mgi_count_results[grep("mt-", Mgi_count_results$Mgi_Gene), ]
Mt_count_results <- Mt_count_results[-38,] #excluding Bhmt-ps1, which is not a mitochondrial gene
Mt_count_results <- rownames_to_column(Mt_count_results) 
Mt_count_results$rowname <- as.integer(Mt_count_results$rowname) #get indexes of mitochondrial genes
Mtless_count_results <-Mgi_count_results[-c(Mt_count_results$rowname),] #remove mitochondrial genes based on indices from Mt_count_results df
Unique_Mtless_count_results <- Mtless_count_results[-which(duplicated(Mtless_count_results$Row.names)),] #remove duplicate rows that came from ensembl-Mgi conversion
row.names(Unique_Mtless_count_results) <- Unique_Mtless_count_results$Row.names #set ensembl gene names back to row names
Uni_Mtles_counts_num <- Unique_Mtless_count_results[,c(2:29)] #get only numeric data (counts), so that we can use colsums
Mtless_libsize <- as.data.frame(colSums(Uni_Mtles_counts_num)) #sum the reads/counts of nonmitochondrial genes
cov_combined_df <- merge(combined_df, Mtless_libsize, by.x = "ShortRun", by.y = "row.names") #add mitochondria-less library size as a metric to sample-level, summative coverage df

### generate coverage for each bam file of interest

# set genomic region of interest

library("GenomicAlignments")
library(limma) #maybe unnecessary?

# setting genomic range

#grin1_range_N <- GRanges("2", IRanges(25291177, 25319253)) #in chrom two, these coord are for grin1 coordinates from NCBI
#grin1_range_E <- GRanges("2", IRanges(25291181, 25319187)) #in chrom two, these coord are for grin1 coordinates from ensembl
#interest_range <- GRanges("2", IRanges(25291181, 25295937)) #in chrom two, these coord are for grin1 exons 18-20 (introns included) from ensembl; mutation may be in intron 18-19
key_coord <- 25294729 #coordinate of key site of interest
subset_range <- 10000000
interest_range <- GRanges("2", IRanges(key_coord-subset_range, key_coord+subset_range))

# for every bam with sample name, generate coverage plot, stored into one dataframe

coord_bam_dir <- "/external/rprshnas01/netdata_kcni/stlab/RamseyMielnik/STAR_results/coord_bams/" #directory of bam files for processing
coverage_df <- matrix(0, nrow = 2*subset_range+1, ncol = length(cov_combined_df[,1])) #initialize a matrix (soon-to-be df) that spans the genomic region of interest
coverage_df <- as.data.frame(coverage_df) 
colnames(coverage_df) <- cov_combined_df[,1] #set sample names as column names
options(warn=1) #to output warning messages for each loop iteration, warnings about some samples may occur, need to keep in mind in case of outlier samples later

for (i in 1:length(cov_combined_df[,1])) {
  
  curr_sample_name = cov_combined_df[i,1] #defining current sample name, file path, etc.
  curr_sample_path = paste0(coord_bam_dir, curr_sample_name, "_1Aligned.sortedByCoord.out.bam") #this is the full file path for the ith sample
  curr_coord_bam <- readGAlignmentPairs(curr_sample_path, param=ScanBamParam(which=interest_range)) #reading in bam file
  print(curr_sample_path) #for troubleshooting only (see which files/samples give you warning messages)
  curr_cov <- coverage(curr_coord_bam) #computing coverage for given bam file
  curr_cov[interest_range] 
  curr_cov_num <- as.numeric(curr_cov[interest_range][[1]]) #converting coverage into a readable format
  curr_cov_norm <- curr_cov_num/cov_combined_df$`colSums(Uni_Mtles_counts_num)`[i] #normalize coverage with mitochondria-less library size
  coverage_df[,i] <- curr_cov_norm #add coverage to result df

  }

# generating median coverage tracks for each (key) cell type and other information tracks (mutation point)

library(matrixStats)
library(tibble)

Layer5_track <- rowMedians(as.matrix(coverage_df[,12:15])) #columns defined by information from cov_combined_df, which specifies which samples are from which source
Gaba_track <- rowMedians(as.matrix(coverage_df[,5:7]))
Layer2_3_track <- rowMedians(as.matrix(coverage_df[,1:4]))
Layer6_track <- rowMedians(as.matrix(coverage_df[,8:11]))
Layer4_track <- rowMedians(as.matrix(coverage_df[,18:21]))
Composite_track <- data.frame(Layer6_track, Layer5_track, Layer2_3_track, Layer4_track, Gaba_track) #combining tracks into one dataframe
Composite_track <- rownames_to_column(Composite_track) #generate genomic coordinates as a track/variable (to plot on x axis) 
colnames(Composite_track)[1] <- "Genomic_coordinates"
Composite_track$Genomic_coordinates <- as.numeric(Composite_track$Genomic_coordinates)
Composite_track$Genomic_coordinates <- Composite_track$Genomic_coordinates+(key_coord-subset_range-1)
Composite_track$Mutation <- 0.000000
Composite_track$Mutation[which(Composite_track$Genomic_coordinates==key_coord)] <- 0.01 #more of a placeholder, will need to adjust based on the range being plotted/quantified
grin1_end_index <- which(Composite_track$Genomic_coordinates=="25291181")
grin1_start_index <- which(Composite_track$Genomic_coordinates=="25319187")
Composite_track_CPM <- Composite_track*1000000

### plotting data

# rough plots of coverage 

plot(Layer5_track, type="h")
plot(Gaba_track, type="h")

# with ggplot

library(ggplot2)

one_tail_range = 5000

#overlap line plots for some bp around insertion point

Composite_track_CPM$Mutation[which(Composite_track$Genomic_coordinates==key_coord)] <- max(Composite_track_CPM[(subset_range+1-one_tail_range):(subset_range+1+one_tail_range),c(3,6)]) #as needed, set insertion point to max of range

ggplot(Composite_track_CPM[(subset_range+1-one_tail_range):(subset_range+1+one_tail_range),], aes(Genomic_coordinates)) + 
  geom_line(aes(y = Layer5_track, colour = "Layer5")) + 
  geom_line(aes(y = Gaba_track, colour = "Gaba")) +
  geom_line(aes(y = Mutation, colour = "Insertion Point")) +
  scale_color_manual(values = c("black", "blue", "red")) +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('Normalized Genomic Accessibility (CPM)') 

#overlap line plots for whole gene

Composite_track_CPM$Mutation[which(Composite_track$Genomic_coordinates==key_coord)] <- max(Composite_track_CPM[c(grin1_end_index:grin1_start_index),c(2:6)]) #as needed, set insertion point to max of range

ggplot(Composite_track_CPM[grin1_end_index:grin1_start_index,], aes(Genomic_coordinates)) + 
  geom_line(aes(y = Layer5_track, colour = "L5")) + 
  geom_line(aes(y = Gaba_track, colour = "GABA")) +
  geom_line(aes(y = Layer6_track, colour = "L6")) +
  geom_line(aes(y = Layer2_3_track, colour = "L2/3")) +
  geom_line(aes(y = Layer4_track, colour = "L4")) +
  geom_line(aes(y = Mutation, colour = "Insertion Point")) +
  #geom_area(aes(y = Exonjunctions, fill = "Ending Exons"), alpha = 0.4) +
  #scale_color_manual(values = c("red", "black", "#3dad44", "#218a28", "#5cd664", "#7cfc84")) +
  scale_color_manual(values = c("red", "black", "green", "yellow", "cyan", "blue")) +
  #scale_fill_manual(values = "green") +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('Frequency') +
  theme_classic()

#area plots for whole gene, for each cell type; for stacking together

Composite_track_CPM$Mutation[which(Composite_track$Genomic_coordinates==key_coord)] <- max(Composite_track_CPM[(grin1_end_index:grin1_start_index),c(2:6)])

L5 <- ggplot(Composite_track_CPM[grin1_end_index:grin1_start_index,], aes(Genomic_coordinates)) + 
  geom_area(aes(y = Layer5_track, fill = "L5"), alpha = 0.4) + #plotting accessibility peaks
  scale_fill_manual(values = "green") + #setting fill colour
  geom_line(aes(y = Mutation, colour = "Insertion Point")) + #plotting insertion point
  scale_color_manual(values = c("black")) + #setting insertion point line colour
  guides(fill = FALSE, colour = FALSE) + #turn off legends
  xlab('Mouse Chromosome 2 Coordinate') + 
  ylab('L5') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) + #removing borders (otherwise there would be some more white space between axes and plots)
  scale_y_continuous(limits = c(0,93), expand = c(0, 0)) + #removing borders
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) #removing x axes ticks, labels, etc.

L4 <- ggplot(Composite_track_CPM[grin1_end_index:grin1_start_index,], aes(Genomic_coordinates)) + 
  geom_area(aes(y = Layer4_track, fill = "L4"), alpha = 0.4) + 
  scale_fill_manual(values = "green") + 
  geom_line(aes(y = Mutation, colour = "Insertion Point")) +
  scale_color_manual(values = c("black")) +
  guides(fill = FALSE, colour = FALSE) +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('L4') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,93), expand = c(0, 0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

GABA <- ggplot(Composite_track_CPM[grin1_end_index:grin1_start_index,], aes(Genomic_coordinates)) + 
  geom_area(aes(y = Gaba_track, fill = "GABA"), alpha = 0.4) + 
  scale_fill_manual(values = "magenta") + 
  geom_line(aes(y = Mutation, colour = "Insertion Point")) +
  scale_color_manual(values = c("black")) +
  guides(fill = FALSE, colour = FALSE) +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('GABA') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,93), expand = c(0, 0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

Xaxis <- ggplot(Composite_track_CPM[grin1_end_index:grin1_start_index,], aes(Genomic_coordinates)) + 
  geom_area(aes(y = Gaba_track, fill = "GABA"), alpha = 0.4) + 
  scale_fill_manual(values = "white") + 
  guides(fill = FALSE, colour = FALSE) +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('   ') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,93), expand = c(0, 0))
  

L6 <- ggplot(Composite_track_CPM[grin1_end_index:grin1_start_index,], aes(Genomic_coordinates)) + 
  geom_area(aes(y = Layer6_track, fill = "L6"), alpha = 0.4) + 
  scale_fill_manual(values = "green") + 
  geom_line(aes(y = Mutation, colour = "Insertion Point")) +
  scale_color_manual(values = c("black")) +
  guides(fill = FALSE, colour = FALSE) +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('L6') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,93), expand = c(0, 0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

L2_3 <- ggplot(Composite_track_CPM[grin1_end_index:grin1_start_index,], aes(Genomic_coordinates)) + 
  geom_area(aes(y = Layer2_3_track, fill = "L2/3"), alpha = 0.4) +
  scale_fill_manual(values = "green") + 
  geom_line(aes(y = Mutation, colour = "Insertion Point")) +
  scale_color_manual(values = c("black")) +
  guides(fill = FALSE, colour = FALSE) +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('L2/3') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,93), expand = c(0, 0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

library(ggpubr)
stacked_plot <- ggarrange(L4, L2_3, L5, L6, GABA, Xaxis, ncol = 1, nrow = 6)
annotate_figure(stacked_plot, fig.lab = "Normalized Genome Accessibility (CPM)", fig.lab.pos = "top")

#area plots for 1k range, for each cell type; for stacking together

Composite_track_CPM$Mutation[which(Composite_track$Genomic_coordinates==key_coord)] <- max(Composite_track_CPM[(subset_range+1-1000):(subset_range+1+1000),c(2:6)])

L5 <- ggplot(Composite_track_CPM[(subset_range+1-1000):(subset_range+1+1000),], aes(Genomic_coordinates)) + 
  geom_area(aes(y = Layer5_track, fill = "L5"), alpha = 0.4) + #plotting accessibility peaks
  scale_fill_manual(values = "green") + #setting fill colour
  geom_line(aes(y = Mutation, colour = "Insertion Point")) + #plotting insertion point
  scale_color_manual(values = c("black")) + #setting insertion point line colour
  guides(fill = FALSE, colour = FALSE) + #turn off legends
  xlab('Mouse Chromosome 2 Coordinate') + 
  ylab('L5') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) + #removing borders (otherwise there would be some more white space between axes and plots)
  scale_y_continuous(limits = c(0,16), expand = c(0, 0)) + #removing borders
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) #removing x axes ticks, labels, etc.

L4 <- ggplot(Composite_track_CPM[(subset_range+1-1000):(subset_range+1+1000),], aes(Genomic_coordinates)) + 
  geom_area(aes(y = Layer4_track, fill = "L4"), alpha = 0.4) + 
  scale_fill_manual(values = "green") + 
  geom_line(aes(y = Mutation, colour = "Insertion Point")) +
  scale_color_manual(values = c("black")) +
  guides(fill = FALSE, colour = FALSE) +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('L4') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,16), expand = c(0, 0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

GABA <- ggplot(Composite_track_CPM[(subset_range+1-1000):(subset_range+1+1000),], aes(Genomic_coordinates)) + 
  geom_area(aes(y = Gaba_track, fill = "GABA"), alpha = 0.4) + 
  scale_fill_manual(values = "magenta") + 
  geom_line(aes(y = Mutation, colour = "Insertion Point")) +
  scale_color_manual(values = c("black")) +
  guides(fill = FALSE, colour = FALSE) +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('GABA') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,16), expand = c(0, 0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

Xaxis <- ggplot(Composite_track_CPM[(subset_range+1-1000):(subset_range+1+1000),], aes(Genomic_coordinates)) +
  geom_area(aes(y = Gaba_track, fill = "GABA"), alpha = 0.4) + 
  scale_fill_manual(values = "white") + 
  guides(fill = FALSE, colour = FALSE) +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('   ') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,16), expand = c(0, 0))

L6 <- ggplot(Composite_track_CPM[(subset_range+1-1000):(subset_range+1+1000),], aes(Genomic_coordinates)) + 
  geom_area(aes(y = Layer6_track, fill = "L6"), alpha = 0.4) + 
  scale_fill_manual(values = "green") + 
  geom_line(aes(y = Mutation, colour = "Insertion Point")) +
  scale_color_manual(values = c("black")) +
  guides(fill = FALSE, colour = FALSE) +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('L6') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,16), expand = c(0, 0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

L2_3 <- ggplot(Composite_track_CPM[(subset_range+1-1000):(subset_range+1+1000),], aes(Genomic_coordinates)) + 
  geom_area(aes(y = Layer2_3_track, fill = "L2/3"), alpha = 0.4) +
  scale_fill_manual(values = "green") + 
  geom_line(aes(y = Mutation, colour = "Insertion Point")) +
  scale_color_manual(values = c("black")) +
  guides(fill = FALSE, colour = FALSE) +
  xlab('Mouse Chromosome 2 Coordinate') +
  ylab('L2/3') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,16), expand = c(0, 0)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

library(ggpubr)
stacked_plot <- ggarrange(L4, L2_3, L5, L6, GABA, Xaxis, ncol = 1, nrow = 6)
annotate_figure(stacked_plot, fig.lab = "Normalized Genome Accessibility (CPM)", fig.lab.pos = "top.right")
  
### quantification/integration of data area-under-the-graph 

cov_combined_df$Integrated_coverage_Ex18.19 <- colSums(coverage_df[which(Composite_track$Genomic_coordinates=="25294438"):which(Composite_track$Genomic_coordinates=="25295937"),]) #Integration of exons 18-19 around 25294729 (index # 3549)
cov_combined_df$Integrated_coverage_50 <- colSums(coverage_df[(subset_range+1-50):(subset_range+1+50),]) #Integration of 50bp around 25294729 (index # 3549)
cov_combined_df$Integrated_coverage_100 <- colSums(coverage_df[(subset_range+1-100):(subset_range+1+100),]) #Integration of 100bp around 25294729 (index # 3549)
cov_combined_df$Integrated_coverage_200 <- colSums(coverage_df[(subset_range+1-200):(subset_range+1+200),]) #Integration of 200bp around 25294729 (index # 3549)
cov_combined_df$Integrated_coverage_500 <- colSums(coverage_df[(subset_range+1-500):(subset_range+1+500),]) #Integration of 500bp around 25294729 (index # 3549)
cov_combined_df$Integrated_coverage_1k <- colSums(coverage_df[(subset_range+1-1000):(subset_range+1+1000),]) #Integration of 1kbp around 25294729 (index # 3549)
cov_combined_df$Integrated_coverage_2.2k <- colSums(coverage_df[(subset_range+1-2200):(subset_range+1+2200),]) #Integration of 2.2kbp before 25294729 (index # 3549) and up until end
cov_combined_df$Integrated_coverage_150 <- colSums(coverage_df[(subset_range+1-150):(subset_range+1+150),]) #Integration of 150bp around 25294729 (index # 3549)

cov_combined_df$Integrated_coverage_5k <- colSums(coverage_df[(subset_range+1-5000):(subset_range+1+5000),]) #Integration of 5kbp around 25294729
cov_combined_df$Integrated_coverage_10k <- colSums(coverage_df[(subset_range+1-10000):(subset_range+1+10000),]) #Integration of 10kbp around 25294729
cov_combined_df$Integrated_coverage_100k <- colSums(coverage_df[(subset_range+1-100000):(subset_range+1+100000),]) #Integration of 100kbp around 25294729
cov_combined_df$Integrated_coverage_1M <- colSums(coverage_df[(subset_range+1-1000000):(subset_range+1+1000000),]) #Integration of 1Mbp around 25294729
cov_combined_df$Integrated_coverage_10M <- colSums(coverage_df) #Integration of 10Mbp around 25294729
cov_combined_df$Integrated_coverage_gene <- colSums(coverage_df[grin1_end_index:grin1_start_index,]) #Integration of 1Mbp around 25294729

cov_combined_df[,c(11:24)] = cov_combined_df[,c(11:24)]*1000000

#tidying and organizing quantifications

library(tidyr)

cov_combined_df$source_name <- as.character(cov_combined_df$source_name)
cov_combined_df[c(16:17),4] <- "Rbp4_L5_Mi_500cell" #separate miseq from hiseq data
cov_combined_df$source_name <- gsub("_500cell", "", cov_combined_df$source_name) #cleaning up unneessary parts of terms
cov_combined_df$source_name <- gsub("Cux2_", "", cov_combined_df$source_name) 
cov_combined_df$source_name <- gsub("Gad2_", "", cov_combined_df$source_name) 
cov_combined_df$source_name <- gsub("Ntsr1_", "", cov_combined_df$source_name)
cov_combined_df$source_name <- gsub("Rbp4_", "", cov_combined_df$source_name)
cov_combined_df$source_name <- gsub("Scnn1a_", "", cov_combined_df$source_name) 
cov_combined_df$source_name <- as.factor(cov_combined_df$source_name)
levels(cov_combined_df$source_name) #organize levels; help with order of group graphing
cov_combined_df$source_name = factor(cov_combined_df$source_name,levels(cov_combined_df$source_name)[c(7,5,3,4,1,6,8,2)]) #as needed/appropriate for organizing facet-wrapped graphs

cov_combined_df_gathered <- gather(data = cov_combined_df, "Int_Range", "Quantification", c(8,9,11:24))  #gather data for facet-wrap plotting
cov_combined_df_gathered$Int_Range <- gsub("Integrated_coverage_", "", cov_combined_df_gathered$Int_Range) #clean the integration range terms
cov_combined_df_gathered$Int_Range <- as.factor(cov_combined_df_gathered$Int_Range)
levels(cov_combined_df_gathered$Int_Range)
cov_combined_df_gathered$Int_Range = factor(cov_combined_df_gathered$Int_Range,levels(cov_combined_df_gathered$Int_Range)[c(15,16,13,10,1,5,9,11,6,14,8,12,3,2,7,4)]) #organize factors for plotting order

# plot the quantifications

ggplot(data = cov_combined_df_gathered, mapping = aes(x = source_name, y = Quantification)) +
         geom_jitter(width = 0) + 
         geom_boxplot(alpha = 0.4) + 
         labs(x = "Cell type", y = "Quantification (CPM, Count, or TPM)") + 
         facet_wrap(~Int_Range, scales = "free") #for ALL cell types, and ALL integration ranges

ggplot(data = cov_combined_df[cov_combined_df$source_name %in% c("GABA", "L5"),], mapping = aes(x = source_name, y = Integrated_coverage_gene)) + 
  geom_jitter(width = 0) + 
  geom_boxplot(alpha = 0.4) + 
  labs(x = "Cell type", y = "Genomic Accessibility (CPM)") #for L5 vs GABA, for integration over gene

ggplot(data = cov_combined_df, mapping = aes(x = source_name, y = Integrated_coverage_gene)) + 
  geom_jitter(width = 0) + 
  geom_boxplot(alpha = 0.4) + 
  labs(x = "Cell type", y = "Genomic Accessibility (CPM)") #for L5 vs GABA, for integration over gene

ggplot(data = cov_combined_df_gathered[cov_combined_df_gathered$source_name %in% c("GABA", "L5", "L2/3", "L4", "L6"),], mapping = aes(x = source_name, y = Quantification)) + 
  geom_boxplot(alpha = 0.4, fill = c("green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta", "green", "green", "green", "green", "magenta")) + 
  geom_jitter(width = 0) + 
  labs(x = "Cell type", y = "Quantification (CPM, Count, or TPM)") + 
  facet_wrap(~Int_Range, scales = "free") #for neuron types, and ALL integration ranges

ggplot(data = cov_combined_df[cov_combined_df$source_name %in% c("GABA", "L5", "L2/3", "L4", "L6"),], mapping = aes(x = source_name, y = Integrated_coverage_gene)) + 
  geom_boxplot(alpha = 0.6, fill = c("green", "green", "green", "green", "magenta")) + 
  geom_jitter(width = 0) + 
  labs(x = "Cell origin", y = "Genomic Accessibility (CPM)") +
  #stat_boxplot(geom ='errorbar') + 
  theme_classic() #for neuron types at gene level

ggplot(data = cov_combined_df[cov_combined_df$source_name %in% c("GABA", "L5", "L2/3", "L4", "L6"),], mapping = aes(x = source_name, y = Integrated_coverage_1k)) + 
  geom_boxplot(alpha = 0.6, fill = c("green", "green", "green", "green", "magenta")) + 
  geom_jitter(width = 0) + 
  labs(x = "Cell origin", y = "Genomic Accessibility (CPM)") +
  #stat_boxplot(geom ='errorbar') + 
  theme_classic() #for neuron types at 1k range around insertion site

# grouping together the GLUT types

levels(cov_combined_df$broad_class) <- c(levels(cov_combined_df$broad_class),"Glutamatergic_Mi") #add new factor for separating out Mi-seq L5 data
cov_combined_df[c(16:17),3] <- as.factor("Glutamatergic_Mi") 
cov_combined_df$broad_class = factor(cov_combined_df$broad_class,levels(cov_combined_df$broad_class)[c(2,1,4,3,5)]) #factor reordered for graphing purposes

ggplot(data = cov_combined_df[cov_combined_df$broad_class %in% c("Glutamatergic", "GABAergic"),], mapping = aes(x = broad_class, y = Integrated_coverage_gene)) + 
  geom_boxplot(alpha = 0.6, fill = c("green", "magenta")) + 
  geom_jitter(width = 0) + 
  labs(x = "Cell type", y = "Genomic Accessibility (CPM)") +
  #stat_boxplot(geom ='errorbar') + 
  theme_classic() #gene-level coverage

ggplot(data = cov_combined_df[cov_combined_df$broad_class %in% c("Glutamatergic", "GABAergic"),], mapping = aes(x = broad_class, y = Integrated_coverage_1k)) + 
  geom_boxplot(alpha = 0.6, fill = c("green", "magenta")) + 
  geom_jitter(width = 0) + 
  labs(x = "Cell type", y = "Genomic Accessibility (CPM)") +
  #stat_boxplot(geom ='errorbar') + 
  theme_classic() #1k around insertion coverage

### Wilcoxon/Mann-whitney-u test

wilcox.test(cov_combined_df$Integrated_coverage_gene[12:15], cov_combined_df$Integrated_coverage_gene[5:7])
wilcox.test(cov_combined_df$Integrated_coverage_1k[12:15], cov_combined_df$Integrated_coverage_1k[5:7])
wilcox.test(cov_combined_df$Integrated_coverage_gene[c(1:4,8:15,18:21)], cov_combined_df$Integrated_coverage_gene[5:7])
wilcox.test(cov_combined_df$Integrated_coverage_1k[c(1:4,8:15,18:21)], cov_combined_df$Integrated_coverage_1k[5:7])

### confirming tack plots with Gviz

library(Gviz) #key
library("rtracklayer")
library("GenomicFeatures") #key
library("GenomicRanges")
library("GenomicAlignments")
library("Rsubread")
library("erer")
library("plyr")
options(ucscChromosomeNames=FALSE)
gtfref <- makeTxDbFromGFF("/external/rprshnas01/kcni/ychen/References/Ensembl/Mus_musculus.GRCm38.98.gtf", format="gtf") #reading gtf, for plotting
grtrack <- GeneRegionTrack(gtfref, showTitle=FALSE, background.title="white", col.axis="black")
gtrack <- GenomeAxisTrack()
#ases_bam <- system.file("/external/rprshnas01/netdata_kcni/stlab/RamseyMielnik/STAR_results/coord_bams/343_1Aligned.sortedByCoord.out.bam", package="Gviz") #maybe useless?
L23_track <- AlignmentsTrack("/external/rprshnas01/netdata_kcni/stlab/RamseyMielnik/STAR_results/coord_bams/343_1Aligned.sortedByCoord.out.bam", type=c("coverage"), legend = TRUE, ylim=c(0,102), fill="green", showTitle=FALSE, background.title="white", col.axis="black")
GABA_track <- AlignmentsTrack("/external/rprshnas01/netdata_kcni/stlab/RamseyMielnik/STAR_results/coord_bams/348_1Aligned.sortedByCoord.out.bam", type=c("coverage"), legend = TRUE, ylim=c(0,102), fill="magenta", showTitle=FALSE, background.title="white", col.axis="black")
L6_track <- AlignmentsTrack("/external/rprshnas01/netdata_kcni/stlab/RamseyMielnik/STAR_results/coord_bams/352_1Aligned.sortedByCoord.out.bam", type=c("coverage"), legend = TRUE, ylim=c(0,102), fill="green", showTitle=FALSE, background.title="white", col.axis="black")
L5_track <- AlignmentsTrack("/external/rprshnas01/netdata_kcni/stlab/RamseyMielnik/STAR_results/coord_bams/356_1Aligned.sortedByCoord.out.bam", type=c("coverage"), legend = TRUE, ylim=c(0,102), fill="green", showTitle=FALSE, background.title="white", col.axis="black")
L4_track <- AlignmentsTrack("/external/rprshnas01/netdata_kcni/stlab/RamseyMielnik/STAR_results/coord_bams/363_1Aligned.sortedByCoord.out.bam", type=c("coverage"), legend = TRUE, ylim=c(0,102), fill="green", showTitle=FALSE, background.title="white", col.axis="black")
#ases_track <- DataTrack(range="/external/rprshnas01/netdata_kcni/stlab/RamseyMielnik/STAR_results/coord_bams/343_1Aligned.sortedByCoord.out.bam", genome="mm10", type="l", name="Coverage", window=-1, chromosome="2")
class(L4_track)

png(
  'GvizStack.png',
  width = 600,
  height = 1800)

plotTracks(list(grtrack, L4_track, L23_track, L5_track, L6_track, GABA_track, gtrack), from=25291181, to=25319187, chromosome = "2", sizes = c(1,2,2,2,2,2,0.33))

dev.off()

#export results to xlsx (this is before correcting for technical duplicates; do the new export in analysis_RamChrom script)

exported_df <- cov_combined_df[,c(1,3,6,5,7:9,24)]
colnames(exported_df)[1] <- "ShortSampleID"
colnames(exported_df)[6] <- "Grin1RawCount"
colnames(exported_df)[8] <- "Grin1Accessibility"

library(xlsx)
write.xlsx(exported_df, "/external/rprshnas01/kcni/ychen/git/RamChromatin/GenomeAccessData.xlsx")
