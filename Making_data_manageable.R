
###### First steps to load data and manage unhandy data #######
###############################################################

#Load sample data
Samples <- readRDS(file = "Mantle-Bcell_list.RDS.gz")

#Loading annotation data file
input_data_csv <- read.csv(file ="sample_annotation.csv", sep = ",")
View(input_data_csv)

#creating a data frame to understand the different names of the original data frame of genes
Name <- c("Chromosome","Start", "End", "Strand", "Symbol", "entrezID", "CpG", "GC", "G", "C", ".bed", ".bed_coverage", "average depth of sequencing coverage", "beta value", "Whole genome bisulfite sequencing")
Name2 <- c("...", "Start position (bp)", "End position (bp", "...", "Corresponding gene symbol (only in gene table)", "Gene ID from entrezGene Database", "# of CpGs in region", "Number of GCs in region", "Number of Gs in region", "Number of Cs in region", "methylation BETA value (Percentage of Methylation)", "sequencing depth at this position aka. the number of reads mapped to this region", "can be defined theoretically as LN/G, where L is the read length, N is the number of reads and G is the haploid genome length.", " from 0 (demethylated) to 1 (fully methylated)", "https://de.wikipedia.org/wiki/Bisulfit-Sequenzierung")
table_names <- data.frame(Name = Name, y= Name2)

#copying "genes" data from general list to create a data frame of genes
Gene_data_frame <- Samples$genes
dim(Gene_data_frame)
View(Gene_data_frame)

#some pre-cleaning up: deleting x and y chromosome specific genes

Gene_data_frame_x_y <-  Gene_data_frame[-which(Gene_data_frame$Chromosome == "chrY") , ]  
Gene_data_frame_x_y <-  Gene_data_frame_x_y[-which(Gene_data_frame_x_y$Chromosome == "chrX") , ]  

View(Gene_data_frame_x_y)
#tidy up the data by spliting up the data to different data frame 
healthy_coverage <- Gene_data_frame_x_y[, 26:30]
cancer_coverage <- Gene_data_frame_x_y[, 21:25]
healthy_beta_values <- Gene_data_frame_x_y[, 11:15]
cancer_beta_values <- Gene_data_frame_x_y[, 16:20]

####Coverage and Beta-Value problem####
#######################################

#plot of mean coverages
mean <- rowMeans(cancer_coverage)
hist(log10(mean), breaks = "fd")

#plot of sd coverages

#finding coverages in threshold ##does not work :/



##mit extraspalte und dann if 



















































