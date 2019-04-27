#creating a variable which you can call for looking up the data in console
input_data <- readRDS(file = "Mantle-Bcell_list.RDS.gz")


#looking up structure of data
names(input_data)
View(input_data)

#copying data from general list to create a data frame of genes
Gene_data_frame<- input_data$genes
dim(Gene_data_frame)
View(Gene_data_frame)

#creating a data frame to understand the different names of the original data frame of genes
Name <- c("Chromosome","Start", "End", "Strand", "Symbol", "entrezID", "CpG", "GC", "G", "C", ".bed", ".bed_coverage", "average depth of sequencing coverage", "beta value", "Whole genome bisulfite sequencing")
Name2 <- c("...", "Start position (bp)", "End position (bp", "...", "Corresponding gene symbol (only in gene table)", "Gene ID from entrezGene Database", "# of CpGs in region", "Number of GCs in region", "Number of Gs in region", "Number of Cs in region", "methylation BETA value (Percentage of Methylation)", "sequencing depth at this position aka. the number of reads mapped to this region", "can be defined theoretically as LN/G, where L is the read length, N is the number of reads and G is the haploid genome length.", " from 0 (demethylated) to 1 (fully methylated)", "https://de.wikipedia.org/wiki/Bisulfit-Sequenzierung")
table_names <- data.frame(Name = Name, y= Name2)

View(table_names)

#overview of another data file
input_data_csv <- read.csv(file ="sample_annotation.csv", sep = ",")
View(input_data_csv)

#tidy up the data
##split up the data to data frame with beta values and covarage
coverage_beta_values_data_frame <- Gene_data_frame[, 11:30]
View(coverage_beta_values_data_frame)

##remove columns with missing values
rmv.rows = apply(coverage_beta_values_data_frame, 1, function(x) {
  sum(is.na(x))
})  # Go through each row and sum up all missing values
coverage_beta_values_data_frame_manipulated = coverage_beta_values_data_frame[-which(rmv.rows > 0), ]  # Removing any row with 1 or more missing values
rm(rmv.rows)
dim(coverage_beta_values_data_frame_manipulated)
dim(coverage_beta_values_data_frame)

#visualize data for better understanding

plot(density(coverage_beta_values_data_frame_manipulated$Bcell_naive_VB_NBC_NC11_41.bed))

pairs(coverage_beta_values_data_frame_manipulated[1:10, 1:ncol(coverage_beta_values_data_frame_manipulated)], pch = 20, cex = 0.2, col = "grey")

heatmap(b[1:50, 1:ncol(coverage_beta_values_data_frame_manipulated)], col = cm.colors(256))

b <- data.matrix(coverage_beta_values_data_frame_manipulated)











