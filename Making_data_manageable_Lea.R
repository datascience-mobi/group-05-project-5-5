

###### First steps to load data and manage unhandy data #######
###############################################################

#Load sample data
Samples <- readRDS(file = "Mantle-Bcell_list.RDS.gz")

#Loading annotation data file
input_data_csv <- read.csv(file = "sample_annotation.csv", sep = ",")
View(input_data_csv)

#creating a data frame to understand the different names of the original data frame of genes
Name <-
  c(
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Symbol",
    "entrezID",
    "CpG",
    "GC",
    "G",
    "C",
    ".bed",
    ".bed_coverage",
    "average depth of sequencing coverage",
    "beta value",
    "Whole genome bisulfite sequencing"
  )
Name2 <-
  c(
    "...",
    "Start position (bp)",
    "End position (bp",
    "...",
    "Corresponding gene symbol (only in gene table)",
    "Gene ID from entrezGene Database",
    "# of CpGs in region",
    "Number of GCs in region",
    "Number of Gs in region",
    "Number of Cs in region",
    "methylation BETA value (Percentage of Methylation)",
    "sequencing depth at this position aka. the number of reads mapped to this region",
    "can be defined theoretically as LN/G, where L is the read length, N is the number of reads and G is the haploid genome length.",
    " from 0 (demethylated) to 1 (fully methylated)",
    "https://de.wikipedia.org/wiki/Bisulfit-Sequenzierung"
  )
table_names <- data.frame(Name = Name, y = Name2)

#copying "genes" data from general list to create a data frame of genes
Gene_data_frame <- Samples$genes
dim(Gene_data_frame)
View(Gene_data_frame)

#some pre-cleaning up: deleting x and y chromosome specific genes

Gene_data_frame_x_y <-
  Gene_data_frame[-which(Gene_data_frame$Chromosome == "chrX"),]
Gene_data_frame_x_y <-
  Gene_data_frame_x_y[-which(Gene_data_frame_x_y$Chromosome == "chrY"),]

#tidy up the data by spliting up the data to different data frame
healthy_coverage <- Gene_data_frame_x_y[, 21:25]
cancer_coverage <- Gene_data_frame_x_y[, 26:30]
healthy_beta_values <- Gene_data_frame_x_y[, 11:15]
cancer_beta_values <- Gene_data_frame_x_y[, 16:20]

####Coverage and Beta-Value problem####
#######################################

#mean and sd value histogram of every gene

mean_cancer_coverage <- rowMeans(cancer_coverage)
hist(log10(mean_cancer_coverage), breaks = "fd")

mean_healthy_coverage <- rowMeans(healthy_coverage)
hist(log10(mean_healthy_coverage), breaks = "fd")

sd_cancer_coverage <- apply(cancer_coverage, 1, sd)
hist(log10(sd_cancer_coverage), breaks = "fd")

sd_healthy_coverage <- apply(healthy_coverage, 1, sd)
hist(log10(sd_healthy_coverage), breaks = "fd")

#include mean and sd column

cancer_coverage <-  cbind(cancer_coverage, mean_cancer_coverage)
healthy_coverage <- cbind(healthy_coverage, mean_healthy_coverage)

cancer_coverage <- cbind(cancer_coverage, sd_cancer_coverage)
healthy_coverage <- cbind(healthy_coverage, sd_healthy_coverage)

####find coverage value for threshold and remove coverages in threshold






















































#deal with NA's
##output of the number of all NA's in healthy genes
sum(is.na(healthy_beta_values))
sum(is.na(cancer_beta_values))

##add a new column with the number of NA's per gene
healthy_beta_values$new=rowSums(is.na(healthy_beta_values))
colnames(healthy_beta_values)[colnames(healthy_beta_values) == 'new'] <- 'Number_of_NA'
cancer_beta_values$new=rowSums(is.na(cancer_beta_values))
colnames(cancer_beta_values)[colnames(cancer_beta_values) == 'new'] <- 'Number_of_NA'

##Histogram of NA's
###Skala muss noch angepasst wereden
hist(healthy_beta_values$Number_of_NA, main = "NA's per Gene in Healthy Samples", breaks = 5, xlab = "Number of NA's", ylab = "Number of Genes", ylim = c(0,1000), col = "seagreen2" )
hist(cancer_beta_values$Number_of_NA, main = "NA's per Gene in MCL Samples", breaks = 5, xlab = "Number of NA's", ylab = "Number of Genes",ylim = c(0,1000) , col = "indianred2")

#set a threshold for the NA values and remove the gene if there are to much NA's
healthy_beta_values <- healthy_beta_values[!(healthy_beta_values$`Number_of_NA`>2), ]
cancer_beta_values <- cancer_beta_values[!(cancer_beta_values$`Number_of_NA`>2), ]

##remove column Number of NA
healthy_beta_values = healthy_beta_values[,  -which( colnames(healthy_beta_values)  %in%  c('Number_of_NA'))]
cancer_beta_values = cancer_beta_values[, -which( colnames(cancer_beta_values) %in% c('Number_of_NA'))]

#add column with mean of each row
healthy_beta_values$new=rowMeans(healthy_beta_values, na.rm = TRUE)
colnames(healthy_beta_values)[colnames(healthy_beta_values) == 'new'] <- 'mean_value'
cancer_beta_values$new=rowMeans(cancer_beta_values, na.rm = TRUE)
colnames(cancer_beta_values)[colnames(cancer_beta_values) == 'new'] <- 'mean_value'

#replace remaining NA's with the mean of the respective gene
...



#remove column mean-value
healthy_beta_values = healthy_beta_values[,  -which( colnames(healthy_beta_values)  %in%  c('mean_value'))]
cancer_beta_values = cancer_beta_values[, -which( colnames(cancer_beta_values) %in% c('mean_value'))]

#test if there are realy no more NA's
sum(is.na(healthy_beta_values))
sum(is.na(cancer_beta_values))

#compare the 2 data frames and remove rows which are only in one of them
 () {
  rownames(healthy_beta_values) != rownames(cancer_beta_values)
  healthy_beta_values <- healthy_beta_values
  cancer_beta_values <- cancer_beta_values
}  