

###### First steps to load data and manage unhandy data #######
###############################################################

#Load sample data
Samples <- readRDS(file = "Mantle-Bcell_list.RDS.gz")

#Loading annotation data file
input_data_csv <- read.csv(file = "sample_annotation.csv", sep = ",")

#Franzi's code for importing sample_annotation.csv: 
#library(readr)
#sample_annotation <- read_csv("Uni/4. Semester/Bioinfo/MantlevsBcell/sample_annotation.csv")
#View(sample_annotation)
#Thus extra code for changing data name for future working on the dataset.
input_data_csv <- sample_annotation

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

#tidy up the data by splitting up the data to different data frame
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

##add a new column with the number of NA's per gene
healthy_beta_values$new=rowSums(is.na(healthy_beta_values))
colnames(healthy_beta_values)[colnames(healthy_beta_values) == 'new'] <- 'Number_of_NA'

##Histogram of NA's
###Skala muss noch angepasst wereden
hist(healthy_beta_values$Number_of_NA, breaks =5, xlab = "Number of NA's", ylab = "per Genes", ylim = c(0,55000), col = "seagreen2" )

#set a threshold for the NA values and remove the gene if there are too many (3) NA's 
healthy_beta_values<- healthy_beta_values[!(healthy_beta_values$`Number_of_NA`>2), ]

##remove column Number of NA
healthy_beta_values = healthy_beta_values[,  -which( colnames(healthy_beta_values)  %in%  c('Number_of_NA'))]

#add column with mean of each row
healthy_beta_values$new=rowMeans(healthy_beta_values, na.rm = TRUE)
colnames(healthy_beta_values)[colnames(healthy_beta_values) == 'new'] <- 'mean_value'

#how many na do we have now?
sum(is.na(healthy_beta_values))

#replace remaining NA's with the mean of the respective gene
#first transposing the data frame because working on columns, e.g getting the mean, is easier than with rows
transposed_healthy_beta_values <- t(healthy_beta_values)

#unneccessary step
transposed_healthy_beta_values [which(is.na(transposed_healthy_beta_values))]

#how many elements do we have in our new data frame?
dim(transposed_healthy_beta_values)

#going through all elements of the (already reduced) data frame and replace NA's with mean
for (i in 1:52814){
  transposed_healthy_beta_values[is.na(transposed_healthy_beta_values[,i]), i] <- mean(transposed_healthy_beta_values[,i], na.rm = TRUE) 
}

#did we eliminate all NA's?
sum(is.na(transposed_healthy_beta_values))

#retranspose the data frame
healthy_beta_values2 <- t(transposed_healthy_beta_values)

#check if the NA's are successfully replaced with the mean with histogram
hist(healthy_beta_values$mean_value, breaks = "fd")
hist(healthy_beta_values2$mean_value, breaks = "fd")

#remove column mean-value
healthy_beta_values = healthy_beta_values[,  -which( colnames(healthy_beta_values)  %in%  c('mean_value'))]

#test if there are really no more NA's
sum(is.na(healthy_beta_values))