

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

#mean and sd value histogram of every gene + quantiles

mean_cancer_coverage <- rowMeans(cancer_coverage)
hist(log10(mean_cancer_coverage), breaks = "fd")
abline(v = log10(quantile(mean_cancer_coverage, probs = seq(0, 1, 0.1), na.rm = TRUE)))

mean_healthy_coverage <- rowMeans(healthy_coverage)
hist(log10(mean_healthy_coverage), breaks = "fd")
abline(v = log10(quantile(mean_healthy_coverage, probs = seq(0, 1, 0.1), na.rm = TRUE)))

sd_cancer_coverage <- apply(cancer_coverage, 1, sd)
hist(log10(sd_cancer_coverage), breaks = "fd")

sd_healthy_coverage <- apply(healthy_coverage, 1, sd)
hist(log10(sd_healthy_coverage), breaks = "fd")

#include mean and sd column to cancer and healthy data set

cancer_coverage <-  cbind(cancer_coverage, mean_cancer_coverage)
healthy_coverage <- cbind(healthy_coverage, mean_healthy_coverage)

cancer_coverage <- cbind(cancer_coverage, sd_cancer_coverage)
healthy_coverage <- cbind(healthy_coverage, sd_healthy_coverage)

####find coverage value for threshold and remove coverages in threshold

sum(cancer_coverage == 0)

threshold <- quantile(mean_cancer_coverage, probs = seq(0.10, 0.10, 0.05), na.rm = TRUE)
threshold2 <- quantile(mean_cancer_coverage, probs = seq(0.90, 0.90, 0.05), na.rm = TRUE)

threshold3 <- quantile(mean_healthy_coverage, probs = seq(0.10, 0.10, 0.05), na.rm = TRUE)
threshold4 <- quantile(mean_healthy_coverage, probs = seq(0.90, 0.90, 0.05), na.rm = TRUE)

##nestled for loops to set every value of cancer coverage in threshold to 0
for (i in 1: 53470 ) {
  
  for(j in 1:5 ) {
    
    if(cancer_coverage[i,j] <= threshold) {
      
      cancer_coverage[i,j] <- 0
    }
    
    if(cancer_coverage[i,j] >= threshold2){
      
      cancer_coverage[i,j] <- 0
    }
    
    if(cancer_coverage[i,j] == 0) {
      
      cancer_coverage[i,j] <- NA
      cancer_beta_values[i,j] <- NA
        
    } 
  }
}

for (i in 1: 53470 ) {
  
  for(j in 1:5 ) {
    
    if(healthy_coverage[i,j] <= threshold3) {
      
      healthy_coverage[i,j] <- 0
    }
    
    if(healthy_coverage[i,j] >= threshold4){
      
      healthy_coverage[i,j] <- 0
    }
    
    if(healthy_coverage[i,j] == 0) {
      
      healthy_coverage[i,j] <- NA
      healthy_beta_values[i,j] <- NA
      
    } 
  }
}
#overview after data cleaning
mean_cancer_coverage <- rowMeans(cancer_coverage)
hist(log10(mean_cancer_coverage), breaks = "fd")

mean_healthy_coverage <- rowMeans(healthy_coverage)
hist(log10(mean_healthy_coverage), breaks = "fd")


rmv.rows <-  apply(cancer_beta_values, 1, function(x) {
  sum(is.na(x))
})  # Go through each row and sum up all missing values
hist(rmv.rows) # The rows where there is atleast 1 missing value
sum(rmv.rows > 2) #Rows which will be removed because of the threshold 3 NA or more


rmv.rows2 <-  apply(healthy_beta_values, 1, function(x) {
  sum(is.na(x))
})  # Go through each row and sum up all missing values
hist(rmv.rows2) # The rows where there is atleast 1 missing value
sum(rmv.rows2 > 2) #Rows which will be removed because of the threshold 3 NA or more


