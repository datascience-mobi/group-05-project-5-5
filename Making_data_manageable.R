






###### First steps to load data and manage unhandy data #######
###############################################################

#Load sample data
Samples <- readRDS(file = "Mantle-Bcell_list.RDS.gz")

#Loading annotation data file
input_data_csv <-
  read.csv(file = "sample_annotation.csv", sep = ",")
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
hist(
  log10(mean_cancer_coverage),
  breaks = "fd",
  main = "Cancer coverage: Mean frequency",
  xlab = "Common logarithm of coverages",
  col = c("black"),
  border = "white"
)
abline(v = log10(quantile(
  mean_cancer_coverage,
  probs = seq(0, 1, 0.1),
  na.rm = TRUE
)),
col = colors(256),
lwd = 2)


mean_healthy_coverage <- rowMeans(healthy_coverage)
hist(
  log10(mean_healthy_coverage),
  breaks = "fd",
  main = "Healthy coverage: Mean frequency",
  xlab = "Common logarithm of coverages",
  col = c("black"),
  border = "white"
)
abline(v = log10(quantile(
  mean_healthy_coverage,
  probs = seq(0, 1, 0.1),
  na.rm = TRUE
)),
col = colors(256),
lwd = 2)

sd_cancer_coverage <- apply(cancer_coverage, 1, sd)
hist(log10(sd_cancer_coverage), breaks = "fd")

sd_healthy_coverage <- apply(healthy_coverage, 1, sd)
hist(log10(sd_healthy_coverage), breaks = "fd")

#include mean and sd column to cancer and healthy data set

#cancer_coverage <-  cbind(cancer_coverage, mean_cancer_coverage)
#healthy_coverage <- cbind(healthy_coverage, mean_healthy_coverage)

#cancer_coverage <- cbind(cancer_coverage, sd_cancer_coverage)
#healthy_coverage <- cbind(healthy_coverage, sd_healthy_coverage)

####find coverage value for threshold and remove coverages in threshold --> don´t loose more than 90% of information

sum(cancer_coverage == 0)

#cancer coverages: lower boundary
threshold <-
  quantile(mean_cancer_coverage,
           probs = seq(0.05, 0.05),
           na.rm = TRUE)

#cancer coverages: upper boundary
threshold2 <-
  quantile(mean_cancer_coverage,
           probs = seq(0.999, 0.999),
           na.rm = TRUE)

#healthy coverages: lower boundary
threshold3 <-
  quantile(mean_healthy_coverage,
           probs = seq(0.05, 0.05),
           na.rm = TRUE)

#healthy coverages: upper boundary
threshold4 <-
  quantile(mean_healthy_coverage,
           probs = seq(0.999, 0.999),
           na.rm = TRUE)

##nestled for loops to set every value of cancer coverage and cancer beta value to NA if they are in threshold
for (i in 1:nrow(cancer_coverage)) {
  for (j in 1:ncol(cancer_coverage)) {
    if (cancer_coverage[i, j] <= threshold) {
      cancer_coverage[i, j] <- 0
    }
    
    if (cancer_coverage[i, j] >= threshold2) {
      cancer_coverage[i, j] <- 0
    }
    
    if (cancer_coverage[i, j] == 0) {
      cancer_coverage[i, j] <- NA
      cancer_beta_values[i, j] <- NA
      
    }
  }
}

##nested for loops to set every value of healthy coverage and cancer beta value to NA if they are in threshold
for (i in 1:nrow(healthy_coverage)) {
  for (j in 1:ncol(healthy_coverage)) {
    if (healthy_coverage[i, j] <= threshold3) {
      healthy_coverage[i, j] <- 0
    }
    
    if (healthy_coverage[i, j] >= threshold4) {
      healthy_coverage[i, j] <- 0
    }
    
    if (healthy_coverage[i, j] == 0) {
      healthy_coverage[i, j] <- NA
      healthy_beta_values[i, j] <- NA
      
    }
  }
}

#overview after data clean up

#cancer
mean_cancer_coverage <- rowMeans(cancer_coverage)
hist(
  log10(mean_cancer_coverage),
  breaks = "fd",
  main = "Cancer coverage: Mean frequency",
  xlab = "Common logarithm of coverages",
  col = c("black"),
  border = "white"
)
abline(v = log10(quantile(
  mean_cancer_coverage,
  probs = seq(0, 1, 0.1),
  na.rm = TRUE
)),
col = colors(256),
lwd = 2)

#healthy
mean_healthy_coverage <- rowMeans(healthy_coverage)
hist(
  log10(mean_healthy_coverage),
  breaks = "fd",
  main = "Healthy coverage: Mean frequency",
  xlab = "Common logarithm of coverages",
  col = c("black"),
  border = "white"
)
abline(v = log10(quantile(
  mean_healthy_coverage,
  probs = seq(0, 1, 0.1),
  na.rm = TRUE
)),
col = colors(256),
lwd = 2)

####

#deal with NA's
##output of the number of all NA's in healthy genes
sum(is.na(healthy_beta_values))
sum(is.na(cancer_beta_values))

##add a new column with the number of NA's per gene
healthy_beta_values$Number_of_NA <-
  rowSums(is.na(healthy_beta_values))
cancer_beta_values$Number_of_NA <-
  rowSums(is.na(cancer_beta_values))

##Histogram of NA's
###maybe we can work out the scales :/
hist(
  healthy_beta_values$Number_of_NA,
  main =  "NA's per Gene in Healthy Samples",
  breaks = 5,
  xlab = "Number of NA's",
  ylab = "Number of Genes",
  col = "seagreen2"
)
hist(
  cancer_beta_values$Number_of_NA,
  main =  "NA's per Gene in MCL Samples",
  breaks = 5,
  xlab = "Number of NA's",
  ylab = "Number of Genes",
  col = "indianred2"
)

#set a threshold for the NA values and remove the gene if there are too much NA's
healthy_beta_values <-
  healthy_beta_values[!(healthy_beta_values$`Number_of_NA` > 2), ]
cancer_beta_values <-
  cancer_beta_values[!(cancer_beta_values$`Number_of_NA` > 2), ]

cancer_beta_values <-
  healthy_beta_values[!(healthy_beta_values$`Number_of_NA` > 2), ]
healthy_beta_values <-
  cancer_beta_values[!(cancer_beta_values$`Number_of_NA` > 2), ]

##remove column Number_of_NA
healthy_beta_values <-
  healthy_beta_values[, -which(colnames(healthy_beta_values)  %in%  c('Number_of_NA'))]
cancer_beta_values <-
  cancer_beta_values[, -which(colnames(cancer_beta_values) %in% c('Number_of_NA'))]

#replace remaining NA's with the mean of the respective gene
#first transposing the data frame because working on columns, e.g getting the mean, is easier than with rows
transposed_healthy_beta_values <- t(healthy_beta_values)
transposed_cancer_beta_values <- t(cancer_beta_values)

#how many elements do we have in our new data frame?
dim(transposed_healthy_beta_values)
dim(transposed_cancer_beta_values)

#going through all elements of the (already reduced) data frame and replace NA's with mean
for (i in 1:ncol(transposed_healthy_beta_values)) {
  transposed_healthy_beta_values[is.na(transposed_healthy_beta_values[, i]), i] <-
    mean(transposed_healthy_beta_values[, i], na.rm = TRUE)
}

for (i in 1:ncol(transposed_cancer_beta_values)) {
  transposed_cancer_beta_values[is.na(transposed_cancer_beta_values[, i]), i] <-
    mean(transposed_cancer_beta_values[, i], na.rm = TRUE)
}

#did we eliminate all NA's?
sum(is.na(transposed_healthy_beta_values))
sum(is.na(transposed_cancer_beta_values))

#retranspose to data frame
healthy_beta_values <- data.frame(t(transposed_healthy_beta_values))
cancer_beta_values <- data.frame(t(transposed_cancer_beta_values))

#check if genes of one data frame are in the other data frame
sum(rownames(healthy_beta_values) == rownames(cancer_beta_values))

#are important genes still included?
important_genes <-
  data.frame(
    c(
      "ENSG00000176887",
      "ENSG00000185551",
      "ENSG00000141510",
      "ENSG00000110092",
      "ENSG00000106546",
      "ENSG00000169855",
      "ENSG00000125398",
      "ENSG00000078399",
      "ENSG00000039068",
      "ENSG00000081377",
      "ENSG00000054598",
      "ENSG00000123689",
      "ENSG00000211445",
      "ENSG00000131981",
      "ENSG00000172005",
      "ENSG00000106236",
      "ENSG00000007372",
      "ENSG00000105825",
      "ENSG00000159445",
      "ENSG00000122691"
    )
  )

##we only have to check one data set because cancer and healthy have the same genes
important_genes_in_data_set <- data.frame()
for (i in 1:nrow(important_genes)) {
  important_genes_in_data_set[i,] <-
    cancer_beta_values[which(row.names(cancer_beta_values) == b[i]),]
}

##look up if there are NAs --> would be bad
sum(is.na(important_genes_in_data_set))

##hurray no NAs





