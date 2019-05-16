#creating a variable which you can call for looking up the data in console
input_data <- readRDS(file = "Mantle-Bcell_list.RDS.gz")


#looking up structure of data
names(input_data)
View(input_data)

#copying data from general list to create a data frame of genes
Gene_data_frame <- input_data$genes
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

##rename columns
rownames(coverage_beta_values_data_frame_manipulated) = c(1:53620)

##look up differences
dim(coverage_beta_values_data_frame_manipulated)
dim(coverage_beta_values_data_frame)

#visualize data for better understanding

plot(density(coverage_beta_values_data_frame_manipulated$Bcell_naive_VB_NBC_NC11_41.bed))

pairs(coverage_beta_values_data_frame_manipulated[1:10, 1:ncol(coverage_beta_values_data_frame_manipulated)], pch = 20, cex = 0.2, col = "grey")

b <- data.matrix(coverage_beta_values_data_frame_manipulated)

heatmap(cor(b[1:5000, 1:10]), col = cm.colors(256))

qqmath(coverage_beta_values_data_frame_manipulated[,1],coverage_beta_values_data_frame_manipulated)

##trying to understand how to make a volcano plot

      with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
      
      # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
      with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
      with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
      with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
      
      # Label points with the textxy function from the calibrate plot
      library(calibrate)
      with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))
        



plot(table(coverage_beta_values_data_frame_manipulated$Bcell_naive_VB_NBC_NC11_41.bed_coverage))


with(coverage_beta_values_data_frame_manipulated, plot(-log10(coverage_beta_values_data_frame_manipulated$Bcell_naive_VB_NBC_NC11_41.bed), -log10(coverage_beta_values_data_frame_manipulated$Bcell_naive_VB_NBC_NC11_41.bed_coverage), pch=20, main="Volcano plot", xlim=c(-2,5), col = "cornflowerblue"))

##clustering: https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
km <- kmeans(coverage_beta_values_data_frame_manipulated$Bcell_naive_VB_NBC_NC11_41.bed, 5, nstart = 10)

wss = sapply(2:7, function(k) {
  kmeans(coverage_beta_values_data_frame_manipulated$Bcell_naive_VB_NBC_NC11_41.bed, centers = k)$tot.withinss
})
plot(2:7, wss, type = "b", pch = 19, xlab = "Number of clusters K", ylab = "Total within-clusters sum of squares")



rownames(coverage_beta_values_data_frame_manipulated) = c(1:53620)

dim(coverage_beta_values_data_frame_manipulated)





##scaling
coverage_beta_values_data_frame_manipulated2 = data.frame(scale(coverage_beta_values_data_frame_manipulated))




hist(log10(coverage_beta_values_data_frame_manipulated[Bcell_naive_VB_S001JP51.bed_coverage[ ,1:4]]) )

hist(input_data$tiling$Bcell_naive_VB_NBC_NC11_41.bed, main = "Normal B cell Tiling window", xlab = "β-Value DNA Methylation", ylab = "Frequency", col = "green")
hist(input_data$promoters$Bcell_naive_VB_NBC_NC11_41.bed, main = "Normal B cell Promoters", xlab = "β-Value DNA Methylation", ylab = "Frequency", col = "green")
hist(input_data$genes$Bcell_naive_VB_NBC_NC11_41.bed, main = "Normal B cell Genes", xlab = "β-Value DNA Methylation", ylab = "Frequency", col = "green")
hist(input_data$cpgislands$Bcell_naive_VB_NBC_NC11_41.bed, main = "Normal B cell CpG Islands", xlab = "β-Value DNA Methylation", ylab = "Frequency", col = "green")

hist(input_data$tiling$cancer_VB_S01FJZA1.bed, main = "Cancer cell Tiling window", xlab = "β-Value DNA Methylation", ylab = "Frequency", col = "red")
hist(input_data$promoters$cancer_VB_S01FJZA1.bed, main = "Cancer cell Promoters", xlab = "β-Value DNA Methylation", ylab = "Frequency", col = "red")
hist(input_data$genes$cancer_VB_S01FJZA1.bed, main = "Cancer cell Genes", xlab = "β-Value DNA Methylation", ylab = "Frequency", col = "red" )
hist(input_data$cpgislands$cancer_VB_S01FJZA1.bed, main = "Cancer cell CpG Islands", xlab = "β-Value DNA Methylation", ylab = "Frequency", col = "red")


mean <- rowMeans(coverage_beta_values_data_frame_manipulated[, 11:15])
sd_row <- sd(coverage_beta_values_data_frame_manipulated[,11:15])
hist(log10(mean), breaks = "scott")
summary(log10(mean))
    