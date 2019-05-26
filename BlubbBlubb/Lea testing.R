#creating a variable which you can call for looking up the data in console
input_data <- readRDS(file = "Mantle-Bcell_list.RDS.gz")

#looking up structure of data
names(input_data)
View(input_data)


#create a new data file wich contains the gene Values
gene_data_frame<- input_data$genes
dim(gene_data_frame)
View(gene_data_frame)

#number and names of the columns before clearing
ncol(healthy_data_frame)
colnames(gene_data_frame)


#create a new data file with net necessary data of all healthy patioents
##bad way to do this
###healthy_data_frame = gene_data_frame[,  -which(colnames(gene_data_frame)  %in%  c("Chromosome", "Start", "End", "Strand", "symbol", "entrezID", "CpG", "GC", "C", "G", "cancer_VB_S01FE8A1.bed", "cancer_VB_S01FF6A1.bed", "cancer_VB_S01FH2A1.bed", "cancer_VB_S01FJZA1.bed", "cancer_VB_S01FKXA1.bed", "cancer_VB_S01FE8A1.bed_coverage", "cancer_VB_S01FF6A1.bed_coverage", "cancer_VB_S01FH2A1.bed_coverage", "cancer_VB_S01FJZA1.bed_coverage", "cancer_VB_S01FKXA1.bed_coverage"))]
##good way to do this
healthy_coverage <- gene_data_frame[, 21:25]
healthy_beta_values <- gene_data_frame[, 11:15]
cancer_coverage <- gene_data_frame[, 26:30]
cancer_beta_values <- gene_data_frame[, 16:20]


#rename the columns of healthy samples


#check the new data frame of healthy samples
ncol(healthy_data_frame)
colnames(healthy_data_frame)
View (healthy_data_frame)







#deal with NA's
##output of the number of all NA's in healthy genes
sum(is.na(healthy_data_frame))

##output of the number of NA's per gene in healthy genes
rowSums(is.na(healthy_data_frame))

##add a new column with the number of NA's per gene
healthy_beta_values$new=rowSums(is.na(healthy_data_frame))
colnames(healthy_beta_values)[colnames(healthy_beta_values) == 'new'] <- 'Number of NA'

###remove a new row
####healthy_beta_values <- healthy_beta_values[, -which(colnames(healthy_beta_values) %in% c("newcolumn"))]


#set a threshold for the NA values and remove the gene if there are to much NA's
if ("Number of NA" > 2){
  
  healthy_beta_values <- healthy_beta_values[c(x), ]
  
}




