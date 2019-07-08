#----data preparation for pca-----------------------------------------------------------------------------------------

## replacing 0 and 1 in beta values with approximate values, thus no -inf and +inf values in m values

cancer_beta_values[cancer_beta_values == 0] <- 0.00001
cancer_beta_values[cancer_beta_values == 1] <- 0.99999
healthy_beta_values[healthy_beta_values == 0] <- 0.00001
healthy_beta_values[healthy_beta_values == 1] <- 0.99999

# transforming beta values to M values and creating a separate dataframe with those values

cancer_m_values <-
  data.frame(log2(cancer_beta_values / (1 - cancer_beta_values)))
healthy_m_values <-
  data.frame(log2(healthy_beta_values / (1 - healthy_beta_values)))

# changing the ending of patients names from .bed to .M for better overview

# changing healthy patients names
colnames(healthy_m_values) <-
  c(
    "Bcell_naive_VB_NBC_NC11_41.M",
    "Bcell_naive_VB_NBC_NC11_83.M",
    "Bcell_naive_VB_S001JP51.M",
    "Bcell_naive_VB_S00DM851.M",
    "Bcell_naive_VB_S01ECGA1.M"
  )

# changing cancer patients names

colnames(cancer_m_values) <- c(
  "cancer_VB_S01FE8A1.M",
  "cancer_VB_S01FF6A1.M",
  "cancer_VB_S01FH2A1.M",
  "cancer_VB_S01FJZA1.M",
  "cancer_VB_S01FKXA1.M"
)


#data reduction with PCA

#merging both m value dataframes (healthy and cancer) into one again for PCA
m_values <- cbind(healthy_m_values, cancer_m_values)
View(m_values)

#Apply PCA on data frame "m_values" with all m values. For that, the matrix needs to be transposed first
#(variables are scaled to have i) standard deviation one and ii) mean zero)
pca_m_values <- prcomp(t(m_values))
summary(pca_m_values)

#visualize pca variation

##choosing number of pc's to work with for batch effect detection (elbow method)

## calculating the variance of each principal component (sdev^2), then calculating the proportion of each
# variance by dividing it with the sum of the variacnes


std_dev <- pca_m_values$sdev
variance <- std_dev ^ 2
prop_var <- variance / sum(variance)

plot(
  prop_var,
  main = "Variance explained by principal components",
  xlab = "Principal Components",
  ylab = "Proportion of Variance Explained",
  type = "b"
)

#visualize pca
#plot(complete_m_values.pca$x[, 1], complete_m_values.pca$x[, 2])


#adding an extra column with the category of sample with which we can color the pc dots in a ggplot according to their sample group
pcs_of_m_values <-
  data.frame(cbind(
    pca_m_values$x,
    Samples = c(
      "Healthy",
      "Healthy",
      "Healthy",
      "Healthy",
      "Healthy",
      "Cancer",
      "Cancer",
      "Cancer",
      "Cancer",
      "Cancer"
    )
  ))

## remove the numbers in ggplot
#generate a ggplot / scatterplot to visualize the Sample points in a coordinate system with x -
#axis = PC1 and y - axis = PC2
install.packages("plotly")
library(plotly)
p <-
  ggplot(pcs_of_m_values, aes(x = PC1, y = PC2, group = Samples)) +
  geom_point (aes(shape = Samples, color = Samples), size = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )
p <- p + scale_colour_manual(values = c("indianred2", "seagreen2"))
p <- ggplotly(p)
p

#find out how much clusters do we need to group samples (obviously 2 would be perfect because healthy/cancer)


#--------------clustering with all genes, not only the significant ones-----------------------------------------------

wss <-  sapply(1:5, function(k) {
  kmeans(x = pca_m_values$x,
         centers = k,
         iter.max = 100)$tot.withinss
})
plot(
  1:5,
  wss,
  type = "b",
  pch = 19,
  xlab = "Number of clusters K",
  ylab = "Total within-clusters sum of squares"
)

##--> it seems we need 2 clusters (kink in the curve)

#find out if healthy and cancer samples are seperated
#variable with two center positions of value of rotated data
centers <-
  kmeans(
    x = t(pcs_of_m_values[1:10]),
    centers = 2,
    iter.max = 100
  )$centers

#testLea Angabe welcher Punkt zu welchem cluster gehört
cluster <-
  kmeans(
    x = t(pcs_of_m_values[1:10]),
    centers = 2,
    iter.max = 100
  )$cluster

##adding an extra column with the category of sample with which we can color the pc dots in a ggplot according to their sample group
centers <-
  data.frame(cbind(
    t(centers),
    Samples = c(
      "Healthy",
      "Healthy",
      "Healthy",
      "Healthy",
      "Healthy",
      "Cancer",
      "Cancer",
      "Cancer",
      "Cancer",
      "Cancer"
    )
  ))

#visualize cluster x1 and x2 and how samples are seperated
#testLea die punkte der beiden cluster unterschiedliche Formen
PC1 <- centers$X1
PC2 <- centers$X2
p_cluster <-
  ggplot(centers, aes(x = PC1, y = PC2, group = Samples)) +
  geom_point (aes(shape = as.factor(cluster), color = Samples), size = 4) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )
p_cluster <- p_cluster + scale_colour_manual(values = c("indianred2", "seagreen2"))
p_cluster <- p_cluster + scale_shape_manual(values = c(16, 17))
p_cluster <- ggplotly(p_cluster)
p_cluster



##############################################################################################################

## wilkoxon, kruskal wallis, Pearson correlation coefficient berechnen und ein permutation test

#changing the class of the columns PC1, PC2 and PC3 of pcs_of_m_values from "factor" to "numeric"
#to be able to work with the elements in the statistical tests
#pcs_of_m_values$PC1 <- as.numeric(as.character(pcs_of_m_values$PC1))
#pcs_of_m_values$PC2 <- as.numeric(as.character(pcs_of_m_values$PC2))
#pcs_of_m_values$PC3 <- as.numeric(as.character(pcs_of_m_values$PC3))

#New matrix with the PC values from 1 to 3 of the samples + categories of interest from the sample annotation
batch_pcs <-
  cbind(pcs_of_m_values[, 1:3], c(input_data_csv[, c(
    "BIOMATERIAL_PROVIDER",
    "BIOMATERIAL_TYPE",
    "cellTypeShort",
    "DONOR_SEX",
    "DISEASE",
    "FIRST_SUBMISSION_DATE",
    "SEQ_RUNS_COUNT",
    "DONOR_AGE"
  )]))

#making values of PC1-3 numeric for the tests

batch_pcs <- within(batch_pcs, {
  PC1 <- as.numeric(as.character(PC1))
})

batch_pcs <- within(batch_pcs, {
  PC2 <- as.numeric(as.character(PC2))
})

batch_pcs <- within(batch_pcs, {
  PC3 <- as.numeric(as.character(PC3))
})



#------p-values---------------------------------------------------------
p_values_matrix <- matrix(nrow = 8, ncol = 3)
colnames(p_values_matrix) <- c("PC1", "PC2", "PC3")
rownames(p_values_matrix) <-
  c(
    "BIOMATERIAL_PROVIDER",
    "BIOMATERIAL_TYPE",
    "DISEASE",
    "DONOR_SEX",
    "SEQ_RUNS_COUNT",
    "DONOR_AGE",
    "SUBMISSION_DATE",
    "CELL_TYPE"
  )

# define function for permutation test before the for-loop, because it is not necessary to define it every time the loop runs
cor.perm <- function (x, y, nperm = 1000)
{
  r.obs <- cor (x = x, y = y)
  p_value <- cor.test (x = x, y = y)$p.value
  #  r.per <- replicate (nperm, expr = cor (x = x, y = sample (y)))
  r.per <-
    sapply (
      1:nperm,
      FUN = function (i)
        cor (x = x, y = sample (y))
    )
  r.per <- c(r.per, r.obs)
  hist (r.per, xlim = c(-1, 1))
  abline (v = r.obs, col = 'red')
  P.per <- sum (abs (r.per) >= abs (r.obs)) / (nperm + 1)
  return (list (
    r.obs = r.obs,
    p_value = p_value,
    P.per = P.per
  ))
}

for (i in 1:ncol(p_values_matrix)) {
  
  # wilcoxon test: biomaterial provider
  ##H0: x and y differ by a location shift of mu=0
  ##alternative (two-sided/greater/less) and exact (value of p value) didn't need to be set
  ##(because are set like our settings by default. But it helps understanding under which conditions the test was made.
  bio_prov_test <- wilcox.test(
    batch_pcs [, i] ~ batch_pcs$BIOMATERIAL_PROVIDER,
    mu = 0,
    alt = "two.sided",
    conf.int = T,
    conf.level = 0.99,
    paired = F,
    exact = T
  )
  p_values_matrix  [1, i] <- as.matrix(bio_prov_test$p.value)
  
  # wilcoxon test: biomaterial type
  bio_type_test <- wilcox.test(
    batch_pcs [, i] ~ batch_pcs$BIOMATERIAL_TYPE,
    mu = 0,
    alt = "two.sided",
    conf.int = T,
    conf.level = 0.99,
    paired = F,
    exact = T
  )
  p_values_matrix  [2, i] <- as.matrix(bio_type_test$p.value)
  
  # wilcoxon test: disease 
  disease_test <- wilcox.test(
    batch_pcs[, i] ~ batch_pcs$DISEASE,
    mu = 0,
    alt = "two.sided",
    conf.int = T,
    conf.level = 0.99,
    paired = F,
    exact = T
  )
  p_values_matrix  [3, i] <- as.matrix(disease_test$p.value)
  
  #wilcoxon test: donor sex
  donor_sex_test <- wilcox.test(
    batch_pcs [, i] ~ batch_pcs$DONOR_SEX,
    mu = 0,
    alt = "two.sided",
    conf.int = T,
    conf.level = 0.99,
    paired = F,
    exact = T
  )
  p_values_matrix  [4, i] <- as.matrix(donor_sex_test$p.value)

  #permutation test: seq run count  
  seq_runs_count_test <-
    cor.perm (x = batch_pcs [, i], y = batch_pcs$SEQ_RUNS_COUNT)
  p_values_matrix  [5, i] <-
    as.matrix(seq_runs_count_test$p_value)
  
  # permutation test: donor age
  ## to do: find a function to apply on multiple columns at once
  donor_age_test <-
    cor.perm (x = batch_pcs [, i],
              y = c(62, 47, 72, 52, 62, 82, 67, 82, 77, 62))
  p_values_matrix  [6, i] <- as.matrix(donor_age_test$p_value)
  
  # kruskal test: submission data
  submission_date_test <-
    kruskal.test(batch_pcs [, i] ~ batch_pcs$FIRST_SUBMISSION_DATE,
                 data = batch_pcs)
  p_values_matrix  [7, i] <-
    as.matrix(submission_date_test$p.value)
  
  #wilcoxon test: cell type
  cell_type_test <- wilcox.test(
    batch_pcs [, i] ~ batch_pcs$cellTypeShort,
    mu = 0,
    alt = "two.sided",
    conf.int = T,
    conf.level = 0.99,
    paired = F,
    exact = T
  )
  p_values_matrix  [8, i] <- as.matrix(cell_type_test$p.value)
}


###creating a heatmap for p_values_matrix to see if the batch/biological effects are significant in the PCs
###and to determine which PC to use for clustering


##ggplot(p_values_matrix, aes(rownames(p_values_matrix), colnames(p_values_matrix), z= p-values)) + geom_tile(aes(fill = Z)) +
##theme_bw() +
##scale_fill_gradient(low="white", high="blue")
install.packages("gplots")
library(gplots)

my_palette <- colorRampPalette(c("indianred2", "seagreen2")) (n = 3)
color_breaks <- c(seq(0, 0.01, length = 2),
                  seq(0.011, 1, length = 2))
heatmap.2(
  p_values_matrix,
  main = "Batch and biological effects",
  trace = "none",
  margins = c(10, 12),
  cexRow = 0.8,
  Rowv = FALSE,
  Colv = FALSE,
  col = my_palette,
  breaks = color_breaks,
  sepwidth = c(0.01, 0.01),
  sepcolor = "black",
  colsep = 1:ncol(p_values_matrix),
  rowsep = 1:nrow(p_values_matrix)
)

#------------ Conclusion from heatmap: biomaterial provider, cell type and disease have significant p-values in PC1
# visualisation of the batch effects/biological factors on the PC1 against PC2 plot

#---------biomaterial provider-------

p_cluster_bio_prov <-
  ggplot(
    centers,
    aes(
      x = PC1,
      y = PC2,
      group = Samples,
      fill = input_data_csv$BIOMATERIAL_PROVIDER
    )
  ) +
  geom_point (aes(shape = as.factor(cluster), color = Samples), size = 4) +
  theme_bw() +
  ggtitle("Biomaterial provider") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )
p_cluster_bio_prov <- p_cluster_bio_prov + scale_fill_manual(values = c("indianred2", "seagreen2"))
p_cluster_bio_prov <- p_cluster_bio_prov + scale_colour_manual(values = c("indianred2", "seagreen2"))
p_cluster_bio_prov <- p_cluster_bio_prov + scale_shape_manual(values = c(16, 17))
p_cluster_bio_prov <- ggplotly(p_cluster_bio_prov)
p_cluster_bio_prov


#---------cell type------------

p_cluster_cell_type <-
  ggplot(centers,
         aes(
           x = PC1,
           y = PC2,
           group = Samples,
           fill = input_data_csv$cellTypeShort
         )) +
  geom_point (aes(shape = as.factor(cluster), color = Samples), size = 4) +
  theme_bw() +
  ggtitle("Cell type") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )
p_cluster_cell_type <- p_cluster_cell_type + scale_fill_manual(values = c("seagreen2", "indianred2"))
p_cluster_cell_type <- p_cluster_cell_type + scale_colour_manual(values = c("indianred2", "seagreen2"))
p_cluster_cell_type <- ggplotly(p_cluster_cell_type)
p_cluster_cell_type

#---------disease---------------

p_cluster_disease <-
  ggplot(centers,
         aes(
           x = PC1,
           y = PC2,
           group = Samples,
           fill = input_data_csv$DISEASE
         )) +
  geom_point (aes(shape = as.factor(cluster), color = Samples), size = 4) +
  theme_bw() +
  ggtitle("Disease") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )
p_cluster_disease <- p_cluster_disease + scale_fill_manual (values = c("indianred2", "seagreen2"))
p_cluster_disease <- p_cluster_disease + scale_colour_manual (values= c("indianred2", "seagreen2"))
p_cluster_disease <- ggplotly(p_cluster_disease)
p_cluster_disease

