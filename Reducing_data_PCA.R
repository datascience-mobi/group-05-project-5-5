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
complete_m_values <- cbind(healthy_m_values, cancer_m_values)
View(complete_m_values)

#Apply PCA on data frame "complete_m_values" with all m values. For that, the matrix needs to be transposed first
#(variables are scaled to have i) standard deviation one and ii) mean zero)
complete_m_values.pca <- prcomp(t(complete_m_values))
summary(complete_m_values.pca)

#visualize pca variation

###why does the xlab not function
## cancer to MCL and healthy to control??????


plot(
  complete_m_values.pca,
  main = "Variance explained through every principal component",
  type = "l"
)

#visualize pca
#plot(complete_m_values.pca$x[, 1], complete_m_values.pca$x[, 2])


#adding an extra column with the category of sample with which we can color the pc dots in a ggplot according to their sample group
pcs_of_m_values <-
  data.frame(cbind(
    complete_m_values.pca$x,
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
#generate a ggplot/scatterplot to visualize the Sample points in a coordinate system with x-axis = PC1 and y-axis = PC2
p <- ggplot(pcs_of_m_values, aes(PC1, PC2, group = Samples)) +
  geom_point (aes(shape = Samples, color = Samples), size = 4)
p + scale_colour_manual(values = c("seagreen2", "indianred1"))

#finding the top 25 most important genes (with the biggest influence). Therefore we will look at the loading scores (saved in "rotation") of the genes on PC1. Because it's not important
#whether it is positive or negative we will look at the absolute values and rank these

loading_scores <- complete_m_values.pca$rotation[, 1]
ranked_gene_loading <- sort(abs(loading_scores), decreasing = TRUE)
top_25_genes <- names(ranked_gene_loading[1:25])
View(top_25_genes)

##loading plots with elbow method

complete_m_values.pca$rotation[top_25_genes, 1]

#find out how much clusters do we need to group samples (obvisiously 2 would be perfekt because healthy/cancer)

##why sapply??

wss <-  sapply(1:5, function(k) {
  kmeans(x = complete_m_values.pca$x,
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

##adding an extra column with the category of sample with which we can color the pc dots in a ggplot according to their sample group
centers <-
  data.frame(cbind(
    t(centers),
    Samples = c(
      "Control",
      "Control",
      "Control",
      "Control",
      "Control",
      "MCL",
      "MCL",
      "MCL",
      "MCL",
      "MCL"
    )
  ))

#visualize cluster x1 and x2 and how samples are seperated
p_cluster <- ggplot(centers, aes(X1, X2, group = Samples)) +
  geom_point (aes(shape = Samples, color = Samples), size = 4)
p_cluster + scale_colour_manual(values = c("seagreen2", "indianred1"))

## wilkoxon, kruskal wallis, Pearson correlation coefficient berechnen und ein permutation test
