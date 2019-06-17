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

#Apply PCA on data frame "complete_m_values" with all m values. For that, the matrix needs to be transposed first
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
#generate a ggplot/scatterplot to visualize the Sample points in a coordinate system with x-axis = PC1 and y-axis = PC2
p <- ggplot(pcs_of_m_values, aes(PC1, PC2, group = Samples)) +
  geom_point (aes(shape = Samples, color = Samples), size = 4)
p + scale_colour_manual(values = c("seagreen2", "indianred2"))

#finding the top 25 most important genes (with the biggest influence). Therefore we will look at the loading scores (saved in "rotation") of the genes on PC1. Because it's not important
#whether it is positive or negative we will look at the absolute values and rank these

loading_scores <- pca_m_values$rotation[, 1]
ranked_gene_loading <- sort(abs(loading_scores), decreasing = TRUE)
top_25_genes <- names(ranked_gene_loading[1:25])
View(top_25_genes)

##loading plots with elbow method

pca_m_values$rotation[top_25_genes, 1]

#find out how much clusters do we need to group samples (obvisiously 2 would be perfekt because healthy/cancer)

##why sapply??

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
p_cluster <- ggplot(centers, aes(X1, X2, group = Samples)) +
  geom_point (aes(shape = Samples, color = Samples), size = 4)
p_cluster + scale_colour_manual(values = c("seagreen2", "indianred2"))

## wilkoxon, kruskal wallis, Pearson correlation coefficient berechnen und ein permutation test


pcs_of_m_values$PC1 <- as.numeric(as.character(pcs_of_m_values$PC1))
pcs_of_m_values$PC2 <- as.numeric(as.character(pcs_of_m_values$PC2))
pcs_of_m_values$PC3 <- as.numeric(as.character(pcs_of_m_values$PC3))

batch_pcs <-
  cbind(pcs_of_m_values[, 1:3], c(input_data_csv[, c(
    "BIOMATERIAL_PROVIDER",
    "BIOMATERIAL_TYPE",
    "SAMPLE_DESC_3",
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

#wilcoxon test
#--------PC1---------

bio_prov_test_pc1 <- wilcox.test(
  batch_pcs$PC1 ~ batch_pcs$BIOMATERIAL_PROVIDER,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)

bio_type_test_pc1 <- wilcox.test(
  batch_pcs$PC1 ~ batch_pcs$BIOMATERIAL_TYPE,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)


disease_test_pc1 <- wilcox.test(
  batch_pcs$PC1 ~ batch_pcs$DISEASE,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)

donor_sex_test_pc1 <- wilcox.test(
  batch_pcs$PC1 ~ batch_pcs$DONOR_SEX,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)

#-------PC2----------

bio_prov_test_pc2 <- wilcox.test(
  batch_pcs$PC2 ~ batch_pcs$BIOMATERIAL_PROVIDER,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)

bio_type_test_pc2 <- wilcox.test(
  batch_pcs$PC2 ~ batch_pcs$BIOMATERIAL_TYPE,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)


disease_test_pc2 <- wilcox.test(
  batch_pcs$PC2 ~ batch_pcs$DISEASE,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)

donor_sex_test_pc2 <- wilcox.test(
  batch_pcs$PC2 ~ batch_pcs$DONOR_SEX,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)

#---------PC3---------

bio_prov_test_pc3 <- wilcox.test(
  batch_pcs$PC3 ~ batch_pcs$BIOMATERIAL_PROVIDER,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)

bio_type_test_pc3 <- wilcox.test(
  batch_pcs$PC3 ~ batch_pcs$BIOMATERIAL_TYPE,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)


disease_test_pc3 <- wilcox.test(
  batch_pcs$PC3 ~ batch_pcs$DISEASE,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)

donor_sex_test_pc3 <- wilcox.test(
  batch_pcs$PC3 ~ batch_pcs$DONOR_SEX,
  mu = 0,
  alt = "two.sided",
  conf.int = T,
  conf.level = 0.95,
  paired = F,
  exact = T
)

# to do: find a function to apply on multiple columns at once


#permutation test

# on seq runs count

#--------PC1-------

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
seq_runs_count_test_pc1 <-
  cor.perm (x = batch_pcs$PC1, y = batch_pcs$SEQ_RUNS_COUNT)


#permutation test on donor age
#working with mean donor ages


donor_age_test_pc1 <-
  cor.perm (x = batch_pcs$PC1,
            y = c(62, 47, 72, 52, 62, 82, 67, 82, 77, 62))

#--------PC2-------


seq_runs_count_test_pc2 <-
  cor.perm (x = batch_pcs$PC2, y = batch_pcs$SEQ_RUNS_COUNT)


#permutation test on donor age


donor_age_test_pc2 <-
  cor.perm (x = batch_pcs$PC2,
            y = c(62, 47, 72, 52, 62, 82, 67, 82, 77, 62))

#--------PC3-------


seq_runs_count_test_pc3 <-
  cor.perm (x = batch_pcs$PC3, y = batch_pcs$SEQ_RUNS_COUNT)


#permutation test on donor age


donor_age_test_pc3 <-
  cor.perm (x = batch_pcs$PC3,
            y = c(62, 47, 72, 52, 62, 82, 67, 82, 77, 62))


#kruskal wallis test

#------PC1-------
#submission date

submission_date_test_pc1 <- kruskal.test(batch_pcs$PC1 ~ batch_pcs$FIRST_SUBMISSION_DATE,
                                         data = batch_pcs)

#kruskal wallis on cell type

cell_type_test_pc1 <- kruskal.test(batch_pcs$PC1 ~ batch_pcs$SAMPLE_DESC_3,
                                   data = batch_pcs)

#------PC2-------

#submission date

submission_date_test_pc2 <- kruskal.test(batch_pcs$PC2 ~ batch_pcs$FIRST_SUBMISSION_DATE,
                                         data = batch_pcs)

#kruskal wallis on cell type

cell_type_test_pc2 <- kruskal.test(batch_pcs$PC2 ~ batch_pcs$SAMPLE_DESC_3,
                                   data = batch_pcs)

#------PC3-------

#submission date

submission_date_test_pc3 <- kruskal.test(batch_pcs$PC3 ~ batch_pcs$FIRST_SUBMISSION_DATE,
                                         data = batch_pcs)

#kruskal wallis on cell type

cell_type_test_pc3 <- kruskal.test(batch_pcs$PC3 ~ batch_pcs$SAMPLE_DESC_3,
                                   data = batch_pcs)


#creating a matrix with p values
p_values_matrix = matrix(
  c(
    bio_prov_test_pc1$p.value ,
    bio_type_test_pc1$p.value,
    disease_test_pc1$p.value,
    donor_sex_test_pc1$p.value,
    seq_runs_count_test_pc1$p_value,
    donor_age_test_pc1$p_value,
    submission_date_test_pc1$p.value,
    cell_type_test_pc1$p.value,
    bio_prov_test_pc2$p.value ,
    bio_type_test_pc2$p.value,
    disease_test_pc2$p.value,
    donor_sex_test_pc2$p.value,
    seq_runs_count_test_pc2$p_value,
    donor_age_test_pc2$p_value,
    submission_date_test_pc2$p.value,
    cell_type_test_pc2$p.value,
    bio_prov_test_pc3$p.value ,
    bio_type_test_pc3$p.value,
    disease_test_pc3$p.value,
    donor_sex_test_pc3$p.value,
    seq_runs_count_test_pc3$p_value,
    donor_age_test_pc3$p_value,
    submission_date_test_pc3$p.value,
    cell_type_test_pc3$p.value
  ),
  nrow = 8,
  ncol = 3
)
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
colnames(p_values_matrix) <- c("PC1", "PC2", "PC3")
