#Working with PC2

##how many genes do we use for clustering

#finding the most important genes (with the biggest influence). Therefore we will look at the loading scores 
#(saved in "rotation") of the genes on PC1. Because it's not important
#whether it is positive or negative we will look at the absolute values and rank these

loading_scores <- pca_m_values$rotation[, 2]
ranked_gene_loading <- sort(abs(loading_scores), decreasing = TRUE)

##loading plots with elbow method

plot(
  ranked_gene_loading,
  main = "Loading scores of genes",
  xlab = "Genes",
  ylab = "loading scores",
  type = "b"
)

#The kink is somewhere between 0 and 20000 genes, so lets zoom in

plot(
  ranked_gene_loading[0 : 20000],
  main = "Loading scores of genes",
  xlab = "Genes",
  ylab = "loading scores",
  type = "b"
)
abline(v = 14000,
col = "red",
lty = 5,
lwd = 2)

#We will work with the top 14000 genes 

#pick out the names of the top 14000 genes
top_14000_genes <- data.frame(ranked_gene_loading[1 : 14000])

#get m values of top 3000 genes of every sample and put them into a new data frame
clustering_data <- rbind(m_values[c(as.list.data.frame(rownames(top_14000_genes))),])

#How many clusters do we need?
wss2 <-  sapply(1:5, function(k) {
  kmeans(x = clustering_data,
         centers = k,
         iter.max = 100)$tot.withinss
})
plot(
  1:5,
  wss2,
  type = "b",
  pch = 19,
  xlab = "Number of clusters K",
  ylab = "Total within-clusters sum of squares"
)

##--> it seems we need 2 clusters (kink in the curve)

#find out if healthy and cancer samples are seperated
#variable with two center positions of value of rotated data
k <-
  kmeans(
    x = t(clustering_data),
    centers = 2,
    iter.max = 100
  )#$centers
View(k)

centers2 <- kmeans(
  x = clustering_data,
  centers = 2,
  iter.max = 100
)$centers

##adding an extra column with the category of sample with which we can color the pc dots in a ggplot according to their sample group
centers2 <-
  as.data.frame(t(data.frame(rbind(
    centers2,
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
  ))))



p_cluster2 <- ggplot(centers2, aes(centers2$`1`, centers2$`2`, group = Samples))
p_cluster2 + geom_point (aes(color = Samples), size = 4) +
  theme_bw() 
p_cluster2 <- p_cluster2 + scale_colour_manual(values = c("seagreen2", "indianred2"))  


#applying t test for each gene between the contron and cancer groups, generating a p value matrix for each gene

#transposing the matrix for the t.test()
transposed_clustering_data <- as.data.frame(t(clustering_data))

p_value_each_gene <-
  sapply(1:ncol(transposed_clustering_data), function(k) {
   t.test(transposed_clustering_data[1:5, k],
           transposed_clustering_data[6:10, k],
           var.equal = T)$p.value
  })

p_value_each_gene <- as.data.frame(p_value_each_gene)

# --------to do: adjust p values for multiple comparisons with p.adjust() Bonferroni-Holm ("BH") method

p_value_each_gene$BH <-  p.adjust(p_value_each_gene$p_value_each_gene, 
                                          method = "BH")
#adding the rownames (gene names) to the matrix

p_value_each_gene$Names <- rownames(clustering_data)

# setting threshold for p-values to 0.05, and keeping the genes which fulfill this condition

p_value_each_gene <- p_value_each_gene[which(p_value_each_gene$BH < 0.05), ]

rownames(p_value_each_gene) <- p_value_each_gene$Names

#leaving only the genes which fulfill the threshhold condition in the clustering_data dataset, which their corresponding m-values

clustering_data <- clustering_data[c(rownames(p_value_each_gene)), ]

