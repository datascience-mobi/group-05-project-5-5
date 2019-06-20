#Working with PC2

##how many genes do we use for clustering

#finding the most important genes (with the biggest influence). Therefore we will look at the loading scores 
#(saved in "rotation") of the genes on PC1. Because it's not important
#whether it is positive or negative we will look at the absolute values and rank these

loading_scores <- pca_m_values$rotation[, 1]
ranked_gene_loading <- sort(abs(loading_scores), decreasing = TRUE)

##loading plots with elbow method

plot(
  ranked_gene_loading,
  main = "Loading scores of genes",
  xlab = "Genes",
  ylab = "loading scores",
  type = "b"
)

#The kink is somewhere between 0 and 10000 genes, so lets zoom in

plot(
  ranked_gene_loading[0 : 10000],
  main = "Loading scores of genes",
  xlab = "Genes",
  ylab = "loading scores",
  type = "b"
)
abline(v = 2000,
col = "red",
lty = 5,
lwd = 2)

#We will work with the top 2000 genes 

#pick out the names of the top 2000 genes
top_2000_genes <- data.frame(ranked_gene_loading[1 : 2000])

#get m values of top 3000 genes of every sample and put them into a new data frame
clustering_data <- rbind(m_values[c(as.list.data.frame(rownames(top_2000_genes))),])

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
centers2 <-
  kmeans(
    x = clustering_data,
    centers = 2,
    iter.max = 100
  )$centers

##adding an extra column with the category of sample with which we can color the pc dots in a ggplot according to their sample group
centers2 <-
  data.frame(rbind(
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
  ))

centers2 <- as.data.frame(t(centers2))

p_cluster2 <- ggplot(centers2, aes(centers2$`1`, centers2$`2`, group = Samples))
p_cluster2 + geom_point (aes(color = Samples), size = 4) +
  theme_bw() 
p_cluster2 <- p_cluster2 + scale_colour_manual(values = c("seagreen2", "indianred2"))  
