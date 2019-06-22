


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

# calculating mean, sd values and plotting a histogramm for each
mean_cancer_m_values <- rowMeans(cancer_m_values)
hist(
  log10(mean_cancer_m_values),
  breaks = "fd",
  main = "Cancer M values: Mean frequency",
  xlab = "Common logarithm of M values",
  col = "indianred2",
  border = "gray20",
  xlim = c(-1, 1)
)
abline(v = log10(quantile(
  mean_cancer_m_values,
  probs = seq(0, 1, 0.1),
  na.rm = TRUE
)),
col = "black",
lty = 5,
lwd = 1)

mean_healthy_m_values <- rowMeans(healthy_m_values)
hist(
  log10(mean_healthy_m_values),
  breaks = "fd",
  main = "Healthy M values: Mean frequency",
  xlab = "Common logarithm of M values",
  col = "seagreen2",
  border = "gray20",
  xlim = c(-1, 1)
)
abline(v = log10(quantile(
  mean_healthy_m_values,
  probs = seq(0, 1, 0.1),
  na.rm = TRUE
)),
col = "black",
lty = 5,
lwd = 1)

sd_cancer_m_values <- apply(cancer_m_values, 1, sd)
hist(
  log10(sd_cancer_m_values),
  breaks = "fd",
  main = "Cancer M values: SD frequency",
  xlab = "Common logarithm of M values",
  col = "indianred2",
  border = "gray20"
)

sd_healthy_m_values <- apply(healthy_m_values, 1, sd)
hist(
  log10(sd_healthy_m_values),
  breaks = "fd",
  col = "seagreen2",
  main = "Healthy M values: SD frequency",
  xlab = "Common logarithm of M values",
  border = "gray20"
)


#include mean and sd values to the original dataframes

extended_cancer_m_values <-
  cbind.data.frame(cancer_m_values, mean_cancer_m_values, sd_cancer_m_values)
extended_healthy_m_values <-
  cbind.data.frame(healthy_m_values, mean_healthy_m_values, sd_healthy_m_values)

##showing mean cancer m-values vs. mean healthy m-values
#at first we have to load packages for visualization
install.packages("tidyverse")
library(tidyverse)
install.packages("ggrepel")
library(ggrepel)

#extracting values of important genes (defined earlier) therefore we can highlight important genes in diagram
#cancer
extended_cancer_m_values_gene <- extended_cancer_m_values[c(
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
),]

#healthy
extended_healthy_m_values_gene <- extended_healthy_m_values[c(
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
),]

#adding names of important genes to our important gene data frame because it is easier to work with
extended_cancer_m_values_gene <-
  cbind(
    extended_cancer_m_values_gene,
    Important_genes = c(
      "SOX11",
      "NR2F2",
      "p53",
      "CCND1",
      "AHR",
      "ROBO1",
      "SOX9",
      "HOXA9",
      "CDH1",
      "CDC14B",
      "FOXC1",
      "G0S2",
      "GPX3",
      "LGALS3",
      "MAL",
      "NPTX2",
      "PAX6",
      "TFPI2",
      "THEM4",
      "TWIST1"
    ) 
  )

#now we plot our data with following features: 
#mean values of cancer and healthy (black),
#highlighting important genes for MCL (red), 
#showing  general trend of points(blue smoothed curve; because of overplotting), 
#straight line through origin (y=x) for better interpreting

#first general overview
ggplot() +
  geom_point(
    mapping = aes(
      x = extended_cancer_m_values$mean_cancer_m_values,
      y = extended_healthy_m_values$mean_healthy_m_values
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  geom_smooth(
    mapping = aes(
      x = extended_cancer_m_values$mean_cancer_m_values,
      y = extended_healthy_m_values$mean_healthy_m_values
    ),
    na.rm = TRUE,
    alpha = 1 / 10,
    colour = "darkgrey"
  ) +
  labs(x = "Mean cancer m-values",
       y = "Mean healthy m-values",
       title = "Comparison of mean m-values") +
  theme_bw() +
  xlim(-12, 9) +
  ylim(-12, 9) 

#now detail visualization
ggplot() +
  geom_point(
    mapping = aes(
      x = extended_cancer_m_values$mean_cancer_m_values,
      y = extended_healthy_m_values$mean_healthy_m_values
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  geom_point(
    mapping = aes(
      x = extended_cancer_m_values_gene$mean_cancer_m_values,
      y = extended_healthy_m_values_gene$mean_healthy_m_values
    ),
    colour = "red",
    size = 2
  ) +
  geom_label_repel(
    aes(
      label = Important_genes,
      x = extended_cancer_m_values_gene$mean_cancer_m_values,
      y = extended_healthy_m_values_gene$mean_healthy_m_values
    ),
    data = extended_cancer_m_values_gene,
    point.padding = 0.5,
    label.size = 0.1,
    segment.colour = "cornflowerblue",
    segment.alpha = 0.9
  ) +
  labs(x = "Mean cancer m-values",
       y = "Mean healthy m-values",
       title = "Comparison of mean m-values") +
  theme_bw() +
  xlim(-5, 5) +
  ylim(-5, 5) +
  geom_abline(
    mapping = NULL,
    data = NULL,
    slope = 1,
    intercept = 0,
    colour = "yellow2"
  )

#showing SD of cancer m-values vs. SD of healthy m-values
#first again overall visualization
ggplot() +
  geom_point(
    mapping = aes(
      x = extended_cancer_m_values$sd_cancer_m_values,
      y = extended_healthy_m_values$sd_healthy_m_values
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  labs(x = "SD cancer m-values",
       y = "SD healthy m-values",
       title = "Comparison of SD of m-values") +
  theme_bw()

#now detail visualization
ggplot() +
  geom_point(
    mapping = aes(
      x = extended_cancer_m_values$sd_cancer_m_values,
      y = extended_healthy_m_values$sd_healthy_m_values
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  geom_point(
    mapping = aes(
      x = extended_cancer_m_values_gene$sd_cancer_m_values,
      y = extended_healthy_m_values_gene$sd_healthy_m_values
    ),
    colour = "red",
    size = 2
  ) +
  geom_label_repel(
    aes(
      label = Important_genes,
      x = extended_cancer_m_values_gene$sd_cancer_m_values,
      y = extended_healthy_m_values_gene$sd_healthy_m_values
    ),
    data = extended_cancer_m_values_gene,
    point.padding = 0.5,
    label.size = 0.1,
    segment.colour = "cornflowerblue",
    segment.alpha = 0.9,
    max.iter = 5000
  ) +
  labs(x = "SD cancer m-values",
       y = "SD healthy m-values",
       title = "Comparison of SD of m-values") +
  xlim(-1, 2.5) +
  ylim(0, 1) +
  theme_bw()

#comparing mean beta values and mean m values (for Rmarkdown: compare to literature and say if it was successfull)
#cancer
ggplot() +
  geom_point(
    mapping = aes(
      x = extended_cancer_m_values$mean_cancer_m_values,
      y = rowMeans(cancer_beta_values)
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  geom_smooth(
    mapping = aes(
      x = extended_cancer_m_values$mean_cancer_m_values,
      y = rowMeans(cancer_beta_values)
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  labs(x = "Mean cancer beta values",
       y = "Mean cancer m values",
       title = "Comparison of mean values") +
  theme_bw()

#healthy
ggplot() +
  geom_point(
    mapping = aes(
      x = extended_healthy_m_values$mean_healthy_m_values,
      y = rowMeans(healthy_beta_values)
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  geom_smooth(
    mapping = aes(
      x = extended_healthy_m_values$mean_healthy_m_values,
      y = rowMeans(healthy_beta_values)
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  labs(x = "Mean healthy beta values",
       y = "Mean healthy m values",
       title = "Comparison of mean values") +
  theme_bw()
