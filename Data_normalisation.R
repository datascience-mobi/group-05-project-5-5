


# transforming beta values to M values and creating a separate dataframe with those values

cancer_m_values <-
  data.frame(log2(cancer_beta_values / (1 - cancer_beta_values)))
healthy_m_values <-
  data.frame(log2(healthy_beta_values / (1 - healthy_beta_values)))

# changing the ending of patients names from .bed to .M for better overview

# changing healthy patients names
names(healthy_m_values)[names(healthy_m_values) == "Bcell_naive_VB_NBC_NC11_41.bed"] <-
  "Bcell_naive_VB_NBC_NC11_41.M"
names(healthy_m_values)[names(healthy_m_values) == "Bcell_naive_VB_NBC_NC11_83.bed"] <-
  "Bcell_naive_VB_NBC_NC11_83.M"
names(healthy_m_values)[names(healthy_m_values) == "Bcell_naive_VB_S001JP51.bed"] <-
  "Bcell_naive_VB_S001JP51.M"
names(healthy_m_values)[names(healthy_m_values) == "Bcell_naive_VB_S00DM851.bed"] <-
  "Bcell_naive_VB_S00DM851.M"
names(healthy_m_values)[names(healthy_m_values) == "Bcell_naive_VB_S01ECGA1.bed"] <-
  "Bcell_naive_VB_S01ECGA1.M"

# changing cancer patients names
names(cancer_m_values)[names(cancer_m_values) == "cancer_VB_S01FE8A1.bed"] <-
  "cancer_VB_S01FE8A1.M"
names(cancer_m_values)[names(cancer_m_values) == "cancer_VB_S01FF6A1.bed"] <-
  "cancer_VB_S01FF6A1.M"
names(cancer_m_values)[names(cancer_m_values) == "cancer_VB_S01FH2A1.bed"] <-
  "cancer_VB_S01FH2A1.M"
names(cancer_m_values)[names(cancer_m_values) == "cancer_VB_S01FJZA1.bed"] <-
  "cancer_VB_S01FJZA1.M"
names(cancer_m_values)[names(cancer_m_values) == "cancer_VB_S01FKXA1.bed"] <-
  "cancer_VB_S01FKXA1.M"

# calculating mean, sd values and plotting a histogramm for each
mean_cancer_m_values <- rowMeans(cancer_m_values)
hist(
  log10(mean_cancer_m_values),
  breaks = "fd",
  main = "Cancer M values: Mean frequency",
  xlab = "Common logarithm of M values",
  col = "indianred1",
  border = "gray20",
  xlim = c(-1, 1)
)
abline(v = log10(quantile(
  mean_cancer_m_values,
  probs = seq(0, 1, 0.1),
  na.rm = TRUE
)),
col = "black",
lwd = 2)

mean_healthy_m_values <- rowMeans(healthy_m_values)
hist(
  log10(mean_healthy_m_values),
  breaks = "fd",
  main = "Healthy M values: Mean frequency",
  xlab = "Common logarithm of coverages",
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
lwd = 2)

sd_cancer_m_values <- apply(cancer_m_values, 1, sd)
hist(
  log10(sd_cancer_m_values),
  breaks = "fd",
  main = "Cancer M values: SD frequency",
  xlab = "Common logarithm of M values",
  col = "indianred1",
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

complete_cancer_m_values <-
  cbind.data.frame(cancer_m_values, mean_cancer_m_values, sd_cancer_m_values)
complete_healthy_m_values <-
  cbind.data.frame(healthy_m_values, mean_healthy_m_values, sd_healthy_m_values)

#showing mean cancer m-values vs. mean healthy m-values
install.packages("tidyverse")
library(tidyverse)

ggplot() +
  geom_point(
    mapping = aes(
      x = complete_cancer_m_values$mean_cancer_m_values,
      y = complete_healthy_m_values$mean_healthy_m_values
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  geom_smooth(
    mapping = aes(
      x = complete_cancer_m_values$mean_cancer_m_values,
      y = complete_healthy_m_values$mean_healthy_m_values
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  labs(x = "Mean cancer m-values",
       y = "Mean healthy m-values",
       title = "Comparison of mean m-values") +
  theme_bw() +
  xlim(-15, 10) + 
  ylim(-15, 10)

#showing SD of cancer m-values vs. SD of healthy m-values
ggplot() +
  geom_point(
    mapping = aes(
      x = complete_cancer_m_values$sd_cancer_m_values,
      y = complete_healthy_m_values$sd_healthy_m_values
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  labs(x = "Mean cancer m-values",
       y = "Mean healthy m-values",
       title = "Comparison of SD of m-values") +
  theme_bw() 

#comparing mean beta values and mean m values 
#cancer
ggplot() +
  geom_point(
    mapping = aes(
      x = complete_cancer_m_values$mean_cancer_m_values,
      y = rowMeans(cancer_beta_values)
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  geom_smooth(
    mapping = aes(
      x = complete_cancer_m_values$mean_cancer_m_values,
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
      x = complete_healthy_m_values$mean_healthy_m_values,
      y = rowMeans(healthy_beta_values)
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  geom_smooth(
    mapping = aes(
      x = complete_healthy_m_values$mean_healthy_m_values,
      y = rowMeans(healthy_beta_values)
    ),
    na.rm = TRUE,
    alpha = 1 / 10
  ) +
  labs(x = "Mean healthy beta values",
       y = "Mean healthy m values",
       title = "Comparison of mean values") +
  theme_bw() 
 
