Are differentially methylated regions within genes associated with mantle cell lymphoma?
========================================================
author: Pascal Lafrenz, Mari Hambardzumyan, Lea Herzel, Franziska Lam
date: July 24. 2019
autosize: true

Project milestones
========================================================

$$~$$




1. Data processing 
2. Data normalization and visualization 
3. Data reduction
4. Regression and interpretation

1. Data processing
========================================================
Initial goal: retail 90% of the information after data processing.
$$~$$

How many genes across 10 samples in total?



```r
dim(Gene_data_frame)
```

```
[1] 56175    30
```
$$~$$

Problem: Methylation differences in sex chromosomes

Solution:

```r
Gene_data_frame_x_y <- 
  Gene_data_frame[-which(Gene_data_frame$Chromosome == "chrX"),]
Gene_data_frame_x_y <- 
  Gene_data_frame_x_y[-which(Gene_data_frame_x_y$Chromosome == "chrY"),]
```

How do the coverages look?
========================================================



![plot of chunk unnamed-chunk-7](final Presentation-figure/unnamed-chunk-7-1.png)![plot of chunk unnamed-chunk-7](final Presentation-figure/unnamed-chunk-7-2.png)


Problem: Unreliable coverages
========================================================
$$~$$

Solution: Setting a threshold

```r
threshold_cancer_lower <-
  quantile(mean_cancer_coverage,
           probs = 0.05,
           na.rm = TRUE)
threshold_cancer_upper <-
  quantile(mean_cancer_coverage,
           probs = 0.999,
           na.rm = TRUE)

threshold_healthy_lower <-
  quantile(mean_healthy_coverage,
           probs = 0.05,
           na.rm = TRUE)
threshold_healthy_upper <-
  quantile(mean_healthy_coverage,
           probs = 0.999,
           na.rm = TRUE)
```

Applying thresholds to the coverages
========================================================
$$~$$


```r
cancer_threshold_function <- function(cancer_coverage) {
  if(cancer_coverage <= threshold_cancer_lower) {
    return("NA")}
  else {return(cancer_coverage)}
  
  if(cancer_coverage >= threshold_cancer_upper) {
    return("NA")}
  else{return(cancer_coverage)}
}

cancer_coverage <- apply(cancer_coverage, MARGIN = c(1,2), FUN = cancer_threshold_function)

cancer_coverage[cancer_coverage == "NA"] <- NA 
cancer_beta_values[cancer_coverage == "NA"] <- NA
```

Problem: NA's in beta-values
========================================================
$$~$$


Solution: 

```r
cancer_beta_values <-
  cancer_beta_values[-which(
    cancer_beta_values$Number_of_NA_cancer >= 3 |
      cancer_beta_values$Number_of_NA_healthy >= 3
  ),]
healthy_beta_values <-
  healthy_beta_values[-which(
    healthy_beta_values$Number_of_NA_cancer >= 3 |
      healthy_beta_values$Number_of_NA_healthy >= 3
  ),]
sum(rownames(healthy_beta_values) != rownames(cancer_beta_values))
```

```
[1] 0
```

Faith of the remaining NA's
========================================================


```r
transposed_cancer_beta_values <- t(cancer_beta_values)
transposed_healthy_beta_values <- t(healthy_beta_values)

for (i in 1:ncol(transposed_cancer_beta_values)) {
  transposed_cancer_beta_values[is.na(transposed_cancer_beta_values[, i]), i] <-
    mean(transposed_cancer_beta_values[, i], na.rm = TRUE)
}
for (i in 1:ncol(transposed_healthy_beta_values)) {
  transposed_healthy_beta_values[is.na(transposed_healthy_beta_values[, i]), i] <-
    mean(transposed_healthy_beta_values[, i], na.rm = TRUE)
}

cancer_beta_values <- data.frame(t(transposed_cancer_beta_values))
healthy_beta_values <- data.frame(t(transposed_healthy_beta_values))
```
Checking if initial goal fulfilled:

```r
dim(cancer_beta_values)/dim(Gene_data_frame_x_y)
```

```
[1] 0.9837105 0.2333333
```

2. Data normalization and visualization
========================================================

Turning beta-values into m-values

```r
cancer_m_values <-
  data.frame(log2(cancer_beta_values / (1 - cancer_beta_values)))
healthy_m_values <-
  data.frame(log2(healthy_beta_values / (1 - healthy_beta_values)))
```



Check if transformation successful
========================================================
![plot of chunk unnamed-chunk-16](final Presentation-figure/unnamed-chunk-16-1.png)



Comparing mean m-values
========================================================
$$~$$



![plot of chunk unnamed-chunk-19](final Presentation-figure/unnamed-chunk-19-1.png)




Reducing Data
========================================================

```r
cancer_beta_values[cancer_beta_values == 0] <- 0.00001
cancer_beta_values[cancer_beta_values == 1] <- 0.99999
healthy_beta_values[healthy_beta_values == 0] <- 0.00001
healthy_beta_values[healthy_beta_values == 1] <- 0.99999

cancer_m_values <-
  data.frame(log2(cancer_beta_values / (1 - cancer_beta_values)))
healthy_m_values <-
  data.frame(log2(healthy_beta_values / (1 - healthy_beta_values)))
```

unn?tig? Da vorher schon mal berechnet

======================================================

```r
colnames(healthy_m_values) <-
  c(
    "Bcell_naive_VB_NBC_NC11_41.M",
    "Bcell_naive_VB_NBC_NC11_83.M",
    "Bcell_naive_VB_S001JP51.M",
    "Bcell_naive_VB_S00DM851.M",
    "Bcell_naive_VB_S01ECGA1.M"
  )

colnames(cancer_m_values) <- c(
  "cancer_VB_S01FE8A1.M",
  "cancer_VB_S01FF6A1.M",
  "cancer_VB_S01FH2A1.M",
  "cancer_VB_S01FJZA1.M",
  "cancer_VB_S01FKXA1.M"
)
```


======================================================











































```
Error in svd(x, nu = 0, nv = k) : infinite or missing values in 'x'
```
