
##define a function to set every value of cancer coverage and cancer beta value to NA if they are in threshold and apply for the entire dataframe

function_cancer <- function(x = cancer_coverage) {
  if (cancer_coverage <= threshold) {
    cancer_coverage[cancer_coverage <= threshold] <- 0
  }
  
  if (cancer_coverage >= threshold2) {
    cancer_coverage[cancer_coverage >= threshold2] <- 0
  }
  
  if (cancer_coverage == 0) {
    cancer_coverage[cancer_coverage == 0] <- NA
    cancer_beta_values[i, j] <- NA
  }}
apply(cancer_coverage, MARGIN = c(1,2), FUN = function_cancer)


##define a function to set every value of healthy coverage and cancer beta value to NA if they are in threshold and apply for the entire dataframe

function_healthy <- function(x = healthy_coverage) {
  if (healthy_coverage <= threshold3) {
    healthy_coverage[healthy_coverage <= threshold3] <- 0
  }
  
  if (healthy_coverage >= threshold4) {
    healthy_coverage[healthy_coverage >= threshold4] <- 0
  }
  
  if (healthy_coverage == 0) {
    healthy_coverage[healthy_coverage == 0] <- NA
    healthy_coverage[i, j] <- NA
  }}
apply(healthy_coverage, MARGIN = c(1,2), FUN = function_healthy)
