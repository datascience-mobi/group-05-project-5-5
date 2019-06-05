## replacing 0 and 1 in beta values with approximate values, thus no -inf and +inf values in m values

function_inf <- function(x) {
  if (x == 1) {
    return(0.99999)
  } else {return(x)}
}
function_minus_inf <- function(x) {
  if (x == 0) {
    return(0.00001)
  } else {return(x)}
}
cancer_beta_values <- as.data.frame(apply(cancer_beta_values, MARGIN = c(1,2), FUN = function_inf))
cancer_beta_values <- as.data.frame(apply(cancer_beta_values, MARGIN = c(1,2), FUN = function_minus_inf))
healthy_beta_values <- as.data.frame(apply(healthy_beta_values, MARGIN = c(1,2), FUN = function_inf))
healthy_beta_values <- as.data.frame(apply(healthy_beta_values, MARGIN = c(1,2), FUN = function_minus_inf))
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