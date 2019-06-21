install.packages("effects")
library(effects)
install.packages("jtools")
library(jtools)
install.packages("sandwich")
library("sandwich")
install.packages("corrplot")
library(corrplot)
install.packages("GGally")
library(GGally)

#hoe much does our variables correlate
correlation_top_genes <- cor(clustering_data)
corrplot(correlation_top_genes)
ggpairs(clustering_data)
#-> pretty heavy correlation

#first creating data frame which has dichotomous outcome variable included, 0 -> healthy, 1 -> cancer
regression_data <- cbind(Health_status = as.factor(c(
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
)), as.data.frame(t(clustering_data)))

#regression model, family: logistic regression, robust: robust standard errors (HC1: stata analysis), confint = TRUE + digits = 3: confidence intervall
regression_model <- glm(formula = Health_status ~ ., family = binomial(link = "logit"), data = regression_data)
regression_model <- glm(formula = Health_status ~ ENSG00000274727 + ENSG00000259694 + ENSG00000211972 + ENSG00000186092 + ENSG00000207492 + ENSG00000272693 + ENSG00000280406 + ENSG00000273514 + ENSG00000276673, family = binomial(link = "logit"), data = regression_data)
Lr <- step(regression_model)

summary(regression_model)
summ(regression_model, robust = "HC1")
#-> many na's and "singularities" : too much correlation but how to solve it???

##overview
plot(allEffects(regression_model))

#predict on data set
predict(regression_model, newdata = regression_data)



