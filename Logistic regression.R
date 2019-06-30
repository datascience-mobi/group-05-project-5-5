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
install.packages("mctest")
library(mctest)
install.packages("corpcor")
library(corpcor)
#how much does our variables correlate
correlation_top_genes <- cor(clustering_data)
corrplot(correlation_top_genes)
ggpairs(as.data.frame(t(clustering_data)))
heatmap(as.matrix(clustering_data))
hist(cor(t(clustering_data)), main = "Genes correlation", xlab = "Correlation value")
#-> pretty heavy correlation

#first creating data frame which has dichotomous outcome variable included, 0 -> healthy, 1 -> cancer
regression_data <- cbind(Health_status = as.numeric(c(
  0,
  0,
  0,
  0,
  0,
  1,
  1,
  1,
  1,
  1
)), as.data.frame(t(clustering_data)))

regression_data <- regression_data[, sample(ncol(regression_data))]
#regression model, family: logistic regression, robust: robust standard errors (HC1: stata analysis), confint = TRUE + digits = 3: confidence intervall
regression_model <- glm(formula = Health_status ~ ., family = binomial(link = "logit"), data = regression_data, maxit = 4)
Lr <- step(regression_model)

summary(regression_model)
summ(regression_model, robust = "HC1")
#-> many na's and "singularities" : too much correlation but how to solve it???

#trying to solve it through cross validation
control <-  which(regression_data$Health_status == "Healthy")
control.train <- sample(control, floor(0.75 * length(control)))
cancer <- which(regression_data$Health_status == "Cancer")
cancer.train <- sample(cancer, floor(0.75 * length(cancer)))
train <- c(control.train, cancer.train)
train.set <- regression_data[train, ]
test.set <- regression_data[-train, ]
regression_model2 <- glm(formula = Health_status ~ ., family = binomial(link = "logit"), data = train.set, maxit = 4)
predict(regression_model2, newdata = regression_data, type = "response")

##overview
plot(allEffects(regression_model))

#predict on data set
predict(regression_model, newdata = regression_data, type = "response")

              