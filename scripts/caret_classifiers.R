rm(list=ls(all = TRUE))
gc(reset=TRUE)

# Data
load("Glioblastoma.RData")
levels(Y.train) <- c("D","A")

# Filter lowly expressed and low variance genes
library(genefilter)
exprs_set <- ExpressionSet(assayData = t(X.train))
exprs_set$survival <- Y.train

n1 <- sum(Y.train == "A")
n2 <- sum(Y.train == "D")
S <- rowSds(exprs(exprs_set))
temp <- rowttests(exprs_set, exprs_set$survival)
d <- temp$dm
p <- temp$p.value
t <- temp$statistic
S_cutoff <- quantile(S, .5)
filter_volcano(d, p, S, n1, n2, alpha = .01, S_cutoff)

# filter
f2 <- ttest(exprs_set$survival, p=0.01) #73 genes
ffun <- filterfun(f2)
selected <- genefilter(exprs(exprs_set), ffun)
selected <- names(selected[selected])

X.train <- X.train[, colnames(X.train) %in% selected]
X.test <- X.test[, colnames(X.test) %in% selected]

combined_datasets <- data.frame(X.train, Y.train)
ggplot(combined_datasets, aes(Y.train, mRNA_TOX, fill = Y.train)) + geom_boxplot()

# Scale microarray values #z-score vs min-max vs Generalized Logistic scaling
# GL scaling using combined test+train sets so that it can generalize better
source("./scripts/General_Logistic_scaling.R")
data_to_scale <- rbind(X.test, X.train)
data_scaled <- glscale(data_to_scale)[[1]]
colnames(data_scaled) <- colnames(data_to_scale)
rownames(data_scaled) <- rownames(data_to_scale)
X.test <- data_scaled[rownames(data_scaled) %in% rownames(X.test),
                      colnames(data_scaled) %in% colnames(X.test)]
X.train <- data_scaled[rownames(data_scaled) %in% rownames(X.train),
                       colnames(data_scaled) %in% colnames(X.train)]

# performance score
library(missForest)
library(Hmisc)
combined_datasets <- data.frame(rbind(X.test, X.train), 
                                rbind(Z.test, Z.train))

table(combined_datasets$performance_score)
iris.imp <- missForest(combined_datasets)
aregImpute(~ performance_score, data = combined_datasets, n.impute = 5, match='closest')

trCtrl <- trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 10,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary)#,
# add strata 70/30
#allowParallel = TRUE)

trainFun <- function(x, y, trCtrl) {
  ret <- train(x,
               y,
               method = "ranger",
               metric = "ROC",
               family = binomial,
               #do.trace = 10,
               trControl = trCtrl,
               preProc = c("center", "scale"),
               strata = y,
               sampsize = rep(min(table(y)), 2)) # sample down
  return(ret)
}


#dummify factors
Z.train$performance_score <- factor(Z.train$performance_score)
dmy <- dummyVars(" ~ .", data = Z.train)
Z.train <- data.frame(predict(dmy, newdata = Z.train))


# age is important ggplot(abc, aes(age, fill = Y.train)) + geom_histogram(aes(y =..density..)) + geom_density(alpha=0.3)
# performance is important when you are old
# males have much higher chance of dying, female random

combined_datasets <- data.frame(X.train, Z.train)
combined_datasets$survival <- Y.train

# use Amelia to compute missing values
# 1. Expression set follows Multivariate Normal Distribution
# 2. Missing data is random in nature (Assumption)
library(Amelia)
amelia_fit <- amelia(combined_datasets, m=5, parallel = "multicore", noms = c("survival"))
iris.imp <- missForest(combined_datasets)

combined_datasets <- combined_datasets[complete.cases(combined_datasets),]
survival_col <- which("survival" == colnames(combined_datasets))


# Classifiers 

library(caret)
library(doParallel)
registerDoParallel(3)

kfold <- 2
trCtrl <- trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 10,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary)#,
# add strata 70/30
#allowParallel = TRUE)

trainFun <- function(x, y, trCtrl) {
  ret <- train(x,
               y,
               method = "ranger",
               metric = "ROC",
               family = binomial,
               #do.trace = 10,
               trControl = trCtrl,
               preProc = c("center", "scale"))
  return(ret)
}

# repeated n times k-fold crossvalidation
# Create 10 equal size folds
CV <- c()
for (n in 1:5) {
  
  folds <- cut(seq(1, dim(combined_datasets)[1]), breaks=kfold, labels=FALSE)
  for (k in seq_len(kfold)) {
    # Segement your data by fold
    testIndexes <- which(folds == k, arr.ind = TRUE)
    test_set <- combined_datasets[testIndexes,]
    train_set <- combined_datasets[-testIndexes,]
    
    # Train
    model <- trainFun(train_set[, -survival_col], train_set[, survival_col], trCtrl)
    
    test_classes <- predict(RF_model, newdata = test_set[, -survival_col])
    CV <- c(CV, confusionMatrix(data = test_classes, test_set$survival))
    
    print(n)
    print(k)
  }
}


# Final Train Predict
# Train
RF_model <- trainFun(combined_datasets[, -survival_col], combined_datasets[, survival_col], trCtrl)
# Predict
combined_test <- cbind(X.test, Z.test[,-which("performance_score" == colnames(Z.test))])
predictions <- predict(RF_model, newdata =  combined_test, type = "prob")
