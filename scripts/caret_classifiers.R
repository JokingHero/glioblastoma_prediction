rm(list=ls(all = TRUE))
gc(reset=TRUE)

library(caret)
library(doParallel)
registerDoParallel(3)

kfold <- 10
trCtrl <- trainControl(method = "repeatedcv",
                       number = 4,
                       repeats = 4,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       allowParallel = TRUE)

trainFun <- function(x, y, trCtrl) {
  ret <- train(x,
               y,
               method = "parRF",
               metric = "ROC",
               family = binomial,
               do.trace = 10,
               trControl = trCtrl)
  return(ret)
}

load("Glioblastoma.RData")

levels(Y.train) <- c("D","A")
combined_datasets <- cbind(X.train, Z.train)
combined_datasets$survival <- Y.train

combined_datasets <- combined_datasets[complete.cases(combined_datasets),]

# k-fold crossvalidation to use for later
# Create 10 equal size folds
# folds <- cut(seq(1, length(uguides)), breaks=10, labels=FALSE)
# 
#   for (k in seq_len(kfold)) {
#     # Segement your data by fold
#     testIndexes <- which(folds == k, arr.ind = TRUE)
#     test_guides <- uguides[testIndexes]
#     train_guides <- uguides[-testIndexes]
#     
#     train_data <- g_names %in% train_guides
#     train_data <- this_sample[train_data, ]

survival_col <- which("survival" == colnames(combined_datasets))
# Train
RF_model <- trainFun(combined_datasets[, -survival_col], combined_datasets[, survival_col], trCtrl)
# Predict
predictions <- predict(RF_model, newdata =  your_new_data, type = "prob")