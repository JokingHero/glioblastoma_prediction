rm(list=ls(all = TRUE))
gc(reset=TRUE)

library(caret)
library(doParallel)
registerDoParallel(3)

kfold <- 10
trCtrl <- trainControl(method = "repeatedcv",
                       number = 2,
                       repeats = 2,
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
               trControl = trCtrl,
               preProc = c("center", "scale"))
  return(ret)
}

load("Glioblastoma.RData")

levels(Y.train) <- c("D","A")
combined_datasets <- cbind(X.train, Z.train[,-which("performance_score" == colnames(Z.train))])
combined_datasets$survival <- Y.train

combined_datasets <- combined_datasets[complete.cases(combined_datasets),]
survival_col <- which("survival" == colnames(combined_datasets))

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
    RF_model <- trainFun(train_set[, -survival_col], train_set[, survival_col], trCtrl)
    
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
