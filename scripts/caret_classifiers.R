rm(list=ls(all = TRUE))
gc(reset=TRUE)

# Data
load("Glioblastoma.RData")
levels(Y.train) <- c("D","A")

# Filter lowly expressed and low variance genes
library(genefilter)
library("Biobase")
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
f2 <- ttest(exprs_set$survival, p=0.05) #73 genes
ffun <- filterfun(f2)
selected <- genefilter(exprs(exprs_set), ffun)
selected <- names(selected[selected])

X.train <- X.train[, colnames(X.train) %in% selected]
X.test <- X.test[, colnames(X.test) %in% selected]

combined_datasets <- data.frame(X.train, Y.train)
ggplot(combined_datasets, aes(Y.train, mRNA_RESP18, fill = Y.train)) + geom_boxplot()

# Scale microarray values #z-score vs min-max vs Generalized Logistic scaling
# GL scaling using combined test+train sets so that it can generalize better
source("./scripts/General_Logistic_scaling.R")
data_to_scale <- rbind(X.test, X.train)
data_scaled <- glscale(data_to_scale)[[1]]
colnames(data_scaled) <- colnames(data_to_scale)
rownames(data_scaled) <- rownames(data_to_scale)
data_scaled <- data_scaled[,complete.cases(t(data_scaled))]
X.test <- data_scaled[rownames(data_scaled) %in% rownames(X.test),
                      colnames(data_scaled) %in% colnames(X.test)]
X.train <- data_scaled[rownames(data_scaled) %in% rownames(X.train),
                       colnames(data_scaled) %in% colnames(X.train)]

#set.seed(42)
rf_sa <- safs(x = data.frame(X.train, Z.train), y = Y.train,
              iters = 500,
              safsControl = safsControl(functions = rfSA,
                                        method = "repeatedcv",
                                        repeats = 5,
                                        improve = 50))
rf_sa

# performance
Z.train$performance_score <- factor(Z.train$performance_score)
levels(Z.train$performance_score) <- c("l", "l", "n", "g", "e")
Z.test$performance_score <- factor(Z.test$performance_score)
levels(Z.test$performance_score) <- c("l", "n", "g", "e")
combined_datasets <- data.frame(rbind(X.train, X.test), 
                                rbind(Z.train, Z.test))

perf_col <- which("performance_score" == colnames(combined_datasets))
train_set <- combined_datasets[complete.cases(combined_datasets), -perf_col]
pred_me <- combined_datasets[complete.cases(combined_datasets), perf_col]
model <- train(train_set,
               pred_me,
               method = "gbm",
               tuneLength = 5,
               trControl = trainControl(method = "repeatedcv",
                                        number = 10,
                                        repeats = 4),
               preProc = c("center", "scale"))
combined_datasets$performance_score[!complete.cases(combined_datasets)] <- predict(model, newdata = combined_datasets[!complete.cases(combined_datasets), -perf_col])
Z.test <- combined_datasets[rownames(combined_datasets) %in% rownames(Z.test),
                      colnames(combined_datasets) %in% colnames(Z.test)]
Z.train <- combined_datasets[rownames(combined_datasets) %in% rownames(Z.train),
                       colnames(combined_datasets) %in% colnames(Z.train)]
#dummify factors
levels(Z.train$performance_score) <- c(40, 60, 80, 100)
Z.train$performance_score <- as.numeric(as.character(Z.train$performance_score))
dmy <- dummyVars(" ~ .", data = Z.train)
Z.train <- data.frame(predict(dmy, newdata = Z.train))

combined_datasets <- data.frame(X.train, Z.train)
combined_datasets$survival <- Y.train
#combined_datasets <-combined_datasets[,-which("performance_score" == colnames(combined_datasets))]
combined_datasets <- combined_datasets[complete.cases(combined_datasets),]
survival_col <- which("survival" == colnames(combined_datasets))

# Classifiers 
library(caret)

# todo divide dataset
x = combined_datasets[,-survival_col]
y = combined_datasets[,survival_col]
model <- train(x,
               y,
               method = "svmRadial",
               metric = "ROC",
               family = binomial,
               tuneLength = 5,
               trControl = trainControl(method = "repeatedcv",
                                        number = 10,
                                        repeats = 4,
                                        classProbs = TRUE,
                                        summaryFunction = twoClassSummary),
               preProc = c("center", "scale"),
               strata = y,
               sampsize = rep(40, 2))
model

#test_classes <- predict(model, newdata = [, -survival_col])
#CV <- c(CV, confusionMatrix(data = test_classes, test_set$survival))

# Final Train Predict
# Train
RF_model <- trainFun(combined_datasets[, -survival_col], combined_datasets[, survival_col], trCtrl)
# Predict
combined_test <- cbind(X.test, Z.test[,-which("performance_score" == colnames(Z.test))])
predictions <- predict(RF_model, newdata =  combined_test, type = "prob")