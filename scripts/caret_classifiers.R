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
f2 <- ttest(exprs_set$survival, p=0.01) #73 genes
ffun <- filterfun(f2)
selected <- genefilter(exprs(exprs_set), ffun)
selected <- names(selected[selected])

X.train <- X.train[, colnames(X.train) %in% selected]
X.test <- X.test[, colnames(X.test) %in% selected]

combined_datasets <- data.frame(X.train, Y.train)
library(ggplot2)
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
# rf_sa <- safs(x = data.frame(X.train, Z.train[,-3]), y = Y.train,
#               iters = 500,
#               safsControl = safsControl(functions = rfSA,
#                                         method = "repeatedcv",
#                                         repeats = 5,
#                                         improve = 50))
# rf_sa
# rf_sa$optVariables #not improving much
# saveRDS(rf_sa, "./rf_sa.rds")
# rf_sa <- readRDS("./rf_sa.rds")

# rf_ga <- gafs(x = data.frame(X.train, Z.train[,-3]), y = Y.train,
#               iters = 600,
#               gafsControl = gafsControl(functions = rfGA,
#                                         method = "repeatedcv",
#                                         repeats = 5))
# rf_ga

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

library(caret)
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
#combined_datasets <- combined_datasets[, colnames(combined_datasets) %in% c(rf_sa$optVariables, "survival")]
survival_col <- which("survival" == colnames(combined_datasets))

# Classifiers 
x = combined_datasets[,-survival_col]
y = combined_datasets[,survival_col]
grid <- expand.grid(C=c(0.10,0.15,0.25,0.35, 0.45, 0.5, 0.55, 0.6), 
                    sigma=c(0.005, 0.0055, 0.006, 0.0065, 0.007))
model <- train(x,
               y,
               method = "svmRadial",
               metric = "ROC",
               family = binomial,
               #tuneLength = 5,
               tuneGrid = grid,
               trControl = trainControl(method = "repeatedcv",
                                        number = 10,
                                        repeats = 4,
                                        classProbs = TRUE,
                                        summaryFunction = twoClassSummary),
               preProc = c("center", "scale"),
               strata = y,
               sampsize = rep(40, 2))
model
#C = 0.5
#sigma = 0.006

# Final Train Predict
model <- train(x,
               y,
               method = "svmRadial",
               metric = "ROC",
               family = binomial,
               tuneGrid =expand.grid(C = 0.5, sigma=0.006),
               trControl = trainControl(method = "none",
                                        classProbs = TRUE,
                                        summaryFunction = twoClassSummary),
               preProc = c("center", "scale"),
               strata = y,
               sampsize = rep(40, 2))

# Final Predict
levels(Z.test$performance_score) <- c(40, 60, 80, 100)
Z.test$performance_score <- as.numeric(as.character(Z.test$performance_score))
dmy <- dummyVars(" ~ .", data = Z.test)
Z.test <- data.frame(predict(dmy, newdata = Z.test))
predictions <- predict(model, data.frame(X.test, Z.test))
final <- data.frame(final_results = predictions)
rownames(final) <- rownames(X.test)
write.csv(final, "final_results.csv", row.names = T)

library(pROC)
probsTrain <- predict(model, x, type = "prob")
rocCurve   <- roc(response = y,
                  predictor = probsTrain[, "D"],
                  levels = rev(levels(y)))
plot(rocCurve, print.thres = "best")

# best threshold 
probsTest <- predict(model, x, type = "prob")
threshold <- 0.768
pred <- factor(ifelse(probsTest[, "D"] > threshold, "D", "A"))
pred <- relevel(pred, "D")   # you may or may not need this; I did
confusionMatrix(pred, y)

# Predict
combined_test <- cbind(X.test, Z.test[,-which("performance_score" == colnames(Z.test))])
predictions <- predict(RF_model, newdata =  combined_test, type = "prob")