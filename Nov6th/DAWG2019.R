ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("doParallel", "randomForest", "mlbench", "caret", "e1071", "phyloseq", "ggplot")
ipak(packages)

## Random Forest classification on microbiome data##
# 1) Tune the mtry parameters using the caret package.
# 2) Run cforest using parameters from 1)
# 3) Calculate conditional variable importance for RF from 2).
# 4) Constructing AUC curve

#multicore processing with doParallel
cl <- makePSOCKcluster(4)
registerDoParallel(cl)

# 0) Preparing data for running random forest
# Since we are not using phyloseq object for random forest, we have to convert it into the seperate objects
# Also, using ASV table by itself takes too much time, so make it as Genus level.
ps.scaled <- scale(otu_table(ps.pruned), center = TRUE, scale = TRUE)  
ps.scaled[is.na(ps.scaled)] <- 0

map <- sample_data(ps.pruned)
tax <- tax_table(ps.pruned)
ps <- merge_phyloseq(sample_data(map), tax_table(tax), otu_table(ps.scaled, taxa_are_rows=FALSE))

genus_level <- taxa_level(ps, "Genus")
predictors <- as.data.frame(otu_table(genus_level))
dim(predictors)

response <- as.factor(sample_data(genus_level)$Sample)
rf.data <- data.frame(response, predictors)
rf.data[,1] # the fist column is defined as variable that we want to test

# 1) Tune the 'mtry' parameters using the caret pacakge
# Reference : https://machinelearningmastery.com/tune-machine-learning-algorithms-in-r/

control <- trainControl(method="repeatedcv", number=10, repeats=3)
mtry <- sqrt(ncol(predictors))
tunegrid <- expand.grid(.mtry=mtry)

# Algorithm Tune (tuneRF)
set.seed(7)
bestmtry <- tuneRF(predictors, response, stepFactor=1.5, improve=1e-5, ntree=1000)
print(bestmtry)

# optimal mtry is decided based on the OBB error rate#

# 2) Run random forest using parameters from 1)
set.seed(10)
samples <- sample(NROW(rf.data), NROW(rf.data) * .7)
data.train <- rf.data[samples, ]
data.test <- rf.data[-samples, ]

ctrl <- trainControl(method="repeatedcv", number=10, repeats=3,
                     summaryFunction=twoClassSummary, 
                     classProbs=T,
                     savePredictions = T)

forest.model <- train(response ~., data.train,method="rf",savePredictions = TRUE,ntree=1000,
                      prox=TRUE,allowParallel=TRUE, metric="ROC", trControl = ctrl,tuneGrid = data.frame(mtry = 51))

print(forest.model)
print(forest.model$finalModel)

result.predicted.prob <- predict(forest.model, data.test, type="prob") # Prediction

# 3) Calculate variable importance for RF from 2).
var_imp <- varImp(forest.model)

ImpMeasure<-data.frame(varImp(forest.model)$importance)
ImpMeasure$Vars<-row.names(ImpMeasure)
ImpMeasure <- ImpMeasure[order(-ImpMeasure$Overall),][1:15,]

ggplot(ImpMeasure, aes(x = reorder(Vars, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying samples\n with Sample type")

##AUC curve
set.seed(11)
result.roc <- roc(data.test$response, result.predicted.prob$stool) # Draw ROC curve.
print(result.roc)
plot(result.roc)




##Try to use it for the real data##
##Dividing data into two group, samples from stool and tracheal
##try to find if there is any successful prediction of mortaility in 1 month (variable "month")
##and what kind of genus of bacteria is highly associated with classification 

stool <- subset_samples(genus_level, Sample=="stool")
tracheal <- subset_samples(genus_level, Sample=="tracheal")

predictors_tracheal <- as.data.frame(otu_table(tracheal))
dim(predictors_tracheal)
response_tracheal <- as.factor(sample_data(tracheal)$month)
rf.data_tracheal_month <- data.frame(response_tracheal, predictors_tracheal)
rf.data_tracheal_month[,1] # the fist column is defined as variable that we want to test

set.seed(12)
bestmtry_tracheal_month <- tuneRF(predictors_tracheal, response_tracheal, stepFactor=1.5, improve=1e-5, ntree=1000)
print(bestmtry)

# seems like the error rate is huge.. Whatever#
set.seed(10)
samples <- sample(NROW(rf.data_tracheal_month), NROW(rf.data_tracheal_month) * .7)
data.train_tracheal <- rf.data_tracheal_month[samples, ]
data.test_tracheal <- rf.data_tracheal_month[-samples, ]

levels(data.train_tracheal$response_tracheal)[levels(data.train_tracheal$response_tracheal)=="0"] <- "zero"
levels(data.train_tracheal$response_tracheal)[levels(data.train_tracheal$response_tracheal)=="1"] <- "one"


ctrl <- trainControl(method="repeatedcv", number=10, repeats=3,
                     summaryFunction=twoClassSummary, 
                     classProbs=T,
                     savePredictions = T)

forest.model_tracheal_month <- train(response_tracheal~., data.train_tracheal,method="rf",savePredictions = TRUE,ntree=1000,
                      prox=TRUE,allowParallel=TRUE, metric="ROC", trControl = ctrl,tuneGrid = data.frame(mtry = 23))

print(forest.model_tracheal_month)
print(forest.model_tracheal_month$finalModel)

result.predicted.prob_tracheal <- predict(forest.model_tracheal_month, data.test_tracheal, type="prob") # Prediction

#Variable importance
var_imp <- varImp(forest.model_tracheal_month)

ImpMeasure<-data.frame(varImp(forest.model_tracheal_month)$importance)
ImpMeasure$Vars<-row.names(ImpMeasure)
ImpMeasure <- ImpMeasure[order(-ImpMeasure$Overall),][1:15,]

ggplot(ImpMeasure, aes(x = reorder(Vars, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying samples\n with Sample type")

##AUC curve
set.seed(11)
result.roc <- roc(data.test_tracheal$response_tracheal, result.predicted.prob_tracheal$zero) # Draw ROC curve.
print(result.roc)
plot(result.roc)


RF_flav_classify <- randomForest( x=esv_table_scaled_flav[,1:(ncol(esv_table_scaled_flav)-1)] , y=esv_table_scaled_flav[ , ncol(esv_table_scaled_flav)] , ntree=501, importance=TRUE, proximities=TRUE )

