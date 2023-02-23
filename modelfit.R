########  Details ########
### Author: Marco Hamins-Puertolas 
### Code to fit xgboost model and subsequently calculate household related metrics
### Reads in anonymized patient data and outputs prediction of infection for each interval
### as well as a household specific dataset that contains information on attack rates
### and infection history within the household from previous intervals

######## Load necessary libraries and create useful functions ########
library(xgboost)
library(tidyverse)
library(caret)
library(data.table)

# function to create importance matrix 
impMatAddVar <- function(ogMat,newMat,fitnum){
  for(i in 1:dim(newMat)[1]){
    loc <- which((ogMat$Feature == newMat$Feature[i])&(ogMat$model == fitnum))
    ogMat$Gain[loc] <- ogMat$Gain[loc] +  newMat$Gain[i]
    ogMat$Cover[loc] <- ogMat$Cover[loc] +  newMat$Cover[i]
    ogMat$Frequency[loc] <- ogMat$Frequency[loc] +  newMat$Frequency[i]
  }
  return(ogMat)
}


######## Set up analysis  ########
curDir <- getwd()

# Read in data
load(file=paste0(curDir,"/Data/housedfCohortData.rData"))

# Subset for training
householdCohortData$prediction <- 0
householdCohortData <- householdCohortData[order(householdCohortData$trainInf,decreasing = F),]
householdCohortDataForFit <- householdCohortData %>% dplyr::select(c("je_pre","d1_pre","d2_pre","d3_pre","d4_pre","max_pre","min_pre","avg_pre",
                                                                     "var_pre","je_post","d1_post","d2_post","d3_post","d4_post","max_post","min_post","var_post",  
                                                                     "avg_post","diff_je","diff_d1","diff_d2","diff_d3","diff_d4","max_diff","min_diff","avg_diff",  
                                                                     "var_diff","div_je","div_d1","div_d2","div_d3","div_d4","max_div","min_div","avg_div",
                                                                     "var_div","timeDiff","sampDate_pre_year","sampDate_post_year","ageAtFollow","ageAtEnr","trainInf","training"))
householdCohortDataTraining <- householdCohortDataForFit[householdCohortDataForFit$training==1,]
householdCohortDataTraining <- householdCohortDataTraining[,-ncol(householdCohortDataTraining)]
householdCohortDataTesting <- householdCohortDataForFit[!householdCohortDataForFit$training==1,]
householdCohortDataTesting <- householdCohortDataTesting[,-ncol(householdCohortDataTesting)]



##### Set up xgboost and run to find best fitting models #####
# XGBoost fit parameters
totatConfirmed <- 90
set.seed(12) # set for reproducibility
totiter <- 500 # total number of iterations run per outer k fold
num2avg <- 100 # total # of best iterations saved
outerkfold <- 4 # number of outer folds for nested k-fold CV
outerCV <- c(sample(rep(1:outerkfold,length=dim(householdCohortDataTraining)[1]-totatConfirmed),replace=F),
             sample(rep(1:outerkfold,length=totatConfirmed),replace=F)) # randomly assign k fold to both +s and -s
innerkfold <- 5 # number of inner folds for nested k-fold CV
relPosWeight <- 30 # relative importance of each positive sample in comparison to each negative

# Create places to save information about the best models
best_param = list()
best_seednumber = rep(1234,num2avg)
best_logloss = rep(Inf,num2avg)
best_logloss_index =  rep(0,num2avg)
best_logloss_k =  rep(0,num2avg)

# Create Dmatrix from training data
df <- xgb.DMatrix(data=data.matrix(householdCohortDataTraining[,-ncol(householdCohortDataTraining)]),
                  label = as.numeric(data.matrix(householdCohortDataTraining[,ncol(householdCohortDataTraining)])),
                  weight=householdCohortDataTraining$trainInf*(relPosWeight-1)+1)

# Run XGBoost -- runs through and finds the best parameter combinations according to logloss
for(k in 1:outerkfold){ # iterate through outer folds
  kfold_householdCohortDataTraining_train <- householdCohortDataTraining[outerCV != k,]
  kfold_householdCohortDataTraining_test <- householdCohortDataTraining[outerCV == k,]
  
  for (iter in 1:totiter) { # run through total iterations
    param <- list(objective = "multi:softprob",
                  eval_metric = c("mlogloss","auc"),
                  num_class = 2,
                  max_depth = sample(3:10, 1),
                  eta = runif(1, .01, .3),
                  gamma = runif(1, 0.0, 0.2), 
                  subsample = runif(1, .6, .9),
                  colsample_bytree = runif(1, .5, .8), 
                  min_child_weight = sample(1:40, 1),
                  max_delta_step = sample(1:10, 1)
    )
    cv.nround = 1000
    cv.nfold = innerkfold
    seed.number = sample.int(10000, 1)[[1]]
    set.seed(seed.number)
    
    # Create downsampled models
    df_down <- downSample(x = kfold_householdCohortDataTraining_train[, -ncol(kfold_householdCohortDataTraining_train)],
                          y = as.factor(kfold_householdCohortDataTraining_train$trainInf))
    df_down_mat <- xgb.DMatrix(data=data.matrix(df_down[,-ncol(kfold_householdCohortDataTraining_train)]),
                               label = as.numeric(data.matrix(df_down[,ncol(kfold_householdCohortDataTraining_train)])),
                               weight=(df_down$trainInf*(relPosWeight-1)+1))
    
    mdcv <- xgb.cv(data=df_down_mat, params = param, nthread=6, 
                   nfold=cv.nfold, nrounds=cv.nround,stratified = TRUE,
                   verbose = F, early_stopping_rounds=8, maximize=FALSE)
    
    min_eval_metric = min(mdcv$evaluation_log$test_mlogloss_mean)
    min_eval_metric_index = which.min(mdcv$evaluation_log$test_mlogloss_mean)
    
    if (sum(min_eval_metric < best_logloss)>0) { # save best models
      maxloc <- min(which(best_logloss == max(best_logloss)))
      best_logloss[maxloc] = min_eval_metric
      best_logloss_index[maxloc] = min_eval_metric_index
      best_seednumber[maxloc] = seed.number
      best_param[[maxloc]] = param
      best_logloss_k[[maxloc]] = k
    }
  }
  
}


##### Get predictions for the testing data ######

df2 <- xgb.DMatrix(data=data.matrix(householdCohortDataTesting[,-ncol(householdCohortDataTesting)]),
                   label = as.numeric(data.matrix(householdCohortDataTesting[,ncol(householdCohortDataTesting)])))
pred2 <- rep(0,dim(df2)[1]*2)

# creating empty parameter and importance mat values 
md <- xgb.train(data=df2, params=param, nrounds=50, nthread=6)
importanceMat <- data.table(Feature = rep(md$feature_names,num2avg),
                            Gain = 0,
                            Cover = 0,
                            Frequency = 0,
                            model = rep(1:num2avg,each=length(md$feature_names)))

for(bestiter in 1:num2avg){
  # train on best models
  nround = best_logloss_index[bestiter]
  set.seed(best_seednumber[bestiter])
  df_down <- downSample(x = householdCohortDataTraining[, -ncol(householdCohortDataTraining)],
                        y = as.factor(householdCohortDataTraining$trainInf))
  df_down_mat <- xgb.DMatrix(data=data.matrix(df_down[,-ncol(householdCohortDataTraining)]),
                             label = as.numeric(data.matrix(df_down[,ncol(householdCohortDataTraining)])),
                             weight=df_down$trainInf*(relPosWeight-1)+1)
  
  md <- xgb.train(data=df_down_mat, params=best_param[[bestiter]], nrounds=nround, nthread=6)
  
  pred2 <- predict(md,df2,"prob") + pred2
  
  # save importance matrix
  tempImp <- xgb.importance(model=md)
  importanceMat <- impMatAddVar(importanceMat,tempImp,bestiter)
}

pred2 <- matrix(pred2/num2avg,nrow=dim(df2)[1],ncol=2,byrow = T)
householdCohortData[!householdCohortData$training==1,]$prediction <- pred2[,2]


##### Get predictions for the training data ######

# create a 10 fold prediction model: use only on the training data
outerCV <- c(sample(rep(1:10,length=dim(householdCohortDataTraining)[1]-90),replace=F),
             sample(rep(1:10,length=90),replace=F))
pred2 <- rep(0,dim(householdCohortDataTraining)[1])

for(bestiter in 1:num2avg){
  for(k in 1:10){
    # train on best model
    nround <- best_logloss_index[bestiter]
    set.seed(best_seednumber[bestiter])
    kfold_tempdata_train <- householdCohortDataTraining[outerCV != k,]
    kfold_tempdata_test <- householdCohortDataTraining[outerCV == k,]
    
    df <- xgb.DMatrix(data=data.matrix(householdCohortDataTraining[,-ncol(householdCohortDataTraining)]),
                      label = as.numeric(data.matrix(householdCohortDataTraining[,ncol(householdCohortDataTraining)])),
                      weight=householdCohortDataTraining$trainInf*(relPosWeight-1)+1)
    
    df_down <- downSample(x = kfold_tempdata_train[, -ncol(kfold_tempdata_train)],
                          y = as.factor(kfold_tempdata_train$trainInf))
    df_down_mat <- xgb.DMatrix(data=data.matrix(df_down[,-ncol(kfold_tempdata_train)]),
                               label = as.numeric(data.matrix(df_down[,ncol(kfold_tempdata_train)])),
                               weight=df_down$trainInf*(relPosWeight-1)+1)
    df_test <- xgb.DMatrix(data=data.matrix(kfold_tempdata_test[,-ncol(kfold_tempdata_test)]),label = as.numeric(data.matrix(kfold_tempdata_test[,ncol(kfold_tempdata_test)])))
    
    md <- xgb.train(data=df_down_mat, params=best_param[[bestiter]], nrounds=nround, nthread=6)
    
    # get predictions
    pred2[outerCV==k] <- predict(md,df_test,"prob")[1:length(pred2[outerCV == k])*2] + pred2[outerCV == k]
  }
}

householdCohortDataTraining$prediction <- pred2/num2avg
householdCohortData[householdCohortData$training==1,]$prediction <- householdCohortDataTraining$prediction
householdCohortData <- householdCohortData %>% rowwise() %>% mutate(prediction_round = ifelse(training,trainInf,round(prediction)==1))


###### Save/export ######
save(householdCohortData,file=paste0(curDir,"/Data/housedfCohortData_update.rData"))

save(importanceMat,file = paste0(curDir,"/Data/importanceMat_update.rData"))


