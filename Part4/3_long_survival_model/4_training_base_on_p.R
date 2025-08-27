#!/usr/bin/Rscript

rm(list = ls())
library(stringr)
library(gridExtra)
library(future)
library(sva)
library(neuralnet)
# library(e1071)
library(pROC)
library(ROCit)
library(caret)
library(doParallel)
library(cancerclass)
library(tidyverse)

setwd("~/help_for_others/DYQ/output/6_long_survival_model/")
load("./csv/CrNorm_normalized_training_data.rdata")
load("./csv/CrNorm_normalized_validation_data.rdata")
load("~/help_for_others/DYQ/output/6_long_survival_model_2/training_C1C6_delta_after_selection/training_feature_selection_based_on_p.rdata")

# sel = sel %>% arrange(desc(V1))
library(DMwR)
table(training$response)
set.seed(123)
training<-SMOTE(response~.,training,perc.over=200,perc.under=150)
table(training$response)


mm=2
sig_ls_final <- sapply(2:10,function(mm){
  if(length(which(lapply(sig_ls,length)==mm))>1){
    return(sig_ls[[which(lapply(sig_ls,length)==mm)[1]]])
  }else{
    return(sig_ls[[which(lapply(sig_ls,length)==mm)]])
  }
}) 


# sig_ls_final <- sig_ls_final[1:50]  

library(parallel)
library(doParallel)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

############建模与验证
mm=3
total_res <- lapply(1:9,function(mm){
  
  sig_this <- sig_ls_final[[mm]]

  training = training %>% dplyr::select(1,any_of(sig_this))
  validation = validation %>% dplyr::select(1,any_of(sig_this))

  
  Grid <- list( 
    svmLinear = expand.grid(C = c(0.01,0.1,1,3)),
    svmLinearWeights = expand.grid(cost = c(0.01,0.1,1,3), weight = c(0.1 ,0.5 ,1 ,2 ,3 ,5 ,10)),
    svmRadial = expand.grid(sigma = c(0.0005 ,0.001 ,0.005 ,0.01),C = c(0.01,0.1,1,3)),
    svmRadialWeights = expand.grid(sigma = c(0.0005 ,0.001 ,0.005 ,0.01),C = c(0.01,0.1,1,3), Weight = c(0.1 ,0.5 ,1 ,2 ,3 ,5 ,10)),
    lda = expand.grid(),
    glmnet = expand.grid(alpha = 1,lambda = seq(0.5, 1, by = 0.1)),
    bayesglm =  expand.grid()
  )
  TuneLength =  list( svmLinear = nrow(Grid[['svmLinear']]),
                      svmLinearWeights = nrow(Grid[['svmLinearWeights']]),
                      svmRadial = nrow(Grid[['svmRadial']]),
                      svmRadialWeights = nrow(Grid[['svmRadialWeights']]) ,
                      lda = 1,
                      glmnet = nrow(Grid[['glmnet']]),
                      bayesglm =  1
  )
  
  
  method = c('svmLinear','svmLinearWeights',
             'svmRadial','svmRadialWeights',
             'lda',
             'glmnet','bayesglm'
             )
  ls_model <- lapply(method,function(m){
    
    set.seed(10086)
    print(m)
    f = 5  # f folds resampling
    r = 5 # r repeats
    n = f*r
    
    # sets random seeds for parallel running for each single resampling f-folds and r-repeats cross-validation
    seeds <- vector(mode = "list", length = n + 1)
    #the number of tuning parameter
    for(i in 1:n) seeds[[i]] <- sample.int(n=10000, TuneLength[[m]])
    
    #for the last model
    seeds[[n+1]]<-sample.int(10000, 1)
    
    
    ctrl <- trainControl(method="repeatedcv",
                         number = f, ## 5-folds cv
                         summaryFunction=twoClassSummary,   # Use AUC to pick the best model
                         classProbs=TRUE,
                         repeats = r, ## 10-repeats cv,
                         seeds = seeds
    )
    
    
    if(m %in% c("lda","qda","bayesglm")){
      model.tune <- train(response ~ .,
                          data = training,
                          method = m,
                          metric="ROC",
                          trControl=ctrl)
    }else{
      model.tune <- train(response ~ .,
                          data = training,
                          method = m,
                          metric="ROC",
                          trControl=ctrl,
                          tuneGrid = Grid[[m]])
    }
    
    print(m)
    model.tune
    
    load("./csv/CrNorm_normalized_training_data.rdata")
    
    prob <- predict(model.tune,training[,-1],type = "prob")
    pre <- predict(model.tune,training[,-1])
    test_set <- data.frame(obs = as.factor(training$response), No = prob[,'No'], Yes = prob[,'Yes'], pred=pre)
    test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="No",0,1)})
    test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="No",0,1)})
    

    roc <- pROC::roc(as.numeric(test_set$obs),as.numeric(test_set$Yes),levels=c(0,1),direction="<")
    # roc_result <- coords(roc, "best")
    auc <- roc %>% auc()
    auc_ci <- ci(roc)
    
    accuracy <- sum(test_set$pred == test_set$obs) / length(test_set$obs)
    accuracy_ci <- binom.test(sum(test_set$pred == test_set$obs), length(test_set$obs))$conf.int
    
    test_set$obs <- factor(test_set$obs,levels = c(1,0))
    test_set$pred <- factor(test_set$pred,levels = c(1,0))
    
    conf_matrix <- confusionMatrix(data = test_set$pred, reference = test_set$obs)
    sensitivity <- conf_matrix$table["1","1"]/sum(conf_matrix$table[,"1"])
    sensitivity_ci <- GenBinomApps::clopper.pearson.ci(conf_matrix$table["1","1"],sum(conf_matrix$table[,"1"]),
                                                       CI = "two.sided", alpha = .05)
    
    specificity <- conf_matrix$table["0","0"]/sum(conf_matrix$table[,"0"])
    specificity_ci <- GenBinomApps::clopper.pearson.ci(conf_matrix$table["0","0"],sum(conf_matrix$table[,"0"]),
                                                       CI = "two.sided", alpha = .05)
    
    PPV <- conf_matrix$table["1","1"]/sum(conf_matrix$table["1",])
    PPV_ci <- GenBinomApps::clopper.pearson.ci(conf_matrix$table["1","1"],sum(conf_matrix$table["1",]),
                                               CI = "two.sided", alpha = .05)
    
    FPV <- conf_matrix$table["0","0"]/sum(conf_matrix$table["0",])
    FPV_ci <- GenBinomApps::clopper.pearson.ci(conf_matrix$table["0","0"],sum(conf_matrix$table["0",]),
                                               CI = "two.sided", alpha = .05)
    
    df_training <- data.frame(group="training",ROC=auc_ci[2], ROC_l=auc_ci[1],ROC_u=auc_ci[3],
                              accuracy=accuracy, acc_l=accuracy_ci[1], acc_u=accuracy_ci[2],
                              sensitivity=sensitivity, sen_l=as.numeric(sensitivity_ci[2]),sen_u=as.numeric(sensitivity_ci[3]),
                              specificity=specificity,spe_l=as.numeric(specificity_ci[2]),spe_u=as.numeric(specificity_ci[3]),
                              PPV=PPV,PPV_l=as.numeric(PPV_ci[2]),PPV_u=as.numeric(PPV_ci[3]),
                              FPV=FPV,FPV_l=as.numeric(FPV_ci[2]),FPV_u=as.numeric(FPV_ci[3]))
    
    
    
    
    
    
    
    prob <- predict(model.tune,validation[,-1],type = "prob")
    pre <- predict(model.tune,validation[,-1])
    test_set <- data.frame(obs = as.factor(validation$response), No = prob[,'No'], Yes = prob[,'Yes'], pred=pre)

    test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="No",0,1)})
    test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="No",0,1)})
    

    roc <- pROC::roc(as.numeric(test_set$obs),as.numeric(test_set$Yes),levels=c(0,1),direction="<")
    # roc_result <- coords(roc, "best")
    auc <- roc %>% auc()
    auc_ci <- ci(roc)
    
    accuracy <- sum(test_set$pred == test_set$obs) / length(test_set$obs)
    accuracy_ci <- binom.test(sum(test_set$pred == test_set$obs), length(test_set$obs))$conf.int
    
    test_set$obs <- factor(test_set$obs,levels = c(1,0))
    test_set$pred <- factor(test_set$pred,levels = c(1,0))
    
    conf_matrix <- confusionMatrix(data = test_set$pred, reference = test_set$obs)
    sensitivity <- conf_matrix$table["1","1"]/sum(conf_matrix$table[,"1"])
    sensitivity_ci <- GenBinomApps::clopper.pearson.ci(conf_matrix$table["1","1"],sum(conf_matrix$table[,"1"]),
                                                       CI = "two.sided", alpha = .05)
    
    specificity <- conf_matrix$table["0","0"]/sum(conf_matrix$table[,"0"])
    specificity_ci <- GenBinomApps::clopper.pearson.ci(conf_matrix$table["0","0"],sum(conf_matrix$table[,"0"]),
                                                       CI = "two.sided", alpha = .05)
    
    PPV <- conf_matrix$table["1","1"]/sum(conf_matrix$table["1",])
    PPV_ci <- GenBinomApps::clopper.pearson.ci(conf_matrix$table["1","1"],sum(conf_matrix$table["1",]),
                                               CI = "two.sided", alpha = .05)
    
    FPV <- conf_matrix$table["0","0"]/sum(conf_matrix$table["0",])
    FPV_ci <- GenBinomApps::clopper.pearson.ci(conf_matrix$table["0","0"],sum(conf_matrix$table["0",]),
                                               CI = "two.sided", alpha = .05)
    
    df_validation <- data.frame(group="validation",ROC=auc_ci[2], ROC_l=auc_ci[1],ROC_u=auc_ci[3],
                                accuracy=accuracy, acc_l=accuracy_ci[1], acc_u=accuracy_ci[2],
                                sensitivity=sensitivity, sen_l=as.numeric(sensitivity_ci[2]),sen_u=as.numeric(sensitivity_ci[3]),
                                specificity=specificity,spe_l=as.numeric(specificity_ci[2]),spe_u=as.numeric(specificity_ci[3]),
                                PPV=PPV,PPV_l=as.numeric(PPV_ci[2]),PPV_u=as.numeric(PPV_ci[3]),
                                FPV=FPV,FPV_l=as.numeric(FPV_ci[2]),FPV_u=as.numeric(FPV_ci[3]))
    
    
    df <- rbind(df_training,df_validation) %>% dplyr::mutate(method=m)%>% dplyr::mutate(num_of_mt=mm)
    
    return_ls <- list()
    return_ls[[1]]<-model.tune
    return_ls[[2]]<-df
    
    return(return_ls)
  } 

  ) 

  save(ls_model,file = paste("~/help_for_others/DYQ/output/6_long_survival_model_2/delta_model_res_base_on_p/model_res/",mm,".rdata",sep = ""))
  
  
})

