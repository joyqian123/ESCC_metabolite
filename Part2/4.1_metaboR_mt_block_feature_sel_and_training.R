#!/usr/bin/Rscript
rm(list = ls())
library(stringr)
library(gridExtra)
library(future)
library(sva)
# library(e1071)
library(pROC)
library(ROCit)
library(caret)
library(doParallel)
library(cancerclass)

dir.create("~/help_for_others/DYQ/output/3_baseline_model/training_mt_block/")  ######
setwd("~/help_for_others/DYQ/output/3_baseline_model/training_mt_block/")
library(MetaboAnalystR)
library(tidyverse)
mSet<-InitDataObjects("conc", "stat", FALSE)


mSet<-Read.TextData(mSet, "~/help_for_others/DYQ/output/3_baseline_model/csv/block_training.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)



mSet$dataSet$filt <- mSet$dataSet$edit <- NULL
preproc <- qs::qread("preproc.qs")
int.mat <- preproc
mSet$dataSet$proc.feat.num <- ncol(int.mat)
qs::qsave(as.data.frame(int.mat), file = "data_proc.qs")
mSet$msgSet$replace.msg <- paste("Zero or missing values were replaced by 1/5 of the min positive value for each variable.")
invisible(gc())

mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)

mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###

training = cbind(mSet$dataSet$cls,mSet$dataSet$norm) %>% dplyr::rename(response=colnames(.)[1])
save(training,file = "../csv/df_mt_block_training_norm.rdata")







#########################
#########################
#####################fold change analysis
# Perform fold-change analysis on uploaded data, unpaired
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)  

# Plot fold-change analysis
mSet<-PlotFC(mSet, "fc_0_", "pdf", 72, width=NA)

# To view fold-change 
mSet$analSet$fc$fc.log

FC_tb = mSet$analSet$fc$fc.log %>% data.frame() %>% rownames_to_column("mt")
colnames(FC_tb)[2]="FC"


#####################T-Test
# Perform T-test (parametric)
mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, pvalType = "fdr", TRUE)

# Plot of the T-test results
mSet<-PlotTT(mSet, imgName = "tt_0_", format = "pdf", dpi = 72, width=NA)
t_tb = mSet$analSet$tt$p.value %>% data.frame() %>% rownames_to_column("mt")
colnames(t_tb)[2]="p"


FC_tb = FC_tb %>% inner_join(t_tb,by = "mt")


# #####################Volcano Plot
# Perform the volcano analysis
mSet<-Volcano.Anal(mSet, paired=FALSE, fcthresh = 2.0, cmpType = 0, nonpar = F, threshp = 0.1,
                   equal.var = TRUE, pval.type = "fdr")
# Create the volcano plot
# mSet<-PlotVolcano(mSet, "volcano_0_", 1, 0, format ="pdf", dpi=72, width=NA)

# 
# 
####################Orthogonal Partial Least Squares - Discriminant Analysis (orthoPLS-DA)
# Perform oPLS-DA analysis
mSet<-OPLSR.Anal(mSet, reg=TRUE)

# Create a 2D oPLS-DA score plot
mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", format = "pdf", dpi=72, width=NA, inx1 = 1,inx2 = 2,
                      reg = 0.95,show = 1,grey.scale = 0)
mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", plotType = "all", format = "pdf", dpi = 72, width=NA)
mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", format = "pdf", dpi = 72, width=NA, type = "vip", feat.nm = "tscore", 
                   feat.num = 15,color.BW = FALSE)
mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", format = "pdf", dpi=72, width=NA)
mSet<-OPLSDA.Permut(mSet, 100)
mSet<-PlotOPLS.Permutation(mSet, "opls_perm_2_", format = "pdf", dpi=72, width=NA)


VIP_tb <- mSet$analSet$oplsda$vip.mat %>% as.data.frame(check.names=FALSE) %>% rownames_to_column("mt") %>% dplyr::select(-V2)



sel <- FC_tb %>% inner_join(VIP_tb,by = "mt")
ggplot(sel)+
  geom_point(aes(x=V1, y = FC))

save(sel,file = "./training_sel.rdata")






########################################################################
#######################################################################
#######################################################################
sel_bl = sel %>% dplyr::filter(V1>0.5) %>% pull(mt)
norm = mSet$dataSet$norm %>% dplyr::select(any_of(sel_bl))
block_dt_train_sel = data.frame(response=mSet$dataSet$cls,norm) %>% rownames_to_column("Sample")

write.csv(block_dt_train_sel,file = "../csv/block_training_after_selection.csv",row.names = FALSE)



dir.create("~/help_for_others/DYQ/output/3_baseline_model/training_mt_block_after_selection/")  #######
setwd("~/help_for_others/DYQ/output/3_baseline_model/training_mt_block_after_selection//")
library(MetaboAnalystR)
library(tidyverse)
mSet<-InitDataObjects("conc", "stat", FALSE)


mSet<-Read.TextData(mSet, "~/help_for_others/DYQ/output/3_baseline_model/csv/block_training_after_selection.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)


mSet$dataSet$filt <- mSet$dataSet$edit <- NULL
preproc <- qs::qread("preproc.qs")
int.mat <- preproc
mSet$dataSet$proc.feat.num <- ncol(int.mat)
qs::qsave(as.data.frame(int.mat), file = "data_proc.qs")
mSet$msgSet$replace.msg <- paste("Zero or missing values were replaced by 1/5 of the min positive value for each variable.")
invisible(gc())

mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)

mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###




#########################
#########################
#####################fold change analysis
# Perform fold-change analysis on uploaded data, unpaired
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)  

# Plot fold-change analysis
mSet<-PlotFC(mSet, "fc_0_", "pdf", 72, width=NA)

# To view fold-change 
mSet$analSet$fc$fc.log

FC_tb = mSet$analSet$fc$fc.log %>% data.frame() %>% rownames_to_column("mt")
colnames(FC_tb)[2]="FC"


#####################T-Test
# Perform T-test (parametric)
mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, pvalType = "fdr", TRUE)

# Plot of the T-test results
mSet<-PlotTT(mSet, imgName = "tt_0_", format = "pdf", dpi = 72, width=NA)
t_tb = mSet$analSet$tt$p.value %>% data.frame() %>% rownames_to_column("mt")
colnames(t_tb)[2]="p"


FC_tb = FC_tb %>% inner_join(t_tb,by = "mt")


# #####################Volcano Plot
# Perform the volcano analysis
mSet<-Volcano.Anal(mSet, paired=FALSE, fcthresh = 2.0, cmpType = 0, nonpar = F, threshp = 0.1,
                   equal.var = TRUE, pval.type = "fdr")
# Create the volcano plot
# mSet<-PlotVolcano(mSet, "volcano_0_", 1, 0, format ="pdf", dpi=72, width=NA)

# 
# 
####################Orthogonal Partial Least Squares - Discriminant Analysis (orthoPLS-DA)
# Perform oPLS-DA analysis
mSet<-OPLSR.Anal(mSet, reg=TRUE)

# Create a 2D oPLS-DA score plot
mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", format = "pdf", dpi=72, width=NA, inx1 = 1,inx2 = 2,
                      reg = 0.95,show = 1,grey.scale = 0)
mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", plotType = "all", format = "pdf", dpi = 72, width=NA)
mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", format = "pdf", dpi = 72, width=NA, type = "vip", feat.nm = "tscore", 
                   feat.num = 15,color.BW = FALSE)
mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", format = "pdf", dpi=72, width=NA)
mSet<-OPLSDA.Permut(mSet, 100)
mSet<-PlotOPLS.Permutation(mSet, "opls_perm_2_", format = "pdf", dpi=72, width=NA)


VIP_tb <- mSet$analSet$oplsda$vip.mat %>% as.data.frame(check.names=FALSE) %>% rownames_to_column("mt") %>% dplyr::select(-V2)













######################
######################
load("~/help_for_others/DYQ/output/3_baseline_model/csv/df_mt_block_training_norm.rdata")
load("~/help_for_others/DYQ/output/3_baseline_model/csv/df_mt_block_validation_norm.rdata")

load("~/help_for_others/DYQ/output/3_baseline_model/training_mt_block//training_sel.rdata")
i=2
sig_ls <- lapply(2:100,function(i){
  print(i)
  for(t in 2:i){
    if(t==2){
      sel_sig_1 = sel %>% arrange(dplyr::desc(V1)) %>% .[1:t,] %>% pull(mt)
      training_try <- training[,colnames(training) %in% c(sel_sig_1)]
      cor_sel=cor(training_try)
      if(table(cor_sel>0.8)["TRUE"]==ncol(cor_sel)){
        sel_sig = sel_sig_1
      }else{
        sel_sig = sel_sig_1[-t]
      }
    }else{
      sel_sig_add <- sel %>% arrange(dplyr::desc(V1)) %>% .[t,] %>% pull(mt)
      sel_sig_1 <- c(sel_sig,sel_sig_add) %>% unlist()
      training_try <- training[,colnames(training) %in% c(sel_sig_1)]
      cor_sel=cor(training_try)
      print(table(cor_sel>0.8)["TRUE"]==ncol(cor_sel))
      if(table(cor_sel>0.8)["TRUE"]==ncol(cor_sel)){
        sel_sig = sel_sig_1
      }else{
        sel_sig = sel_sig
      }
    }
  }
  return(sel_sig)
})


sig_ls_final <- sapply(2:20,function(mm){
  if(length(which(lapply(sig_ls,length)==mm))>1){
    return(sig_ls[[which(lapply(sig_ls,length)==mm)[1]]])
  }else{
    return(sig_ls[[which(lapply(sig_ls,length)==mm)]])
  }
}) 


# sig_ls_final <- sig_ls_final[1:50]  
save(sig_ls_final,file = "~/help_for_others/DYQ/output/3_baseline_model/mt_block_model_res/sig_ls_final.rdata")








############
library(parallel)
library(doParallel)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

total_res <- lapply(1:19,function(mm){
  
  sig_this <- sig_ls_final[[mm]]
  
  training <- training[,colnames(training) %in% c('response', sig_this)]
  validation  <- validation[,colnames(validation) %in% c('response',sig_this)]
  
  Grid <- list( #nb = expand.grid(fL =  c(0,0.5,1,1.5,2.0), usekernel = TRUE, adjust = c(0.5,0.75,1,1.25,1.5)),
    svmLinear = expand.grid(C = c( 1 ,3 ,5 ,10 ,20)),
    svmLinearWeights = expand.grid(cost = c( 1 ,3 ,5 ,10 ,20), weight = c(0.1 ,0.5 ,1 ,2 ,3 ,5 ,10)),
    svmRadial = expand.grid(sigma = c(0.0005 ,0.001 ,0.005 ,0.01 ,0.05),C = c( 1 ,3 ,5 ,10 ,20)),
    svmRadialWeights = expand.grid(sigma = c(0.0005 ,0.001 ,0.005 ,0.01 ,0.05),C = c( 1 ,3 ,5 ,10 ,20), Weight = c(0.1 ,0.5 ,1 ,2 ,3 ,5 ,10)),
    lda = expand.grid(),
    glmnet = expand.grid(alpha = seq(0.1, 1, by = 0.1),lambda = seq(0.1, 1, by = 0.1)),
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
  ##model training with different algorithms
  ls_model <- lapply(method,function(m){
    set.seed(10086)
  
    print(m)
    f = 5  # f folds resampling
    r = 10 # r repeats
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
    
    
    prob <- predict(model.tune,training[,-1],type = "prob")
    pre <- predict(model.tune,training[,-1])
    test_set <- data.frame(obs = as.factor(training$response), NR = prob[,'NR'], R = prob[,'R'], pred=pre)
    test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="NR",0,1)})
    test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="NR",0,1)})
    
    # class(test_set$pred)
    # auc_training <- twoClassSummary(test_set, lev = levels(test_set$obs))
    # auc_training
    
    roc <- pROC::roc(as.numeric(test_set$obs),as.numeric(pre),levels=c(0,1),direction="<")
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
    test_set <- data.frame(obs = as.factor(validation$response), NR = prob[,'NR'], R = prob[,'R'], pred=pre)
    # class(test_set$pred)
    # auc_validation <- twoClassSummary(test_set, lev = levels(test_set$obs))
    # auc_validation
    
    test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="NR",0,1)})
    test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="NR",0,1)})
    
    # class(test_set$pred)
    # auc_training <- twoClassSummary(test_set, lev = levels(test_set$obs))
    # auc_training
    
    roc <- pROC::roc(as.numeric(test_set$obs),as.numeric(pre),levels=c(0,1),direction="<")
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
    
    
    df <- rbind(df_training,df_validation) %>% dplyr::mutate(method=m) %>% dplyr::mutate(num_of_mt=mm)
    
    return_ls <- list()
    return_ls[[1]]<-model.tune
    return_ls[[2]]<-df
    
    return(return_ls)
  } 
  
  
  ) 
  
  
  save(ls_model,file = paste("~/help_for_others/DYQ/output/3_baseline_model/mt_block_model_res/model_res/",mm,".rdata",sep = ""))
  
  
})




