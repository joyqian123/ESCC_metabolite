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
library(tidyverse)
setwd("~/help_for_others/DYQ/output//")

# select_method="bayesglm"
load("~/help_for_others/DYQ/output/6_long_survival_model_2/delta_model_res_base_on_p/model_res/8.rdata")
model = ls_model[[4]]


setwd("~/help_for_others/DYQ/output/6_long_survival_model/")
load("./csv/CrNorm_normalized_training_data.rdata")
load("./csv/CrNorm_normalized_validation_data.rdata")

model.tune = model[[1]]


setwd("~/help_for_others/DYQ/output/6_long_survival_model_2//")
###################################################################
###########################################################
library(caret)
library(colorspace)
# hcl_palettes(plot = TRUE)
cols <- diverge_hcl(5,"Blue-Red 2")[c(1,5)]
cols <- sequential_hcl(5,"Purples")[c(3)]
cols <- sequential_hcl(5,"RdPu")[c(4)]
draw_confusion_matrix <- function(cm,title_name) {
  layout(matrix(c(1,1,2,2)))
  par(mar=c(2,2,2,2))
  plot(c(120, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title(title_name, cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#5F8151')
  text(195, 435, 'NRnresponse', cex=1.2)
  rect(250, 430, 340, 370, col='#c8cfc4')
  text(295, 435, 'Response', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#c8cfc4')
  rect(250, 305, 340, 365, col='#5F8151')
  text(140, 400, 'NRnresponse', cex=1.2, srt=90)
  text(140, 335, 'Response', cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics
  plot(c(120, 345), c(0,110), type = "n", MAIN=NULL,xlab="", ylab="", xaxt='n', yaxt='s')
  rect(150, 0, 240, as.numeric(100*cm$byClass[3]), col='#A69DCC')
  rect(250, 0, 340, as.numeric(100*cm$byClass[4]), col='#FFC0C0')
  text(195, as.numeric(100*cm$byClass[3])/2, 
       scales::percent(round(cm$byClass[3],4),0.01), cex=1.6, font=2, col='white')
  text(295, as.numeric(100*cm$byClass[4])/2, 
       scales::percent(round(cm$byClass[4],4),0.01), cex=1.6, font=2, col='white')
  text(195, as.numeric(100*cm$byClass[3])+5,"Negative Predictive Value" , cex=1.2, font=2, col='black')
  text(295, as.numeric(100*cm$byClass[4])+5,"Positive Predictive Value" , cex=1.2, font=2, col='black')
} 





###################################################################
###################################################################
auc_plot <- function(rocobj){
  auc<-auc(rocobj)[1]
  auc_low<-pROC::ci(rocobj,of="auc")[1]
  auc_high<-pROC::ci(rocobj,of="auc")[3]
  ciobj <- ci.se(rocobj,specificities=seq(0, 1, 0.01))
  data_ci<-ciobj[1:101,1:3]
  data_ci<-as.data.frame(data_ci)
  x=as.numeric(rownames(data_ci))
  data_ci<-data.frame(x,data_ci)
  ggroc(rocobj,
        color="#D33F6A",
        size=1,
        legacy.axes = T 
  )+
    theme_classic()+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),       
                 colour='grey', 
                 linetype = 'dotdash'
    )+
    geom_text(aes(x=0.4,y=0.2,label=paste("AUC=",round(auc,3),", 95%CI: (",round(auc_low,3),","
                                          ,round(auc_high,3),")",sep = "")),size=4)
}


prob <- predict(model.tune,training[,-1],type = "prob")
pre <- predict(model.tune,training[,-1])
test_set <- data.frame(obs = as.factor(training$response), NR = prob[,'No'], R = prob[,'Yes'], pred=pre)
test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="No",0,1)})
test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="No",0,1)})

result_matrix<-confusionMatrix(table(test_set$pred, test_set$obs))
pdf("./plot/cm_training.pdf",width = 5,height = 7)
draw_confusion_matrix(result_matrix,"Training Cohort")
dev.off()




prob <- predict(model.tune,validation[,-1],type = "prob")
pre <- predict(model.tune,validation[,-1])
test_set <- data.frame(obs = as.factor(validation$response), NR = prob[,'No'], R = prob[,'Yes'], pred=pre)

test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="No",0,1)})
test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="No",0,1)})

pdf("./plot/cm_validation.pdf",width = 5,height = 7)
result_matrix<-confusionMatrix(table(test_set$pred, test_set$obs))
draw_confusion_matrix(result_matrix,"Validation Cohort")
dev.off()




#################################################################
prob <- predict(model.tune,training[,-1],type = "prob")
pre <- predict(model.tune,training[,-1])
test_set <- data.frame(obs = as.factor(training$response), NR = prob[,'No'], R = prob[,'Yes'], pred=pre)
test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="No",0,1)})
test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="No",0,1)})

pdf("./plot/roc_training.pdf",width = 4,height = 3)
rocobj<-roc(test_set$obs, test_set$R)
auc_plot(rocobj)
dev.off()



# 
load("../1_preclean/data/data_preclean.rdata")
clinic_need = clinic_in %>% dplyr::select(SUBJID,AVAL1,CNSR1,AVAL2,CNSR2) %>% distinct(.,.keep_all = T)
clinic_need$CNSR1 = ifelse(clinic_need$CNSR1==1,0,1)
clinic_need$CNSR2 = ifelse(clinic_need$CNSR2==1,0,1)
                           
                           
clinic_training = cbind(rownames(training),test_set)
colnames(clinic_training)[1]="Sample"
clinic_training$Sample = sapply(1:nrow(clinic_training),function(kk){
  if(nchar(clinic_training$Sample[kk])==4){
    return(paste("0",clinic_training$Sample[kk],sep = ""))
  }else{
    return(clinic_training$Sample[kk])
  }
})
clinic_training = clinic_training %>% 
  dplyr::rename(SUBJID=Sample) %>% inner_join(clinic_need[,c("SUBJID","AVAL1","CNSR1","AVAL2","CNSR2")],by = "SUBJID")

library(survminer) 
library(survival) 
data<-clinic_training
fit<-survfit(Surv(AVAL1, CNSR1) ~ pred, data = data)
summary(fit)
pdf("./plot/pfs_training.pdf",width = 5,height = 5)
ggsurvplot(fit, data = data,
           surv.median.line = "hv", 
           conf.int = FALSE, 
           risk.table = TRUE, 
           pval = TRUE,
           palette = "jco", 
           xlab = "Follow up time(m)", 
           legend = c(0.8,0.9), 
           legend.title = "", 
           legend.labs = c("High", "Low"), 
           xlim=c(0,35),
           break.x.by = 5)  
dev.off()


fit<-survfit(Surv(AVAL2, CNSR2) ~ pred, data = data)
summary(fit)
pdf("./plot/os_training.pdf",width = 5,height = 5)
ggsurvplot(fit, data = data,
           surv.median.line = "hv", 
           conf.int = FALSE, 
           risk.table = TRUE, 
           pval = TRUE,
           palette = "jco",  
           xlab = "Follow up time(m)", 
           legend = c(0.8,0.9), 
           legend.title = "", 
           legend.labs = c("High", "Low"), 
           xlim=c(0,40),
           break.x.by = 5)  
dev.off()

tt=20
pred_cutoff_train_PFS = lapply(6:24,function(tt){
  print(tt)
  tmp = data %>% dplyr::mutate(long=ifelse(AVAL1>tt,1,ifelse(CNSR1==1,0,NA))) %>% na.omit()
  library(pROC)
  roc <- pROC::roc(as.numeric(tmp$long),as.numeric(tmp$R),levels=c(0,1),direction="<")
  auc <- roc %>% auc()
  auc_ci <- ci(roc)
  best_cutoff <- pROC::coords(roc) %>% data.frame() %>% dplyr::mutate(value=specificity+sensitivity-1) %>% 
    arrange(desc(value))
  best_cf = best_cutoff$threshold[1]
  
  tmp$pred = ifelse(tmp$R>best_cf,1,0)
  
  accuracy <- sum(tmp$pred == tmp$long) / length(tmp$long)
  accuracy_ci <- binom.test(sum(tmp$pred == tmp$long), length(tmp$long))$conf.int
  
  return(data.frame(time_cf=tt,auc=auc,auc_l=auc_ci[1],auc_h=auc_ci[3],
                    accur=accuracy,accur_l=accuracy_ci[1],accur_h=accuracy_ci[2]))
}) %>% do.call(rbind,.) %>% dplyr::mutate(dataset="train")




tt=20
pred_cutoff_train_OS = lapply(18:36,function(tt){
  print(tt)
  tmp = data %>% dplyr::mutate(long=ifelse(AVAL2>tt,1,ifelse(CNSR2==1,0,NA))) %>% na.omit()
  library(pROC)
  roc <- pROC::roc(as.numeric(tmp$long),as.numeric(tmp$R),levels=c(0,1),direction="<")
  auc <- roc %>% auc()
  auc_ci <- ci(roc)
  best_cutoff <- pROC::coords(roc) %>% data.frame() %>% dplyr::mutate(value=specificity+sensitivity-1) %>% 
    arrange(desc(value))
  best_cf = best_cutoff$threshold[1]
  
  tmp$pred = ifelse(tmp$R>best_cf,1,0)
  
  accuracy <- sum(tmp$pred == tmp$long) / length(tmp$long)
  accuracy_ci <- binom.test(sum(tmp$pred == tmp$long), length(tmp$long))$conf.int
  
  return(data.frame(time_cf=tt,auc=auc,auc_l=auc_ci[1],auc_h=auc_ci[3],
                    accur=accuracy,accur_l=accuracy_ci[1],accur_h=accuracy_ci[2]))
}) %>% do.call(rbind,.) %>% dplyr::mutate(dataset="train")





#################################################################
prob <- predict(model.tune,validation[,-1],type = "prob")
pre <- predict(model.tune,validation[,-1])
test_set <- data.frame(obs = as.factor(validation$response), NR = prob[,'No'], R = prob[,'Yes'], pred=pre)
test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="No",0,1)})
test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="No",0,1)})

pdf("./plot/roc_validation.pdf",width = 4,height = 3)
rocobj<-roc(test_set$obs, test_set$R)
auc_plot(rocobj)
dev.off()




clinic_validation = cbind(rownames(validation),test_set)
colnames(clinic_validation)[1]="Sample"
clinic_validation$Sample = sapply(1:nrow(clinic_validation),function(kk){
  if(nchar(clinic_validation$Sample[kk])==4){
    return(paste("0",clinic_validation$Sample[kk],sep = ""))
  }else{
    return(clinic_validation$Sample[kk])
  }
})
clinic_validation = clinic_validation %>% 
  dplyr::rename(SUBJID=Sample) %>% inner_join(clinic_need[,c("SUBJID","AVAL1","CNSR1","AVAL2","CNSR2")],by = "SUBJID")

library(survminer) 
library(survival) 
data<-clinic_validation
fit<-survfit(Surv(AVAL1, CNSR1) ~ pred, data = data)
summary(fit)
pdf("./plot/pfs_validation.pdf",width = 5,height = 5)
ggsurvplot(fit, data = data,
           surv.median.line = "hv", 
           conf.int = FALSE, 
           risk.table = TRUE, 
           pval = TRUE,
           palette = "jco",  
           xlab = "Follow up time(m)", 
           legend = c(0.8,0.9), 
           legend.title = "", 
           legend.labs = c("High", "Low"), 
           xlim=c(0,35),
           break.x.by = 5)  
dev.off()


fit<-survfit(Surv(AVAL2, CNSR2) ~ pred, data = data)
summary(fit)
pdf("./plot/os_validation.pdf",width = 5,height = 5)
ggsurvplot(fit, data = data,
           surv.median.line = "hv", 
           conf.int = FALSE, 
           risk.table = TRUE, 
           pval = TRUE,
           palette = "jco", 
           xlab = "Follow up time(m)", 
           legend = c(0.8,0.9), 
           legend.title = "", 
           legend.labs = c("High", "Low"), 
           xlim=c(0,40),
           break.x.by = 5)  
dev.off()


tt=20
pred_cutoff_val_PFS = lapply(6:24,function(tt){
  print(tt)
  tmp = data %>% dplyr::mutate(long=ifelse(AVAL1>tt,1,ifelse(CNSR1==1,0,NA))) %>% na.omit()
  library(pROC)
  roc <- pROC::roc(as.numeric(tmp$long),as.numeric(tmp$R),levels=c(0,1),direction="<")
  # roc_result <- coords(roc, "best")
  auc <- roc %>% auc()
  auc_ci <- ci(roc)
  best_cutoff <- pROC::coords(roc) %>% data.frame() %>% dplyr::mutate(value=specificity+sensitivity-1) %>% 
    arrange(desc(value))
  best_cf = best_cutoff$threshold[1]
  
  tmp$pred = ifelse(tmp$R>best_cf,1,0)
  
  accuracy <- sum(tmp$pred == tmp$long) / length(tmp$long)
  accuracy_ci <- binom.test(sum(tmp$pred == tmp$long), length(tmp$long))$conf.int
  
  return(data.frame(time_cf=tt,auc=auc,auc_l=auc_ci[1],auc_h=auc_ci[3],
                    accur=accuracy,accur_l=accuracy_ci[1],accur_h=accuracy_ci[2]))
}) %>% do.call(rbind,.) %>% dplyr::mutate(dataset="validation")


tt=20
pred_cutoff_val_OS = lapply(18:36,function(tt){
  print(tt)
  tmp = data %>% dplyr::mutate(long=ifelse(AVAL2>tt,1,ifelse(CNSR2==1,0,NA))) %>% na.omit()
  library(pROC)
  roc <- pROC::roc(as.numeric(tmp$long),as.numeric(tmp$R),levels=c(0,1),direction="<")
  # roc_result <- coords(roc, "best")
  auc <- roc %>% auc()
  auc_ci <- ci(roc)
  best_cutoff <- pROC::coords(roc) %>% data.frame() %>% dplyr::mutate(value=specificity+sensitivity-1) %>% 
    arrange(desc(value))
  best_cf = best_cutoff$threshold[1]
  
  tmp$pred = ifelse(tmp$R>best_cf,1,0)
  
  accuracy <- sum(tmp$pred == tmp$long) / length(tmp$long)
  accuracy_ci <- binom.test(sum(tmp$pred == tmp$long), length(tmp$long))$conf.int
  
  return(data.frame(time_cf=tt,auc=auc,auc_l=auc_ci[1],auc_h=auc_ci[3],
                    accur=accuracy,accur_l=accuracy_ci[1],accur_h=accuracy_ci[2]))
}) %>% do.call(rbind,.) %>% dplyr::mutate(dataset="validation")




PFS = rbind(pred_cutoff_train_PFS,pred_cutoff_val_PFS)

pdf("./plot/PFS_accuracy_month_cf.pdf",width = 5,height = 3)
ggplot()+
  geom_point(data = pred_cutoff_train_PFS,aes(x=as.numeric(time_cf)-0.15,y=accur,col=dataset))+
  geom_line(data = pred_cutoff_train_PFS,aes(x=as.numeric(time_cf)-0.15,y=accur,col=dataset))+
  geom_point(data = pred_cutoff_val_PFS,aes(x=as.numeric(time_cf)+0.15,y=accur,col=dataset))+
  geom_line(data = pred_cutoff_val_PFS,aes(x=as.numeric(time_cf)+0.15,y=accur,col=dataset))+
  scale_color_manual(values = c("#C71585","#4B0082"))+
  geom_hline(yintercept = 0.6,lwd=0.5,lty=2,col="black")+
  xlab("time_cf")+
  ylab("accuracy")+
  scale_x_continuous(breaks = seq(6,24,1),labels = unique(pred_cutoff_train_PFS$time_cf))+
  theme_bw()
dev.off()



pdf("./plot/OS_accuracy_month_cf.pdf",width = 5,height = 3)
ggplot()+
  geom_point(data = pred_cutoff_train_OS,aes(x=as.numeric(time_cf)-0.15,y=accur,col=dataset))+
  geom_line(data = pred_cutoff_train_OS,aes(x=as.numeric(time_cf)-0.15,y=accur,col=dataset))+
  geom_point(data = pred_cutoff_val_OS,aes(x=as.numeric(time_cf)+0.15,y=accur,col=dataset))+
  geom_line(data = pred_cutoff_val_OS,aes(x=as.numeric(time_cf)+0.15,y=accur,col=dataset))+
  scale_color_manual(values = c("#C71585","#4B0082"))+
  geom_hline(yintercept = 0.6,lwd=0.5,lty=2,col="black")+
  xlab("time_cf")+
  ylab("accuracy")+
  # xlim(15,30)+
  scale_x_continuous(breaks = seq(18,36,1),labels = unique(pred_cutoff_train_OS$time_cf))+
  theme_bw()
dev.off()




pdf("./plot/PFS_AUCy_month_cf.pdf",width = 5,height = 3)
ggplot()+
  geom_point(data = pred_cutoff_train_PFS,aes(x=as.numeric(time_cf)-0.15,y=auc,col=dataset))+
  geom_line(data = pred_cutoff_train_PFS,aes(x=as.numeric(time_cf)-0.15,y=auc,col=dataset))+
  geom_point(data = pred_cutoff_val_PFS,aes(x=as.numeric(time_cf)+0.15,y=auc,col=dataset))+
  geom_line(data = pred_cutoff_val_PFS,aes(x=as.numeric(time_cf)+0.15,y=auc,col=dataset))+
  scale_color_manual(values = c("#C71585","#4B0082"))+
  geom_hline(yintercept = 0.6,lwd=0.5,lty=2,col="black")+
  xlab("time_cf")+
  ylab("AUCy")+
  scale_x_continuous(breaks = seq(6,24,1),labels = unique(pred_cutoff_train_PFS$time_cf))+
  theme_bw()
dev.off()



pdf("./plot/OS_AUCy_month_cf.pdf",width = 5,height = 3)
ggplot()+
  geom_point(data = pred_cutoff_train_OS,aes(x=as.numeric(time_cf)-0.15,y=auc,col=dataset))+
  geom_line(data = pred_cutoff_train_OS,aes(x=as.numeric(time_cf)-0.15,y=auc,col=dataset))+
  geom_point(data = pred_cutoff_val_OS,aes(x=as.numeric(time_cf)+0.15,y=auc,col=dataset))+
  geom_line(data = pred_cutoff_val_OS,aes(x=as.numeric(time_cf)+0.15,y=auc,col=dataset))+
  scale_color_manual(values = c("#C71585","#4B0082"))+
  geom_hline(yintercept = 0.6,lwd=0.5,lty=2,col="black")+
  xlab("time_cf")+
  ylab("AUCy")+
  scale_x_continuous(breaks = seq(18,36,1),labels = unique(pred_cutoff_train_OS$time_cf))+
  theme_bw()
dev.off()



output = rbind(pred_cutoff_train_OS,pred_cutoff_val_OS,pred_cutoff_train_PFS,pred_cutoff_val_PFS)
output = output %>% dplyr::mutate(auc_1 = str_c(round(auc,3)," (",round(auc_l,3),",",round(auc_h,3),")")) %>% 
  dplyr::mutate(accuracy_1 = str_c(round(accur,3)," (",round(accur_l,3),",",round(accur_h,3),")")) %>% 
  dplyr::select(time_cf,dataset,auc_1,accuracy_1) %>% 
  dplyr::rename(auc=auc_1,accuracy=accuracy_1)


write.csv(output,file = "./plot/output.csv")

