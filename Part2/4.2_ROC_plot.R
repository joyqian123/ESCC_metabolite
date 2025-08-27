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
library(readxl)
library(readr)
setwd("~/help_for_others/DYQ/output//")


kk=15
load(paste("~/help_for_others/DYQ/output/3_baseline_model/mt_block_model_res/model_res/",kk,".rdata",sep = ""))

model.tune = ls_model[[3]][[1]]


setwd("~/help_for_others/DYQ/output/3_baseline_model//")
load("./csv/df_mt_block_training_norm.rdata")
load("./csv/df_mt_block_validation_norm.rdata")

sig_this <- model.tune$trainingData %>% colnames(.) %>% .[-1]

training = training %>% dplyr::select(1,any_of(sig_this))
validation = validation %>% dplyr::select(1,any_of(sig_this))
colnames(training)[1]="response"
colnames(validation)[1]="response"











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
  # 计算AUC值
  auc<-auc(rocobj)[1]
  # AUC的置信区间
  auc_low<-pROC::ci(rocobj,of="auc")[1]
  auc_high<-pROC::ci(rocobj,of="auc")[3]
  # 计算置信区间
  ciobj <- ci.se(rocobj,specificities=seq(0, 1, 0.01))
  data_ci<-ciobj[1:101,1:3]
  data_ci<-as.data.frame(data_ci)
  x=as.numeric(rownames(data_ci))
  data_ci<-data.frame(x,data_ci)
  ggroc(rocobj,
        color="#D33F6A",
        size=1,
        legacy.axes = T #
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
test_set <- data.frame(obs = as.factor(training$response), NR = prob[,'NR'], R = prob[,'R'], pred=pre)
test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="NR",0,1)})
test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="NR",0,1)})

result_matrix<-confusionMatrix(table(test_set$pred, test_set$obs))
pdf("./plot/cm_training.pdf",width = 5,height = 7)
draw_confusion_matrix(result_matrix,"Training Cohort")
dev.off()




prob <- predict(model.tune,validation[,-1],type = "prob")
pre <- predict(model.tune,validation[,-1])
test_set <- data.frame(obs = as.factor(validation$response), NR = prob[,'NR'], R = prob[,'R'], pred=pre)

test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="NR",0,1)})
test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="NR",0,1)})

pdf("./plot/cm_validation.pdf",width = 5,height = 7)
result_matrix<-confusionMatrix(table(test_set$pred, test_set$obs))
draw_confusion_matrix(result_matrix,"Validation Cohort")
dev.off()










#################################################################
prob <- predict(model.tune,training[,-1],type = "prob")
pre <- predict(model.tune,training[,-1])
test_set <- data.frame(obs = as.factor(training$response), NR = prob[,'NR'], R = prob[,'R'], pred=pre)
test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="NR",0,1)})
test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="NR",0,1)})

pdf("./plot/roc_training.pdf",width = 4,height = 3)
rocobj<-roc(test_set$obs, test_set$R)
auc_plot(rocobj)
dev.off()

# 
load("../1_preclean/data/data_preclean.rdata")
clinic_1 = clinic_in %>% dplyr::filter(Sample %in% rownames(training)) %>% 
  dplyr::select(SUBJID,6,7,9,10,11:13,15:18,20:21,R) %>% dplyr::select(SUBJID,R,everything()) %>% 
  distinct(SUBJID,.keep_all = T) 
multi_roc = lapply(3:ncol(clinic_1),function(kk){
  tmp = clinic_1 %>% dplyr::select(R,kk) %>% dplyr::mutate(R=ifelse(R=="R",1,0)) %>% 
    dplyr::rename(pred=colnames(.)[2]) %>% dplyr::mutate(pred1=as.numeric(as.factor(pred)))
  if(colnames(clinic_1)[kk]=="AGEGROUP"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###<65——R
  }
  if(colnames(clinic_1)[kk]=="SEX"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###女——R
  }
  if(colnames(clinic_1)[kk]=="ETHNIC"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###汉族——R
  }
  if(colnames(clinic_1)[kk]=="BLWTGR1"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###>=60 kg——R
  }
  if(colnames(clinic_1)[kk]=="BLBMGR1"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###<20 kg/m^2——R
  }
  if(colnames(clinic_1)[kk]=="BLECOG"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###1——R
  }
  if(colnames(clinic_1)[kk]=="LMYN"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###否——R
  }
  if(colnames(clinic_1)[kk]=="PDL1GR1"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###PDL1>=1%——R
  }
  if(colnames(clinic_1)[kk]=="PDL1GR2"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))   ###PDL1>=5%——R
  }
  if(colnames(clinic_1)[kk]=="PDL1GR3"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))   ###PDL1>=10%——R
  }
  if(colnames(clinic_1)[kk]=="MHTMENUM"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))   ###1——R
  }
  if(colnames(clinic_1)[kk]=="SUSMFREQ"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))   ###1——R
  }
  if(colnames(clinic_1)[kk]=="SUDRFREQ"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))   ###1——R
  }
  return(r)
})





names(multi_roc)=colnames(clinic_1)[3:ncol(clinic_1)]
multi_roc[[length(multi_roc)+1]] = rocobj
names(multi_roc)[length(multi_roc)] = "metabolites_score"
# ggroc(multi_roc)+
#   theme_bw()
col_merge = c("#478bbe","#e6966a","#7c6aa4","#83c5db","#aedacb","#aeb0d8","#d75821","#f7bfa4",
              "#e6b119","#b5c574","#b65741","#4eb69e","lightgrey","#D33F6A")

pdf("./plot/multi_auc_training.pdf",width = 5,height = 4)
ggroc(multi_roc,legacy.axes = T)+ #
  theme_classic()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),        #
               colour='grey', 
               linetype = 'dotdash'
  )+
  scale_color_manual(values = col_merge)+
  annotate("text",x=0.75,y=0.52,size=4,col=col_merge[1],label=paste(names(multi_roc)[1],"_AUC = ",round(multi_roc[[1]]$auc,3)))+
  annotate("text",x=0.75,y=0.48,size=4,col=col_merge[2],label=paste(names(multi_roc)[2],"_AUC = ",round(multi_roc[[2]]$auc,3)))+
  annotate("text",x=0.75,y=0.44,size=4,col=col_merge[3],label=paste(names(multi_roc)[3],"_AUC = ",round(multi_roc[[3]]$auc,3)))+
  annotate("text",x=0.75,y=0.40,size=4,col=col_merge[4],label=paste(names(multi_roc)[4],"_AUC = ",round(multi_roc[[4]]$auc,3)))+
  annotate("text",x=0.75,y=0.36,size=4,col=col_merge[5],label=paste(names(multi_roc)[5],"_AUC = ",round(multi_roc[[5]]$auc,3)))+
  annotate("text",x=0.75,y=0.32,size=4,col=col_merge[6],label=paste(names(multi_roc)[6],"_AUC = ",round(multi_roc[[6]]$auc,3)))+
  annotate("text",x=0.75,y=0.28,size=4,col=col_merge[7],label=paste(names(multi_roc)[7],"_AUC = ",round(multi_roc[[7]]$auc,3)))+
  annotate("text",x=0.75,y=0.24,size=4,col=col_merge[8],label=paste(names(multi_roc)[8],"_AUC = ",round(multi_roc[[8]]$auc,3)))+
  annotate("text",x=0.75,y=0.20,size=4,col=col_merge[9],label=paste(names(multi_roc)[9],"_AUC = ",round(multi_roc[[9]]$auc,3)))+
  annotate("text",x=0.75,y=0.16,size=4,col=col_merge[10],label=paste(names(multi_roc)[10],"_AUC = ",round(multi_roc[[10]]$auc,3)))+
  annotate("text",x=0.75,y=0.12,size=4,col=col_merge[11],label=paste(names(multi_roc)[11],"_AUC = ",round(multi_roc[[11]]$auc,3)))+
  annotate("text",x=0.75,y=0.08,size=4,col=col_merge[12],label=paste(names(multi_roc)[12],"_AUC = ",round(multi_roc[[12]]$auc,3)))+
  annotate("text",x=0.75,y=0.04,size=4,col=col_merge[13],label=paste(names(multi_roc)[13],"_AUC = ",round(multi_roc[[13]]$auc,3)))+
  annotate("text",x=0.75,y=0.56,size=4,col=col_merge[14],label=paste(names(multi_roc)[14],"_AUC = ",round(multi_roc[[14]]$auc,3)))
dev.off()




clinic_need = clinic_in %>% dplyr::select(SUBJID,AVAL1,CNSR1,AVAL2,CNSR2) %>% distinct(.,.keep_all = T)
clinic_need$CNSR1 = ifelse(clinic_need$CNSR1==1,0,1)
clinic_need$CNSR2 = ifelse(clinic_need$CNSR2==1,0,1)
                           
                           
clinic_training = cbind(rownames(training),test_set)
colnames(clinic_training)[1]="Sample"
clinic_training$Sample = str_remove(clinic_training$Sample,"C[0-9]D[0-9]")
clinic_training$Sample = sapply(1:nrow(clinic_training),function(kk){
  if(nchar(clinic_training$Sample[kk])==4){
    return(paste("0",clinic_training$Sample[kk],sep = ""))
  }else{
    return(clinic_training$Sample[kk])
  }
})
clinic_training = clinic_training %>% 
  dplyr::rename(SUBJID=Sample) %>% inner_join(clinic_need[,c("SUBJID","AVAL1","CNSR1","AVAL2","CNSR2")],by = "SUBJID")

library(survminer) #
library(survival) #
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
           break.x.by = 5))
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


roc <- pROC::roc(as.numeric(test_set$obs),as.numeric(test_set$R),levels=c(0,1),direction="<")
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

df_train <- data.frame(group="training",ROC=auc_ci[2], ROC_l=auc_ci[1],ROC_u=auc_ci[3],
                      accuracy=accuracy, acc_l=accuracy_ci[1], acc_u=accuracy_ci[2],
                      sensitivity=sensitivity, sen_l=as.numeric(sensitivity_ci[2]),sen_u=as.numeric(sensitivity_ci[3]),
                      specificity=specificity,spe_l=as.numeric(specificity_ci[2]),spe_u=as.numeric(specificity_ci[3]),
                      PPV=PPV,PPV_l=as.numeric(PPV_ci[2]),PPV_u=as.numeric(PPV_ci[3]),
                      FPV=FPV,FPV_l=as.numeric(FPV_ci[2]),FPV_u=as.numeric(FPV_ci[3]))



clinic_training_1 = clinic_training %>% dplyr::select(SUBJID,R,6:9) %>% dplyr::rename(metabo_score=R) %>% 
  inner_join(clinic_1,by = "SUBJID") %>% dplyr::select(-R) %>% 
  dplyr::select(1,3:6,everything())



##############################################
PFS_tb = clinic_training_1 %>% dplyr::select(AVAL1,CNSR1,6:19)
PFS_tb$SUSMFREQ = factor(PFS_tb$SUSMFREQ,levels = c("从不吸烟","已经戒烟","正在吸烟"),ordered = T) %>% as.numeric()
table(PFS_tb$SUSMFREQ)
PFS_tb$SUDRFREQ = factor(PFS_tb$SUDRFREQ,levels = c("从不饮酒","已经戒酒","正在饮酒"),ordered = T) %>% as.numeric()
table(PFS_tb$SUDRFREQ)

multicox <- coxph(Surv(AVAL1,CNSR1) ~ ., data = PFS_tb) 
multisum <- summary(multicox)

item <- colnames(PFS_tb)[3:ncol(PFS_tb)]
HR <- multisum$coefficients[,2]
L95CI <- multisum$conf.int[,3]
H95CI <- multisum$conf.int[,4]
pvalue <- multisum$coefficients[,5]
multiresult <- data.frame(item=item,
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue)

multiresult = multiresult %>% dplyr::mutate(type="multi-cox")

uniresult <- data.frame()  
for(i in colnames(PFS_tb[,3:ncol(PFS_tb)])){   
  unicox <- coxph(Surv(AVAL1,CNSR1) ~ PFS_tb[,i], data = PFS_tb)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
    uniresult <- rbind(uniresult,
                       cbind(item=i,
                             HR=unisum$coefficients[,2],
                             L95CI=unisum$conf.int[,3],
                             H95CI=unisum$conf.int[,4],
                             pvalue=unisum$coefficients[,5]
                       ))
}   

uniresult = uniresult %>% dplyr::mutate(type="uni-cox")

plot = multiresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
p1=ggplot(plot)+
  geom_point(aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("PFS_multivariant-Cox")+
  scale_x_log10()+
  theme_bw()

plot = uniresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
plot$HR = as.numeric(plot$HR)
plot$L95CI = as.numeric(plot$L95CI)
plot$H95CI = as.numeric(plot$H95CI)
p2=ggplot(plot)+
  geom_point(aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("PFS_univariant-Cox")+
  scale_x_log10()+
  theme_bw()

##############################################
OS_tb = clinic_training_1 %>% dplyr::select(AVAL2,CNSR2,6:19)
OS_tb$SUSMFREQ = factor(OS_tb$SUSMFREQ,levels = c("从不吸烟","已经戒烟","正在吸烟"),ordered = T) %>% as.numeric()
table(OS_tb$SUSMFREQ)
OS_tb$SUDRFREQ = factor(OS_tb$SUDRFREQ,levels = c("从不饮酒","已经戒酒","正在饮酒"),ordered = T) %>% as.numeric()
table(OS_tb$SUDRFREQ)

multicox <- coxph(Surv(AVAL2,CNSR2) ~ ., data = OS_tb) 
multisum <- summary(multicox)

item <- colnames(OS_tb)[3:ncol(OS_tb)]
HR <- multisum$coefficients[,2]
L95CI <- multisum$conf.int[,3]
H95CI <- multisum$conf.int[,4]
pvalue <- multisum$coefficients[,5]
multiresult <- data.frame(item=item,
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue)

multiresult = multiresult %>% dplyr::mutate(type="multi-cox")

uniresult <- data.frame()  
for(i in colnames(OS_tb[,3:ncol(OS_tb)])){   
  unicox <- coxph(Surv(AVAL2,CNSR2) ~ OS_tb[,i], data = OS_tb)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  uniresult <- rbind(uniresult,
                     cbind(item=i,
                           HR=unisum$coefficients[,2],
                           L95CI=unisum$conf.int[,3],
                           H95CI=unisum$conf.int[,4],
                           pvalue=unisum$coefficients[,5]
                     ))
}   

uniresult = uniresult %>% dplyr::mutate(type="uni-cox")

plot = multiresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
p3=ggplot(plot)+
  geom_point(aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("OS_multivariant-Cox")+
  scale_x_log10()+
  theme_bw()

plot = uniresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
plot$HR = as.numeric(plot$HR)
plot$L95CI = as.numeric(plot$L95CI)
plot$H95CI = as.numeric(plot$H95CI)
p4=ggplot(plot)+
  geom_point(aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("OS_univariant-Cox")+
  scale_x_log10()+
  theme_bw()


pdf("./plot/cox_coef_training.pdf",width = 10,height = 8)
cowplot::plot_grid(p1,p2,p3,p4,nrow = 2)
dev.off()









#################################################################
prob <- predict(model.tune,validation[,-1],type = "prob")
pre <- predict(model.tune,validation[,-1])
test_set <- data.frame(obs = as.factor(validation$response), NR = prob[,'NR'], R = prob[,'R'], pred=pre)
test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="NR",0,1)})
test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="NR",0,1)})

pdf("./plot/roc_validation.pdf",width = 4,height = 3)
rocobj<-roc(test_set$obs, test_set$R)
auc_plot(rocobj)
dev.off()



load("../1_preclean/data/data_preclean.rdata")
clinic_1 = clinic_in %>% dplyr::filter(Sample %in% rownames(validation)) %>% 
  dplyr::select(SUBJID,6,7,9,10,11:13,15:18,20:21,R) %>% dplyr::select(SUBJID,R,everything()) %>% 
  distinct(SUBJID,.keep_all = T) 
multi_roc = lapply(3:ncol(clinic_1),function(kk){
  tmp = clinic_1 %>% dplyr::select(R,kk) %>% dplyr::mutate(R=ifelse(R=="R",1,0)) %>% 
    dplyr::rename(pred=colnames(.)[2]) %>% dplyr::mutate(pred1=as.numeric(as.factor(pred)))
  if(colnames(clinic_1)[kk]=="AGEGROUP"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###<65——R
  }
  if(colnames(clinic_1)[kk]=="SEX"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###女——R
  }
  if(colnames(clinic_1)[kk]=="ETHNIC"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###汉族——R
  }
  if(colnames(clinic_1)[kk]=="BLWTGR1"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###>=60 kg——R
  }
  if(colnames(clinic_1)[kk]=="BLBMGR1"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###<20 kg/m^2——R
  }
  if(colnames(clinic_1)[kk]=="BLECOG"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###1——R
  }
  if(colnames(clinic_1)[kk]=="LMYN"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###否——R
  }
  if(colnames(clinic_1)[kk]=="PDL1GR1"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###PDL1>=1%——R
  }
  if(colnames(clinic_1)[kk]=="PDL1GR2"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))   ###PDL1>=5%——R
  }
  if(colnames(clinic_1)[kk]=="PDL1GR3"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))   ###PDL1>=10%——R
  }
  if(colnames(clinic_1)[kk]=="MHTMENUM"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))   ###1——R
  }
  if(colnames(clinic_1)[kk]=="SUSMFREQ"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))   ###1——R
  }
  if(colnames(clinic_1)[kk]=="SUDRFREQ"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))   ###1——R
  }
  return(r)
})

names(multi_roc)=colnames(clinic_1)[3:ncol(clinic_1)]
multi_roc[[length(multi_roc)+1]] = rocobj
names(multi_roc)[length(multi_roc)] = "metabolites_score"
# ggroc(multi_roc)+
#   theme_bw()
col_merge = c("#478bbe","#e6966a","#7c6aa4","#83c5db","#aedacb","#aeb0d8","#d75821","#f7bfa4",
              "#e6b119","#b5c574","#b65741","#4eb69e","lightgrey","#D33F6A")

pdf("./plot/multi_auc_validation.pdf",width = 5,height = 4)
ggroc(multi_roc,legacy.axes = T)+ # FALSE时 横坐标为1-0 specificity；TRUE时 横坐标为0-1 1-specificity
  theme_classic()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),        # 绘制对角线
               colour='grey', 
               linetype = 'dotdash'
  )+
  scale_color_manual(values = col_merge)+
  annotate("text",x=0.75,y=0.52,size=4,col=col_merge[1],label=paste(names(multi_roc)[1],"_AUC = ",round(multi_roc[[1]]$auc,3)))+
  annotate("text",x=0.75,y=0.48,size=4,col=col_merge[2],label=paste(names(multi_roc)[2],"_AUC = ",round(multi_roc[[2]]$auc,3)))+
  annotate("text",x=0.75,y=0.44,size=4,col=col_merge[3],label=paste(names(multi_roc)[3],"_AUC = ",round(multi_roc[[3]]$auc,3)))+
  annotate("text",x=0.75,y=0.40,size=4,col=col_merge[4],label=paste(names(multi_roc)[4],"_AUC = ",round(multi_roc[[4]]$auc,3)))+
  annotate("text",x=0.75,y=0.36,size=4,col=col_merge[5],label=paste(names(multi_roc)[5],"_AUC = ",round(multi_roc[[5]]$auc,3)))+
  annotate("text",x=0.75,y=0.32,size=4,col=col_merge[6],label=paste(names(multi_roc)[6],"_AUC = ",round(multi_roc[[6]]$auc,3)))+
  annotate("text",x=0.75,y=0.28,size=4,col=col_merge[7],label=paste(names(multi_roc)[7],"_AUC = ",round(multi_roc[[7]]$auc,3)))+
  annotate("text",x=0.75,y=0.24,size=4,col=col_merge[8],label=paste(names(multi_roc)[8],"_AUC = ",round(multi_roc[[8]]$auc,3)))+
  annotate("text",x=0.75,y=0.20,size=4,col=col_merge[9],label=paste(names(multi_roc)[9],"_AUC = ",round(multi_roc[[9]]$auc,3)))+
  annotate("text",x=0.75,y=0.16,size=4,col=col_merge[10],label=paste(names(multi_roc)[10],"_AUC = ",round(multi_roc[[10]]$auc,3)))+
  annotate("text",x=0.75,y=0.12,size=4,col=col_merge[11],label=paste(names(multi_roc)[11],"_AUC = ",round(multi_roc[[11]]$auc,3)))+
  annotate("text",x=0.75,y=0.08,size=4,col=col_merge[12],label=paste(names(multi_roc)[12],"_AUC = ",round(multi_roc[[12]]$auc,3)))+
  annotate("text",x=0.75,y=0.04,size=4,col=col_merge[13],label=paste(names(multi_roc)[13],"_AUC = ",round(multi_roc[[13]]$auc,3)))+
  annotate("text",x=0.75,y=0.56,size=4,col=col_merge[14],label=paste(names(multi_roc)[14],"_AUC = ",round(multi_roc[[14]]$auc,3)))
dev.off()



clinic_validation = cbind(rownames(validation),test_set)
colnames(clinic_validation)[1]="Sample"
clinic_validation$Sample = str_remove(clinic_validation$Sample,"C[0-9]D[0-9]")
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


roc <- pROC::roc(as.numeric(test_set$obs),as.numeric(test_set$R),levels=c(0,1),direction="<")
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

df_vali <- data.frame(group="validation",ROC=auc_ci[2], ROC_l=auc_ci[1],ROC_u=auc_ci[3],
                       accuracy=accuracy, acc_l=accuracy_ci[1], acc_u=accuracy_ci[2],
                       sensitivity=sensitivity, sen_l=as.numeric(sensitivity_ci[2]),sen_u=as.numeric(sensitivity_ci[3]),
                       specificity=specificity,spe_l=as.numeric(specificity_ci[2]),spe_u=as.numeric(specificity_ci[3]),
                       PPV=PPV,PPV_l=as.numeric(PPV_ci[2]),PPV_u=as.numeric(PPV_ci[3]),
                       FPV=FPV,FPV_l=as.numeric(FPV_ci[2]),FPV_u=as.numeric(FPV_ci[3]))



clinic_validation_1 = clinic_validation %>% dplyr::select(SUBJID,R,6:9) %>% dplyr::rename(metabo_score=R) %>% 
  inner_join(clinic_1,by = "SUBJID") %>% dplyr::select(-R) %>% 
  dplyr::select(1,3:6,everything()) %>% na.omit()





##############################################
PFS_tb = clinic_validation_1 %>% dplyr::select(AVAL1,CNSR1,6:19)
PFS_tb$SUSMFREQ = factor(PFS_tb$SUSMFREQ,levels = c("从不吸烟","已经戒烟","正在吸烟"),ordered = T) %>% as.numeric()
table(PFS_tb$SUSMFREQ)
PFS_tb$SUDRFREQ = factor(PFS_tb$SUDRFREQ,levels = c("从不饮酒","已经戒酒","正在饮酒"),ordered = T) %>% as.numeric()
table(PFS_tb$SUDRFREQ)

colnames(PFS_tb)
library(coxphf)
multicox <- coxphf(Surv(AVAL1,CNSR1) ~ metabo_score + AGEGROUP + SEX + ETHNIC + BLWTGR1 + BLBMGR1 + BLECOG + LMYN + PDL1GR1 + PDL1GR2 + PDL1GR3 + MHTMENUM + SUSMFREQ + SUDRFREQ, data = PFS_tb) 
multisum <- summary(multicox)

item <- colnames(PFS_tb)[3:ncol(PFS_tb)]
HR <- multisum$coefficients %>% exp()
L95CI <- multisum$ci.lower
H95CI <- multisum$ci.upper
pvalue <- multisum$prob
multiresult <- data.frame(item=item,
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue)
multiresult$pvalue[is.na(multiresult$pvalue)]=1

multiresult = multiresult %>% dplyr::mutate(type="multi-cox")

uniresult <- data.frame()  
for(i in colnames(PFS_tb[,3:ncol(PFS_tb)])){   
  unicox <- coxph(Surv(AVAL1,CNSR1) ~ PFS_tb[,i], data = PFS_tb)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  uniresult <- rbind(uniresult,
                     cbind(item=i,
                           HR=unisum$coefficients[,2],
                           L95CI=unisum$conf.int[,3],
                           H95CI=unisum$conf.int[,4],
                           pvalue=unisum$coefficients[,5]
                     ))
}   

uniresult = uniresult %>% dplyr::mutate(type="uni-cox")

plot = multiresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
p1=ggplot()+
  geom_point(data = plot,aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=0.01,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=10,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("PFS_multivariant-Cox")+
  scale_x_log10()+
  theme_bw()

plot = uniresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
plot$HR = as.numeric(plot$HR)
plot$L95CI = as.numeric(plot$L95CI)
plot$H95CI = as.numeric(plot$H95CI)
p2=ggplot(plot)+
  geom_point(data = plot,aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=0.01,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=10,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("PFS_univariant-Cox")+
  scale_x_log10()+
  theme_bw()





##############################################
OS_tb = clinic_validation_1 %>% dplyr::select(AVAL2,CNSR2,6:19)
OS_tb$SUSMFREQ = factor(OS_tb$SUSMFREQ,levels = c("从不吸烟","已经戒烟","正在吸烟"),ordered = T) %>% as.numeric()
table(OS_tb$SUSMFREQ)
OS_tb$SUDRFREQ = factor(OS_tb$SUDRFREQ,levels = c("从不饮酒","已经戒酒","正在饮酒"),ordered = T) %>% as.numeric()
table(OS_tb$SUDRFREQ)

multicox <- coxph(Surv(AVAL2,CNSR2) ~ ., data = OS_tb) 
multisum <- summary(multicox)

item <- colnames(OS_tb)[3:ncol(OS_tb)]
HR <- multisum$coefficients[,2]
L95CI <- multisum$conf.int[,3]
H95CI <- multisum$conf.int[,4]
pvalue <- multisum$coefficients[,5]
multiresult <- data.frame(item=item,
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue)

multiresult = multiresult %>% dplyr::mutate(type="multi-cox")

uniresult <- data.frame()  
for(i in colnames(OS_tb[,3:ncol(OS_tb)])){   
  unicox <- coxph(Surv(AVAL2,CNSR2) ~ OS_tb[,i], data = OS_tb)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  uniresult <- rbind(uniresult,
                     cbind(item=i,
                           HR=unisum$coefficients[,2],
                           L95CI=unisum$conf.int[,3],
                           H95CI=unisum$conf.int[,4],
                           pvalue=unisum$coefficients[,5]
                     ))
}   

uniresult = uniresult %>% dplyr::mutate(type="uni-cox")

plot = multiresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
p3=ggplot(plot)+
  geom_point(data = plot,aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=0.01,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=10,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("OS_multivariant-Cox")+
  scale_x_log10()+
  theme_bw()

plot = uniresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
plot$HR = as.numeric(plot$HR)
plot$L95CI = as.numeric(plot$L95CI)
plot$H95CI = as.numeric(plot$H95CI)
p4=ggplot(plot)+
  geom_point(data = plot,aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=0.01,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=10,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("OS_univariant-Cox")+
  scale_x_log10()+
  theme_bw()


pdf("./plot/cox_coef_validation.pdf",width = 10,height = 8)
cowplot::plot_grid(p1,p2,p3,p4,nrow = 2)
dev.off()
















# #####################################################
# #####################################################
load("~/help_for_others/DYQ/output/3_baseline_model/csv/df_mt_block_chemo_norm.rdata")
load("~/help_for_others/DYQ/output/3_baseline_model/csv/df_mt_block_surg_norm.rdata")
load("~/help_for_others/DYQ/output/3_baseline_model/csv/df_mt_block_rad_norm.rdata")


prob <- predict(model.tune,chemo[,-1],type = "prob")
pre <- predict(model.tune,chemo[,-1])
test_set <- data.frame(obs = as.factor(chemo$response), NR = prob[,'NR'], R = prob[,'R'], pred=pre)
test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="NR",0,1)})
test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="NR",0,1)})

person_pred = cbind(rownames(chemo),test_set)
save(person_pred,file = "../7_JM/pred_unknown/R_NR_pred_res.rdata")


result_matrix<-confusionMatrix(table(test_set$pred, test_set$obs))
pdf("./plot/cm_chemo.pdf",width = 5,height = 7)
draw_confusion_matrix(result_matrix,"Test Cohort")
dev.off()

pdf("./plot/roc_chemo.pdf",width = 4,height = 3)
rocobj<-roc(test_set$obs, test_set$R)
auc_plot(rocobj)
dev.off()


chemo_test = cbind(rownames(chemo),test_set)
colnames(chemo_test)[1]="Sample"

chemo_test = chemo_test %>% dplyr::mutate(t_or_f = ifelse(obs==pred,"T","F"))
table(chemo_test$t_or_f)
load("~/help_for_others/DYQ//output/1_preclean/validation_chemo//data_preclean.rdata")
chemo_test = chemo_test %>% inner_join(clinic_in[,c("Sample","PFS","status","OS","status2")],by = "Sample")


clinic_info = read_excel("../../data/validation_data/final_clinical_data.xlsx",sheet=4)
clinic_1 = clinic_in[,c("SUBJID","Sample","R")] %>% dplyr::filter(Sample %in% rownames(chemo)) %>% 
  inner_join(clinic_info,by = "SUBJID") %>% dplyr::select(-Sample) %>% 
  dplyr::select(SUBJID,R,everything()) %>% 
  distinct(SUBJID,.keep_all = T) 
colnames(clinic_1)[10:11]=c("SUSMFREQ","SUDRFREQ")
multi_roc = lapply(3:ncol(clinic_1),function(kk){
  tmp = clinic_1 %>% dplyr::select(R,kk) %>% dplyr::mutate(R=ifelse(R=="R",1,0)) %>% 
    dplyr::rename(pred=colnames(.)[2]) %>% dplyr::mutate(pred1=as.numeric(as.factor(pred)))
  if(colnames(clinic_1)[kk]=="AGEGROUP"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###<65——R
  }
  if(colnames(clinic_1)[kk]=="SEX"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###女——R
  }
  if(colnames(clinic_1)[kk]=="ETHNIC"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###汉族——R
  }
  if(colnames(clinic_1)[kk]=="BLWTGR1"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###>=60 kg——R
  }
  if(colnames(clinic_1)[kk]=="BLBMGR1"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###<20 kg/m^2——R
  }
  if(colnames(clinic_1)[kk]=="BLECOG"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))  ###1——R
  }
  if(colnames(clinic_1)[kk]=="LMYN"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###否——R
  }
  if(colnames(clinic_1)[kk]=="PDL1GR1"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))  ###PDL1>=1%——R
  }
  if(colnames(clinic_1)[kk]=="PDL1GR2"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))   ###PDL1>=5%——R
  }
  if(colnames(clinic_1)[kk]=="PDL1GR3"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))   ###PDL1>=10%——R
  }
  if(colnames(clinic_1)[kk]=="MHTMENUM"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(0,1))   ###1——R
  }
  if(colnames(clinic_1)[kk]=="SUSMFREQ"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))   ###1——R
  }
  if(colnames(clinic_1)[kk]=="SUDRFREQ"){
    r = pROC::roc(as.numeric(tmp$R),as.numeric(tmp$pred1),levels=c(1,0))   ###1——R
  }
  return(r)
})

names(multi_roc)=colnames(clinic_1)[3:ncol(clinic_1)]
multi_roc[[length(multi_roc)+1]] = rocobj
names(multi_roc)[length(multi_roc)] = "metabolites_score"
# ggroc(multi_roc)+
#   theme_bw()
col_merge = c("#e6966a","#478bbe","#83c5db","#aedacb","#aeb0d8","#d75821",
              # "#f7bfa4","#e6b119","#b5c574",
              "#b65741","#4eb69e","lightgrey","#D33F6A")

pdf("./plot/multi_auc_chemo.pdf",width = 5,height = 4)
ggroc(multi_roc,legacy.axes = T)+ 
  theme_classic()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),        
               colour='grey', 
               linetype = 'dotdash'
  )+
  scale_color_manual(values = col_merge)+
  annotate("text",x=0.75,y=0.52,size=4,col=col_merge[1],label=paste(names(multi_roc)[1],"_AUC = ",round(multi_roc[[1]]$auc,3)))+
  annotate("text",x=0.75,y=0.48,size=4,col=col_merge[2],label=paste(names(multi_roc)[2],"_AUC = ",round(multi_roc[[2]]$auc,3)))+
  annotate("text",x=0.75,y=0.44,size=4,col=col_merge[3],label=paste(names(multi_roc)[3],"_AUC = ",round(multi_roc[[3]]$auc,3)))+
  annotate("text",x=0.75,y=0.40,size=4,col=col_merge[4],label=paste(names(multi_roc)[4],"_AUC = ",round(multi_roc[[4]]$auc,3)))+
  annotate("text",x=0.75,y=0.36,size=4,col=col_merge[5],label=paste(names(multi_roc)[5],"_AUC = ",round(multi_roc[[5]]$auc,3)))+
  annotate("text",x=0.75,y=0.32,size=4,col=col_merge[6],label=paste(names(multi_roc)[6],"_AUC = ",round(multi_roc[[6]]$auc,3)))+
  annotate("text",x=0.75,y=0.28,size=4,col=col_merge[7],label=paste(names(multi_roc)[7],"_AUC = ",round(multi_roc[[7]]$auc,3)))+
  annotate("text",x=0.75,y=0.24,size=4,col=col_merge[8],label=paste(names(multi_roc)[8],"_AUC = ",round(multi_roc[[8]]$auc,3)))+
  annotate("text",x=0.75,y=0.20,size=4,col=col_merge[9],label=paste(names(multi_roc)[9],"_AUC = ",round(multi_roc[[9]]$auc,3)))+
  annotate("text",x=0.75,y=0.16,size=4,col=col_merge[10],label=paste(names(multi_roc)[10],"_AUC = ",round(multi_roc[[10]]$auc,3)))
dev.off()





fit<-survfit(Surv(PFS, status) ~ pred, data = chemo_test)
summary(fit)
pdf("./plot/pfs_chemo.pdf",width = 5,height = 5)
ggsurvplot(fit, data = chemo_test,
           surv.median.line = "hv",
           conf.int = FALSE, 
           risk.table = TRUE, 
           pval = TRUE,
           palette = "jco", 
           xlab = "Follow up time(m)", 
           legend = c(0.8,0.9), 
           legend.title = "", 
           legend.labs = c("High", "Low"), 签
           xlim=c(0,40),
           break.x.by = 5) 
dev.off()


fit<-survfit(Surv(OS, status2) ~ pred, data = chemo_test)
summary(fit)
pdf("./plot/OS_chemo.pdf",width = 5,height = 5)
ggsurvplot(fit, data = chemo_test,
           surv.median.line = "hv", 
           conf.int = FALSE, 
           risk.table = TRUE,
           pval = TRUE,
           palette = "jco",  
           xlab = "Follow up time(m)", 
           legend = c(0.8,0.9), 
           legend.title = "", 
           legend.labs = c("High", "Low"), 
           xlim=c(0,55),
           break.x.by = 5)  
dev.off()




roc <- pROC::roc(as.numeric(test_set$obs),as.numeric(test_set$R),levels=c(0,1),direction="<")
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

df_chemo <- data.frame(group="chemo",ROC=auc_ci[2], ROC_l=auc_ci[1],ROC_u=auc_ci[3],
                       accuracy=accuracy, acc_l=accuracy_ci[1], acc_u=accuracy_ci[2],
                       sensitivity=sensitivity, sen_l=as.numeric(sensitivity_ci[2]),sen_u=as.numeric(sensitivity_ci[3]),
                       specificity=specificity,spe_l=as.numeric(specificity_ci[2]),spe_u=as.numeric(specificity_ci[3]),
                       PPV=PPV,PPV_l=as.numeric(PPV_ci[2]),PPV_u=as.numeric(PPV_ci[3]),
                       FPV=FPV,FPV_l=as.numeric(FPV_ci[2]),FPV_u=as.numeric(FPV_ci[3]))



##############################################
PFS_tb = chemo_test %>% dplyr::select(Sample,R,PFS,status) %>% inner_join(clinic_in[,c("SUBJID","Sample")],by = "Sample") %>% 
  inner_join(clinic_1[,-2],by = "SUBJID") %>% dplyr::select(-Sample) %>% 
  dplyr::rename(metabo_score=R) %>% 
  dplyr::select(SUBJID,PFS,status,everything())
PFS_tb$SUSMFREQ = factor(PFS_tb$SUSMFREQ,levels = c("从不吸烟","已经戒烟","正在吸烟"),ordered = T) %>% as.numeric()
table(PFS_tb$SUSMFREQ)
PFS_tb$SUDRFREQ = factor(PFS_tb$SUDRFREQ,levels = c("从不饮酒","已经戒酒","正在饮酒"),ordered = T) %>% as.numeric()
table(PFS_tb$SUDRFREQ)

colnames(PFS_tb)
library(coxphf)
multicox <- coxphf(Surv(PFS,status) ~ metabo_score + AGEGROUP + SEX + BLWTGR1 + BLBMGR1 + BLECOG + LMYN + MHTMENUM + SUSMFREQ + SUDRFREQ, data = PFS_tb) 
multisum <- summary(multicox)

item <- colnames(PFS_tb)[4:ncol(PFS_tb)]
HR <- multisum$coefficients %>% exp()
L95CI <- multisum$ci.lower
H95CI <- multisum$ci.upper
pvalue <- multisum$prob
multiresult <- data.frame(item=item,
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue)
multiresult$pvalue[is.na(multiresult$pvalue)]=1

multiresult = multiresult %>% dplyr::mutate(type="multi-cox")

uniresult <- data.frame()  
for(i in colnames(PFS_tb[,4:ncol(PFS_tb)])){   
  unicox <- coxph(Surv(PFS,status) ~ PFS_tb[,i], data = PFS_tb)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  uniresult <- rbind(uniresult,
                     cbind(item=i,
                           HR=unisum$coefficients[,2],
                           L95CI=unisum$conf.int[,3],
                           H95CI=unisum$conf.int[,4],
                           pvalue=unisum$coefficients[,5]
                     ))
}   

uniresult = uniresult %>% dplyr::mutate(type="uni-cox")

plot = multiresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
p1=ggplot()+
  geom_point(data = plot,aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=0.01,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=10,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("PFS_multivariant-Cox")+
  scale_x_log10()+
  theme_bw()

plot = uniresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
plot$HR = as.numeric(plot$HR)
plot$L95CI = as.numeric(plot$L95CI)
plot$H95CI = as.numeric(plot$H95CI)
p2=ggplot(plot)+
  geom_point(data = plot,aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=0.01,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=10,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("PFS_univariant-Cox")+
  scale_x_log10()+
  theme_bw()





##############################################
OS_tb = chemo_test %>% dplyr::select(Sample,R,OS,status2) %>% inner_join(clinic_in[,c("SUBJID","Sample")],by = "Sample") %>% 
  inner_join(clinic_1[,-2],by = "SUBJID") %>% dplyr::select(-Sample) %>% 
  dplyr::rename(metabo_score=R) %>% 
  dplyr::select(SUBJID,OS,status2,everything())
OS_tb$SUSMFREQ = factor(OS_tb$SUSMFREQ,levels = c("从不吸烟","已经戒烟","正在吸烟"),ordered = T) %>% as.numeric()
table(OS_tb$SUSMFREQ)
OS_tb$SUDRFREQ = factor(OS_tb$SUDRFREQ,levels = c("从不饮酒","已经戒酒","正在饮酒"),ordered = T) %>% as.numeric()
table(OS_tb$SUDRFREQ)

multicox <- coxph(Surv(OS,status2) ~ metabo_score + AGEGROUP + SEX + BLWTGR1 + BLBMGR1 + BLECOG + LMYN + MHTMENUM + SUSMFREQ + SUDRFREQ, data = OS_tb) 
multisum <- summary(multicox)

item <- colnames(OS_tb)[4:ncol(OS_tb)]
HR <- multisum$coefficients[,2]
L95CI <- multisum$conf.int[,3]
H95CI <- multisum$conf.int[,4]
pvalue <- multisum$coefficients[,5]
multiresult <- data.frame(item=item,
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue)

multiresult = multiresult %>% dplyr::mutate(type="multi-cox")

uniresult <- data.frame()  
for(i in colnames(OS_tb[,4:ncol(OS_tb)])){   
  unicox <- coxph(Surv(OS,status2) ~ OS_tb[,i], data = OS_tb)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  uniresult <- rbind(uniresult,
                     cbind(item=i,
                           HR=unisum$coefficients[,2],
                           L95CI=unisum$conf.int[,3],
                           H95CI=unisum$conf.int[,4],
                           pvalue=unisum$coefficients[,5]
                     ))
}   

uniresult = uniresult %>% dplyr::mutate(type="uni-cox")

plot = multiresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
p3=ggplot(plot)+
  geom_point(data = plot,aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=0.01,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=10,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("OS_multivariant-Cox")+
  scale_x_log10()+
  theme_bw()

plot = uniresult
plot$pvalue = as.numeric(plot$pvalue)
plot = plot %>% arrange(desc(HR))
plot$item = factor(plot$item,levels = unique(plot$item),ordered = T)
plot$HR = as.numeric(plot$HR)
plot$L95CI = as.numeric(plot$L95CI)
plot$H95CI = as.numeric(plot$H95CI)
p4=ggplot(plot)+
  geom_point(data = plot,aes(x=HR,y=item,size = -log(pvalue),col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=H95CI,y=item,yend=item,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=0.01,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=10,y=item,yend=item,col = log(HR)),
               arrow = arrow(length = unit(0.05, "inches"), type = "closed"))+
  
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,L95CI>0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  
  geom_segment(data = subset(plot,L95CI<0.01 & H95CI<10),aes(x=H95CI,xend=H95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  geom_segment(data = subset(plot,H95CI>10 & L95CI>0.01),aes(x=L95CI,xend=L95CI,y=as.numeric(item)-0.1,yend=as.numeric(item)+0.1,col = log(HR)))+
  scale_color_gradient2(midpoint = 0,low = "navy", mid = "lightgrey", high = "red")+
  geom_vline(xintercept = 1,lty=2)+
  ggtitle("OS_univariant-Cox")+
  scale_x_log10()+
  theme_bw()


pdf("./plot/cox_coef_chemo.pdf",width = 10,height = 8)
cowplot::plot_grid(p1,p2,p3,p4,nrow = 2)
dev.off()


















#############################################################3
prob <- predict(model.tune,surg[,-1],type = "prob")
pre <- predict(model.tune,surg[,-1])
test_set <- data.frame(obs = as.factor(surg$response), NR = prob[,'NR'], R = prob[,'R'], pred=pre)
test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="NR",0,1)})
test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="NR",0,1)})

result_matrix<-confusionMatrix(table(test_set$pred, test_set$obs))
pdf("./plot/cm_sugery.pdf",width = 5,height = 7)
draw_confusion_matrix(result_matrix,"Test Cohort")
dev.off()

pdf("./plot/roc_sugery.pdf",width = 4,height = 3)
rocobj<-roc(test_set$obs, test_set$R)
auc_plot(rocobj)
dev.off()


surg_test = cbind(rownames(surg),test_set)
colnames(surg_test)[1]="Sample"

surg_test = surg_test %>% dplyr::mutate(t_or_f = ifelse(obs==pred,"T","F"))
table(surg_test$t_or_f)

# roc <- pROC::roc(as.numeric(test_set$obs),as.numeric(pre),levels=c(0,1),direction="<")
# # roc_result <- coords(roc, "best")
# auc <- roc %>% auc()
# auc_ci <- ci(roc)

roc <- pROC::roc(as.numeric(test_set$obs),as.numeric(test_set$R),levels=c(0,1),direction="<")
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

df_surg <- data.frame(group="surg",ROC=auc_ci[2], ROC_l=auc_ci[1],ROC_u=auc_ci[3],
                      accuracy=accuracy, acc_l=accuracy_ci[1], acc_u=accuracy_ci[2],
                      sensitivity=sensitivity, sen_l=as.numeric(sensitivity_ci[2]),sen_u=as.numeric(sensitivity_ci[3]),
                      specificity=specificity,spe_l=as.numeric(specificity_ci[2]),spe_u=as.numeric(specificity_ci[3]),
                      PPV=PPV,PPV_l=as.numeric(PPV_ci[2]),PPV_u=as.numeric(PPV_ci[3]),
                      FPV=FPV,FPV_l=as.numeric(FPV_ci[2]),FPV_u=as.numeric(FPV_ci[3]))









#############################################################3
prob <- predict(model.tune,rad[,-1],type = "prob")
pre <- predict(model.tune,rad[,-1])
test_set <- data.frame(obs = as.factor(rad$response), NR = prob[,'NR'], R = prob[,'R'], pred=pre)
test_set$obs <- sapply(1:nrow(test_set),function(u){ifelse(test_set$obs[u]=="NR",0,1)})
test_set$pred <- sapply(1:nrow(test_set),function(u){ifelse(test_set$pred[u]=="NR",0,1)})

result_matrix<-confusionMatrix(table(test_set$pred, test_set$obs))
pdf("./plot/cm_rad.pdf",width = 5,height = 7)
draw_confusion_matrix(result_matrix,"Test Cohort")
dev.off()

pdf("./plot/roc_rad.pdf",width = 4,height = 3)
rocobj<-roc(test_set$obs, test_set$R)
auc_plot(rocobj)
dev.off()

rad_test = cbind(rownames(rad),test_set)
colnames(rad_test)[1]="Sample"

rad_test = rad_test %>% dplyr::mutate(t_or_f = ifelse(obs==pred,"T","F"))
table(rad_test$t_or_f)

# roc <- pROC::roc(as.numeric(test_set$obs),as.numeric(pre),levels=c(0,1),direction="<")
# # roc_result <- coords(roc, "best")
# auc <- roc %>% auc()
# auc_ci <- ci(roc)

roc <- pROC::roc(as.numeric(test_set$obs),as.numeric(test_set$R),levels=c(0,1),direction="<")
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

df_rad <- data.frame(group="rad",ROC=auc_ci[2], ROC_l=auc_ci[1],ROC_u=auc_ci[3],
                     accuracy=accuracy, acc_l=accuracy_ci[1], acc_u=accuracy_ci[2],
                     sensitivity=sensitivity, sen_l=as.numeric(sensitivity_ci[2]),sen_u=as.numeric(sensitivity_ci[3]),
                     specificity=specificity,spe_l=as.numeric(specificity_ci[2]),spe_u=as.numeric(specificity_ci[3]),
                     PPV=PPV,PPV_l=as.numeric(PPV_ci[2]),PPV_u=as.numeric(PPV_ci[3]),
                     FPV=FPV,FPV_l=as.numeric(FPV_ci[2]),FPV_u=as.numeric(FPV_ci[3]))







###############################################
df = rbind(df_train,df_vali,df_chemo,df_surg,df_rad)
df$group = factor(df$group,levels = df$group,ordered = T)

save(df,file = "./plot/all_predict_res.rdata")

library(Seurat)
pdf("./plot/accuracy.pdf",width = 4,height = 4)
ggplot(df)+
  geom_segment(aes(x=as.numeric(group)-0.1,xend=as.numeric(group)+0.1,y=acc_l,yend=acc_l))+
  geom_segment(aes(x=as.numeric(group)-0.1,xend=as.numeric(group)+0.1,y=acc_u,yend=acc_u))+
  geom_segment(aes(x=as.numeric(group),xend=as.numeric(group),y=acc_l,yend=acc_u))+
  geom_point(aes(x=as.numeric(group),y=accuracy))+
  geom_text(aes(x = as.numeric(group), y = acc_u + 0.05, label = sprintf("%.3f", accuracy))) +
  scale_x_continuous(
    breaks = as.numeric(unique(df$group)),
    labels = unique(df$group)
  ) +
  ylim(0.4,1.1)+
  theme_bw()+
  RotatedAxis()
dev.off()





df_1 = df %>% dplyr::mutate(accuracy_1=str_c(round(accuracy,3),"\n(",round(acc_l,3),"-",round(acc_u,3),")")) %>% 
  dplyr::mutate(sensitivity_1 = str_c(round(sensitivity,3),"\n(",round(sen_l,3),"-",round(sen_u,3),")")) %>% 
  dplyr::mutate(specificity_1 = str_c(round(specificity,3),"\n(",round(spe_l,3),"-",round(spe_u,3),")")) %>% 
  dplyr::mutate(PPV_1 = str_c(round(PPV,3),"\n(",round(PPV_l,3),"-",round(PPV_u,3),")")) %>% 
  dplyr::mutate(FPV_1 = str_c(round(FPV,3),"\n(",round(FPV_l,3),"-",round(FPV_u,3),")")) %>% 
  dplyr::select(group,accuracy_1,sensitivity_1,specificity_1,PPV_1,FPV_1) 
colnames(df_1)[2:6] = str_remove(colnames(df_1)[2:6],"_1") 
colnames(df_1)[2:4] = colnames(df_1)[2:4] %>% str_to_title()
table(training$response)
table(validation$response)
table(chemo$response)
table(surg$response)
table(rad$response)
df_1$group = c("Training (43 vs 126)","Validation (18 vs 54)","Test cohort 1 (23 vs 67)",
               "Test cohort 2 (29 vs 79)","Test cohort 3 (16 vs 38)")
# df_1$NR = c(43,18,24,12,33)
# df_1$R = c(126,54,67,29,88)
# df_1$list = sapply(1:nrow(df_1),function(kk){
#   return(list(c(NR=df_1$NR[kk],R=df_1$R[kk])))
# })

library(gt)
library(tidyverse)
library(glue)
library(gtExtras)

# table(surg_test$obs)

gt_bl = df_1 %>% 
  gt(rowname_col = "group") %>%
  tab_header(
    title = md("**Prediction Result of Different Cohort**")
  ) %>%
  tab_style(
    style = list(
      cell_text(font = "Oswald",
                align = "center",
                weight = "normal")
    ),
    locations = list(
      cells_body(everything())
    )
  )%>%
  tab_style(
    style = list(
      cell_text(
        font = "Oswald",
        align = "center",  
        weight = "bold"   
      )
    ),
    locations = list(
      cells_column_labels(everything())  
    )
  ) %>%
  tab_style(
    style = list(
      cell_text(
        font = "Oswald",
        align = "center",  
        weight = "bold"    
      )
    ),
    locations = list(
      cells_stub()  
    )
  ) 

library(webshot2)
gt_bl %>% gtsave("tab_1.html", path = "./plot/")
