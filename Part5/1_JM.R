#!/usr/bin/Rscript
rm(list = ls())
setwd("~/help_for_others/DYQ/output/")
library(stringr)
library(gridExtra)
library(future)
# library(sva)
# library(e1071)
library(pROC)
library(ROCit)
library(caret)
library(doParallel)
library(cancerclass)
library(tidyverse)
load("./1_preclean/data/data_preclean.rdata")
setwd("./7_JM_1/")
dir.create("csv")
setwd("./csv/")

clinic_in$CNSR1 <- sapply(1:nrow(clinic_in),function(uu){
  if(clinic_in$CNSR1[uu]==0){
    return(1)
  }else{
    return(0)
  }
})


save(data_in,clinic_in,file = "data_preclean_rm_sample.rdata")
sample_all = clinic_in$Sample
data_in_for_csv = data_in %>% dplyr::select(Compounds,any_of(sample_all)) %>% 
  column_to_rownames("Compounds") %>% t() %>% data.frame(.,check.names = FALSE) %>% 
  rownames_to_column("Sample") %>% inner_join(clinic_in[,c("Sample","cycle")],by = "Sample") %>% 
  dplyr::select(Sample,cycle,everything())

write.csv(data_in_for_csv,file = "~/help_for_others/DYQ/output/7_JM_1/csv/all_df.csv",row.names = FALSE)



##########################################
setwd("~/help_for_others/DYQ/output/7_JM_1/")
dir.create("pre_metabo")
setwd("./pre_metabo/")
# write.csv(N_df,"all_response_df.csv",row.names = FALSE)

if(T){
  library(MetaboAnalystR)
  mSet<-InitDataObjects("conc", "stat", FALSE)
  mSet<-Read.TextData(mSet, "../csv/all_df.csv", "rowu", "disc");
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet)  ####替换0值或缺失值
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "NULL", "LogNorm", "AutoNorm", "S10T0", ratio=FALSE, ratioNum=20)
  # mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###按代谢物作图
  # mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###按样本作图
  norm=mSet$dataSet$norm %>% rownames_to_column("Sample")  
  save(norm,file = "norm.rdata")
}


load("norm.rdata")




library(survminer) 
library(survival) 





load("../../../data/blood_time.rdata")
blood_time = blood_time %>% dplyr::mutate(cycle=ifelse(Cycle=="筛选期(D-7~D-1)","C1D1",Cycle))

clinic_in_merge = clinic_in %>% inner_join(blood_time[,c("SUBJID","cycle","Time")],by = c("SUBJID","cycle")) %>% 
  dplyr::mutate(measure_time=ifelse(Time<0,0,Time)) %>% dplyr::mutate(measure_time=measure_time/30.44)
clinic_in = clinic_in_merge

data <- norm %>% inner_join(clinic_in[,c("Sample","SUBJID","measure_time","AVAL1","CNSR1")],.,by = "Sample") 


sd_tb = apply(data[,-c(1:5)],2,var) %>% data.frame()
# hist(sd_tb$.)
# quantile(sd_tb$.,0.2)
hvg = sd_tb %>% dplyr::filter(`.`>quantile(sd_tb$`.`,0.2)) %>% rownames(.)
data = data %>% dplyr::select(1:5,any_of(hvg))


#######################################################
if(T){
  library(survival)
  library(survminer)
  #########dor
  colnames(data)[1:30]
  # data_in = data_2_1_dor
  # table(sample_3_1$cycle)
  data_in = data
  kk=6
  
  res <- lapply(6:ncol(data_in),function(kk){
    print(kk)
    a = data_in %>% dplyr::select(1:5,kk) %>% dplyr::rename(metabo=colnames(.)[6])
    a$measure_time <- as.numeric(a$measure_time)
    a$start_time <- a$measure_time
    tt=12
    
    a1 <- lapply(1:length(unique(a$SUBJID)),function(tt){
      kk=a %>% dplyr::filter(SUBJID==unique(a$SUBJID)[tt]) %>% arrange(measure_time)
      kk$stop_time <- sapply(1:nrow(kk),function(uu){
        if(uu<nrow(kk)){
          return(min(kk$start_time[uu+1],kk$AVAL1[uu]))
        }else{
          return(max(kk$AVAL1[uu],kk$start_time[uu]))
        }
      })
      kk$start_time <- sapply(1:nrow(kk),function(uu){
        if(uu<nrow(kk)){
          return(kk$start_time[uu])
        }else{
          return(min(kk$AVAL1[uu],kk$start_time[uu]))
        }
      })
      kk$event <- sapply(1:nrow(kk),function(uu){
        if(kk$stop_time[uu]<kk$AVAL1[uu]){
          return(0)
        }else{
          return(kk$CNSR1[uu])
        }
      })
      return(kk)
    }) %>% do.call(rbind,.) %>% distinct(.,.keep_all = T)
    
    table(a1$start_time<a1$stop_time)
    a1 = a1[!(a1$start_time>a1$stop_time),]
    td.Cox <- coxph(Surv(start_time, stop_time, event) ~ metabo,
                    data = a1)
    td=summary(td.Cox)
    return(data.frame(metabo=colnames(data_in)[kk],coef=td$conf.int["metabo","exp(coef)"],
                      l95=td$conf.int["metabo","lower .95"],u95=td$conf.int["metabo","upper .95"],
                      p=td$coefficients[,"Pr(>|z|)"][1]))
    
  }) %>% do.call(rbind,.)
  
  # sig_mt = res %>% dplyr::filter(p<0.05)
  res$padj = p.adjust(res$p,method="fdr")
  save(res,data_in,file = "../joint_res_1/joint_tb_and_res_1.rdata")
}




#########################################################################
colnames(data)[1:30]
data_in = data
kk=6

library(survival)
library(survminer)
# library(parallel)
library(JMbayes2)

library(JMbayes2)
library(parallel)

kk=1104
library(parallel)
# mc<-getOption("mc.cores",10)
res <- lapply(6:ncol(data_in),function(kk){
  print(kk)
  a = data_in %>% dplyr::select(1:5,any_of(kk)) %>% dplyr::rename(metabo=colnames(.)[6])
  library(nlme)
  fit<-try(
    lmefit<-lme(metabo ~ measure_time,
                random= ~ measure_time|SUBJID,data=a,
                control = lmeControl(niterEM = 10000, optimMethod = "L-BFGS-B", opt = "nlminb",step.max = 0.5)),
    silent = TRUE)
  if ('try-error' %in% class(fit)){
    
    return(data.frame(metabo=colnames(data_in)[kk],coef=NA,
                      se=NA,
                      p=NA))
  }else{
    
    lmefit
    s=summary(lmefit)
    
    a.id <- a[!duplicated(a$SUBJID), ]
    coxfit <- coxph(Surv(AVAL1, CNSR1) ~ 1,
                    data = a.id, x = TRUE)
    
    library(JM)
    fit<-try(
      jmod<-jointModel(lmefit,coxfit,timeVar="measure_time",method = "piecewise-PH-aGH"),
      silent = TRUE)
    if ('try-error' %in% class(fit)){
      return(data.frame(metabo=colnames(data_in)[kk],coef=NA,
                        se=NA,
                        p=NA))
    }else{
      
      js=summary(jmod)
      return(data.frame(metabo=colnames(data_in)[kk],coef=js$`CoefTable-Event`["Assoct","Value"],
                        se=js$`CoefTable-Event`["Assoct","Std.Err"],
                        p=js$`CoefTable-Event`["Assoct","p-value"]))
    }
  }
}) %>% do.call(rbind,.) %>% dplyr::mutate(l95=coef-1.96*se,u95=coef+1.96*se)

res$padj = p.adjust(res$p,method="fdr")

save(res,data_in,file = "../joint_res_2/joint_tb_and_res_2.rdata")
