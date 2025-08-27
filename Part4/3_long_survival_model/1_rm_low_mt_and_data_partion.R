rm(list = ls())
library(tidyverse)
setwd("~/help_for_others/DYQ/output/6_long_survival_model/")  

mSet = read.csv("~/help_for_others/DYQ/output/5_delta_metabo_response_14m_discovering/csv/C1C6_pair_df.csv",check.names = FALSE)
mSet = mSet %>% dplyr::mutate(ID=str_remove(Sample,"C[1|6]D1")) %>% dplyr::select(ID,everything())
mSet = mSet %>% 
  dplyr::select(-ID)

write.csv(mSet,"./csv/rm_sample_C1C6_df.csv",row.names = FALSE)


rm(list = ls())
dir.create("./C1C6_delta_metabo")
setwd("./C1C6_delta_metabo/")
# .libPaths()
library(MetaboAnalystR)
library(tidyverse)
library(caret)

mSet<-InitDataObjects("conc", "stat", FALSE)
mSet<-Read.TextData(mSet, "~/help_for_others/DYQ/output/6_long_survival_model//csv/rm_sample_C1C6_df.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet)  
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE,ratioNum = 20)

mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)  
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   

norm = mSet$dataSet$norm
info = mSet$dataSet$cls

range(norm)


norm_1 = norm %>% rownames_to_column("ID") %>% dplyr::mutate(cycle=str_extract(ID,"C[1|6]D1")) %>% 
  dplyr::mutate(patient=str_remove(ID,"C[1|6]D1")) %>% dplyr::mutate(is_long_PFS=info) %>% 
  dplyr::select(ID,patient,cycle,is_long_PFS,everything())

norm_2 = norm_1 %>% pivot_longer(cols = !c(1:4),names_to = "mt",values_to="value")

norm_2$is_long_PFS = as.character(norm_2$is_long_PFS)
norm_3_C1 = norm_2 %>% dplyr::select(-ID) %>% dplyr::filter(cycle=="C1D1") %>% arrange(patient,mt) %>% dplyr::rename(C1D1=value)
norm_3_C6 = norm_2 %>% dplyr::select(-ID) %>% dplyr::filter(cycle=="C6D1") %>% arrange(patient,mt) %>% dplyr::rename(C6D1=value)

norm_3 = cbind(norm_3_C1[,c(1,3,4,5)],norm_3_C6[,5])
norm_3$delta = norm_3$C6D1-norm_3$C1D1







################
base_expr = norm_3 %>% group_by(mt) %>% dplyr::summarize(mean_C1=mean(C1D1),mean_C6=mean(C6D1)) %>% 
  dplyr::mutate(mean_all=(mean_C1+mean_C6)/2)
hist(base_expr$mean_all,breaks = 1000)

preserve_expr = base_expr %>% dplyr::filter(mean_all > quantile(base_expr$mean_all,0.05))


norm_3_process = norm_3 %>% dplyr::select(1,2,3,6) %>% pivot_wider(id_cols = c(1,2),names_from = "mt",values_from = "delta")
colnames(norm_3_process)[1]="Sample"
norm_3_process$Sample = sapply(1:nrow(norm_3_process),function(kk){
  print(kk)
  if(nchar(norm_3_process$Sample[kk])==4){
    return(str_c(0,norm_3_process$Sample[kk],sep = ""))
  }else{
    return(norm_3_process$Sample[kk])
  }
})


norm_3_process = norm_3_process %>% dplyr::select(1:2,any_of(preserve_expr$mt))







#################################################

  set.seed(10086)
  trainIndex <- createDataPartition(norm_3_process$is_long_PFS, p = .7,  ## 80% training set; 20% validation set
                                    list = FALSE, 
                                    times = 1)
  training <- norm_3_process[trainIndex,]
  validation  <- norm_3_process[-trainIndex,]
  
  
  setwd("~/help_for_others/DYQ/output/6_long_survival_model/csv/")
  save(training,validation,trainIndex,file = "./training_validation_df.rdata")
  write.csv(training,file = "training_delta_df.csv",row.names = FALSE)
  write.csv(validation,file = "validation_delta_df.csv",row.names = FALSE)
  write.csv(norm_3_process,file = "all_delta_df.csv",row.names = FALSE)

  
  
