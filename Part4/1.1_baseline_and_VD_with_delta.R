rm(list = ls())
load("~/help_for_others/DYQ/output/1_preclean/data/data_preclean.rdata")
table(clinic_in$cycle)



clinic_1 = clinic_in %>% dplyr::filter(cycle %in% c("C1D1","C6D1"))
data_1 = data_in %>% dplyr::select(1:4,any_of(clinic_1$Sample))

setwd("~/help_for_others/DYQ/output/")
dir.create("8_compare_with_mt_change_in_R_and_NR")


setwd("./8_compare_with_mt_change_in_R_and_NR/csv/")

data_for_csv = data_in %>% dplyr::select(-c(2:4)) %>% column_to_rownames("Compounds") %>% t() %>% 
  data.frame(.,check.names = FALSE) %>% rownames_to_column("Sample") %>% 
  inner_join(clinic_1[,c("Sample","cycle")],by = "Sample") %>% 
  dplyr::select(Sample,cycle,everything())

write.csv(data_for_csv,file = "./data_all_C1C6.csv",row.names = FALSE)



setwd("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/")
dir.create("metabo_process")
setwd("./metabo_process/")
dir.create("./all_C1C6")
setwd("./all_C1C6/")


# .libPaths()
library(MetaboAnalystR)
library(tidyverse)

mSet<-InitDataObjects("conc", "stat", FALSE)
mSet<-Read.TextData(mSet, "../../csv/data_all_C1C6.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet)  ####替换0值或缺失值
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE,ratioNum = 20)

mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###按代谢物作图
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###按样本作图

norm = mSet$dataSet$norm
info = mSet$dataSet$meta.info %>% unlist()

norm_1 = norm %>% rownames_to_column("Sample") %>% 
  inner_join(clinic_1[,c("Sample","SUBJID","cycle")],by = "Sample") %>% 
  dplyr::select(Sample,SUBJID,cycle,everything())
norm_2 = norm_1 %>% pivot_longer(cols = !c(1:3),names_to = "mt",values_to="value")

norm_3_C1 = norm_2 %>% dplyr::select(-Sample) %>% dplyr::filter(cycle=="C1D1") %>% 
  arrange(SUBJID,mt) %>% dplyr::rename(C1D1=value)
norm_3_C6 = norm_2 %>% dplyr::select(-Sample) %>% dplyr::filter(cycle=="C6D1") %>% 
  arrange(SUBJID,mt) %>% dplyr::rename(C6D1=value)
  
norm_3 = norm_3_C1[,c(1,3,4)] %>% inner_join(norm_3_C6[,c(1,3,4)],by = c("SUBJID","mt"))
norm_3$delta = norm_3$C6D1-norm_3$C1D1
norm_3_process = norm_3 %>% dplyr::select(1,2,5) %>% pivot_wider(id_cols = c(1),names_from = "mt",values_from = "delta")

  
  
clinic_need = clinic_1[,c(2,6,7,9,11,12,13,14,18,19,20,21,25,26,35)] %>% distinct(.,.keep_all = T) %>% 
  dplyr::filter(SUBJID %in% norm_3_process$SUBJID)


vd_res_sum = lapply(5:20,function(t){
  print(t)
  mt_expr = norm_3_process %>% column_to_rownames("SUBJID")

  clinic_info_1 = clinic_need %>% dplyr::select(1:12) %>% column_to_rownames("SUBJID")
  clinic_info_2 = clinic_need %>% dplyr::select(1,13:14) %>% column_to_rownames("SUBJID") %>% 
    dplyr::mutate(CNSR1=ifelse(CNSR1==1,0,1)) %>% 
    dplyr::mutate(PFS_larger=ifelse(AVAL1>t,"yes",ifelse(CNSR1==1,"no","unknown"))) %>% 
    dplyr::select(PFS_larger)
  clinic_info_3 = clinic_need %>% dplyr::select(1,15) %>% column_to_rownames("SUBJID")

  source("~/help_for_others/DYQ/code/8_compare_with_mt_change_in_R_and_NR/VD.R")
  
  vd.vars <- c(colnames(clinic_info_1),colnames(clinic_info_3),colnames(clinic_info_2))
  meta.data <- cbind(clinic_info_1,clinic_info_3,clinic_info_2)
  table(meta.data$PFS_larger,meta.data$R)
  
  unknown = meta.data %>% dplyr::filter(PFS_larger=="unknown") %>% rownames(.)
  
  ras.data <- mt_expr %>% rownames_to_column("ID") %>% dplyr::filter(!ID %in% unknown) %>% column_to_rownames("ID")
  meta.data = meta.data %>% rownames_to_column("ID") %>% dplyr::filter(!ID %in% unknown) %>% column_to_rownames("ID")
  
  vd.res_def <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = 5)
  
  return(vd.res_def)
})

save(vd_res_sum,file = "../../VD/vd_res_sum_delta.rdata")


