rm(list = ls())
setwd("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/")
library(tidyverse)
dir.create("metabo_process")
setwd("./metabo_process/")
dir.create("./all_C1C6")
setwd("./all_C1C6/")

load("~/help_for_others/DYQ/output/1_preclean/data/data_preclean.rdata")
table(clinic_in$cycle)



clinic_1 = clinic_in %>% dplyr::filter(cycle %in% c("C1D1","C6D1"))


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
# norm_2 = norm_1 %>% pivot_longer(cols = !c(1:3),names_to = "mt",values_to="value")

norm_2_C1 = norm_1%>% dplyr::filter(cycle=="C1D1")  %>% dplyr::select(-c(Sample,cycle)) 
norm_2_C6 = norm_1 %>% dplyr::filter(cycle=="C6D1")  %>% dplyr::select(-c(Sample,cycle)) 


clinic_C1 = clinic_1[,c(2,6,7,9,11,12,13,14,18,19,20,21,25,26,35)] %>% distinct(.,.keep_all = T) %>% 
  dplyr::filter(SUBJID %in% norm_2_C1$SUBJID)

clinic_C6 = clinic_1[,c(2,6,7,9,11,12,13,14,18,19,20,21,25,26,35)] %>% distinct(.,.keep_all = T) %>% 
  dplyr::filter(SUBJID %in% norm_2_C6$SUBJID)
# library(vegan)
# data(mite)#物种数据
# data(mite.env) #环境数据1
# data(mite.pcnm) #环境数据2

vd_res_sum = lapply(4:20,function(t){
  print(t)
  mt_expr = norm_2_C1 %>% column_to_rownames("SUBJID")
  # table(is.na(mt_expr))
  # table(is.infinite(unlist(mt_expr)))
  # mt_expr = exp(mt_expr)
  
  clinic_info_1 = clinic_C1 %>% dplyr::select(1:12) %>% column_to_rownames("SUBJID")
  clinic_info_2 = clinic_C1 %>% dplyr::select(1,13:14) %>% column_to_rownames("SUBJID") %>% 
    dplyr::mutate(CNSR1=ifelse(CNSR1==1,0,1)) %>% 
    dplyr::mutate(PFS_larger=ifelse(AVAL1>t,"yes",ifelse(CNSR1==1,"no","unknown"))) %>% 
    dplyr::select(PFS_larger)
  clinic_info_3 = clinic_C1 %>% dplyr::select(1,15) %>% column_to_rownames("SUBJID")
  # table(is.na(clinic_info_1))
  # table(is.na(clinic_info_3))
  
  source("~/help_for_others/DYQ/code/8_compare_with_mt_change_in_R_and_NR/VD.R")
  
  vd.vars <- c(colnames(clinic_info_1),colnames(clinic_info_3),colnames(clinic_info_2))
  meta.data <- cbind(clinic_info_1,clinic_info_3,clinic_info_2)
  table(meta.data$PFS_larger,meta.data$R)
  
  unknown = meta.data %>% dplyr::filter(PFS_larger=="unknown") %>% rownames(.)
  
  ras.data <- mt_expr %>% rownames_to_column("ID") %>% dplyr::filter(!ID %in% unknown) %>% column_to_rownames("ID")
  meta.data = meta.data %>% rownames_to_column("ID") %>% dplyr::filter(!ID %in% unknown) %>% column_to_rownames("ID")
  
  vd.res_def <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = 5)
  
  # vpa1<- varpart(exp(mt_expr), clinic_info_1, clinic_info_2,transfo="hel")
  # vpa1
  # plot(vpa1, bg=2:3)
  
  return(vd.res_def)
})

save(vd_res_sum,file = "../../VD/vd_res_sum_C1.rdata")







vd_res_sum = lapply(4:20,function(t){
  print(t)
  mt_expr = norm_2_C6 %>% column_to_rownames("SUBJID")
  # table(is.na(mt_expr))
  # table(is.infinite(unlist(mt_expr)))
  # mt_expr = exp(mt_expr)
  
  clinic_info_1 = clinic_C6 %>% dplyr::select(1:12) %>% column_to_rownames("SUBJID")
  clinic_info_2 = clinic_C6 %>% dplyr::select(1,13:14) %>% column_to_rownames("SUBJID") %>% 
    dplyr::mutate(CNSR1=ifelse(CNSR1==1,0,1)) %>% 
    dplyr::mutate(PFS_larger=ifelse(AVAL1>t,"yes",ifelse(CNSR1==1,"no","unknown"))) %>% 
    dplyr::select(PFS_larger)
  clinic_info_3 = clinic_C6 %>% dplyr::select(1,15) %>% column_to_rownames("SUBJID")
  # table(is.na(clinic_info_1))
  # table(is.na(clinic_info_3))
  
  source("~/help_for_others/DYQ/code/8_compare_with_mt_change_in_R_and_NR/VD.R")
  
  vd.vars <- c(colnames(clinic_info_1),colnames(clinic_info_3),colnames(clinic_info_2))
  meta.data <- cbind(clinic_info_1,clinic_info_3,clinic_info_2)
  table(meta.data$PFS_larger,meta.data$R)
  
  unknown = meta.data %>% dplyr::filter(PFS_larger=="unknown") %>% rownames(.)
  
  ras.data <- mt_expr %>% rownames_to_column("ID") %>% dplyr::filter(!ID %in% unknown) %>% column_to_rownames("ID")
  meta.data = meta.data %>% rownames_to_column("ID") %>% dplyr::filter(!ID %in% unknown) %>% column_to_rownames("ID")
  
  vd.res_def <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = 5)
  
  # vpa1<- varpart(exp(mt_expr), clinic_info_1, clinic_info_2,transfo="hel")
  # vpa1
  # plot(vpa1, bg=2:3)
  
  return(vd.res_def)
})

save(vd_res_sum,file = "../../VD/vd_res_sum_C6.rdata")







load("../../VD/vd_res_sum_C6.rdata")

PFS_longer = lapply(1:length(vd_res_sum),function(kk){
  tmp = vd_res_sum[[kk]] %>% dplyr::mutate(PFS_cf=kk+3) %>% 
    dplyr::select(PFS_larger,PFS_cf) %>% rownames_to_column("mt")
}) %>% do.call(rbind,.)
PFS_longer$PFS_cf = factor(PFS_longer$PFS_cf,levels = 4:20,ordered = T)
mt_sel = PFS_longer %>% dplyr::filter(PFS_larger>0.05) %>% pull(mt) %>% unique()
PFS_longer = PFS_longer %>% dplyr::filter(mt %in% mt_sel) %>% 
  dplyr::filter(PFS_cf>=5) %>% 
  dplyr::mutate(mt_time="C6")

PFS_longer_C6 = PFS_longer


load("../../VD/vd_res_sum_C1.rdata")

PFS_longer = lapply(1:length(vd_res_sum),function(kk){
  tmp = vd_res_sum[[kk]] %>% dplyr::mutate(PFS_cf=kk+3) %>% 
    dplyr::select(PFS_larger,PFS_cf) %>% rownames_to_column("mt")
}) %>% do.call(rbind,.)
PFS_longer$PFS_cf = factor(PFS_longer$PFS_cf,levels = 4:20,ordered = T)
mt_sel = PFS_longer %>% dplyr::filter(PFS_larger>0.05) %>% pull(mt) %>% unique()
PFS_longer = PFS_longer %>% dplyr::filter(mt %in% mt_sel) %>% 
  dplyr::filter(PFS_cf>=5) %>% 
  dplyr::mutate(mt_time="C1")

PFS_longer_C1 = PFS_longer

load("../../VD/vd_res_sum_delta.rdata")

PFS_longer = lapply(1:length(vd_res_sum),function(kk){
  tmp = vd_res_sum[[kk]] %>% dplyr::mutate(PFS_cf=kk+3) %>% 
    dplyr::select(PFS_larger,PFS_cf) %>% rownames_to_column("mt")
}) %>% do.call(rbind,.)
PFS_longer$PFS_cf = factor(PFS_longer$PFS_cf,levels = 4:20,ordered = T)
mt_sel = PFS_longer %>% dplyr::filter(PFS_larger>0.05) %>% pull(mt) %>% unique()
PFS_longer = PFS_longer %>% dplyr::filter(mt %in% mt_sel) %>% 
  dplyr::filter(PFS_cf>=5) %>% 
  dplyr::mutate(mt_time="delta")

PFS_longer_delta = PFS_longer

PFS_longer = rbind(PFS_longer_C1,PFS_longer_C6,PFS_longer_delta)


library(ggridges)
library(wesanderson)
pdf("../../../8_compare_with_mt_change_in_R_and_NR/plot/VD_with_only_C1_or_C6_or_delta.pdf",width = 8,height = 4)
ggplot(PFS_longer, aes(x=PFS_cf, y=PFS_larger))+
  geom_line(aes(group=mt,col=mt_time),lwd=0.1)+
  scale_color_manual(values = wes_palette("GrandBudapest1",3))+
  # geom_smooth()+
  facet_grid(~mt_time)+
  theme_bw()+
  theme(strip.background = element_blank())
dev.off()
