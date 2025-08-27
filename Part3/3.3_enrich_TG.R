rm(list = ls())
setwd("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/")
library(tidyverse)
dir.create("metabo_process")
setwd("./metabo_process/")
dir.create("./all_C1C6")
setwd("./all_C1C6/")

load("~/help_for_others/DYQ/output/1_preclean/data/data_preclean.rdata")
table(clinic_in$cycle)

load("../../../1_preclean/data/data_C1_preclean.rdata")
mt_info = data_in_C1 %>% dplyr::select(Compounds,`Class I`,`Class II`) %>% 
  dplyr::rename(mt=Compounds)
class_ID_1 = mt_info %>% dplyr::select(mt,`Class I`) %>% dplyr::rename(class=colnames(.)[2])
class_ID_2 = mt_info %>% dplyr::select(mt,`Class II`) %>% dplyr::rename(class=colnames(.)[2])



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

mt_expr = norm_2_C1 %>% column_to_rownames("SUBJID")

load("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/VD/mt_enricher_TG.rdata")

mt_ls_need = TG_unsature_ls

plot_res_unsature = lapply(1:length(mt_ls_need),function(kk){
  
  norm_sel = mt_expr %>% dplyr::select(any_of(mt_ls_need[[kk]])) %>% rownames_to_column("SUBJID")
  
  library(RobustRankAggreg)
  library(survival)
  library(survminer)
  method_choose="mean"
  glist<-list()
  i=2
  for (i in 2:ncol(norm_sel)){
    data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
    glist[[i-1]]<-data$SUBJID
    names(glist)[i-1]<-colnames(norm_sel)[i]
  }
  library(RobustRankAggreg)
  r<-rankMatrix(glist,full = TRUE)
  RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
  
  RRA = RRA %>% dplyr::rename(SUBJID=Name) %>% inner_join(clinic_C1[,c("SUBJID","R","AVAL1","CNSR1")],by = "SUBJID") %>% 
    dplyr::mutate(class=names(mt_ls_need)[kk])
  

  model =glm(as.factor(R)~Score,data = RRA,family = binomial)
  s=summary(model)
  confint(model)
  return(data.frame(class=names(mt_ls_need)[kk],beita=s$coefficients["Score",1],l_b=confint(model)[2,1],u_b=confint(model)[2,2]))
  # return(RRA)
  
}) %>% do.call(rbind,.)


mt_ls_need = TG_sature_ls

plot_res_sature = lapply(1:length(mt_ls_need),function(kk){
  
  norm_sel = mt_expr %>% dplyr::select(any_of(mt_ls_need[[kk]])) %>% rownames_to_column("SUBJID")
  
  library(RobustRankAggreg)
  library(survival)
  library(survminer)
  method_choose="mean"
  glist<-list()
  i=2
  for (i in 2:ncol(norm_sel)){
    data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
    glist[[i-1]]<-data$SUBJID
    names(glist)[i-1]<-colnames(norm_sel)[i]
  }
  library(RobustRankAggreg)
  r<-rankMatrix(glist,full = TRUE)
  RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
  
  RRA = RRA %>% dplyr::rename(SUBJID=Name) %>% inner_join(clinic_C1[,c("SUBJID","R","AVAL1","CNSR1")],by = "SUBJID") %>% 
    dplyr::mutate(class=names(mt_ls_need)[kk])
  

  model =glm(as.factor(R)~Score,data = RRA,family = binomial)
  s=summary(model)
  confint(model)
  return(data.frame(class=names(mt_ls_need)[kk],beita=s$coefficients["Score",1],l_b=confint(model)[2,1],u_b=confint(model)[2,2]))
  # return(RRA)
  
}) %>% do.call(rbind,.)



plot_res = rbind(plot_res_sature,plot_res_unsature) %>% 
  dplyr::mutate(cutoff=str_extract(class,"[0-9]+$")) %>% 
  dplyr::mutate(cutoff=as.numeric(cutoff)) %>% 
  arrange(cutoff) %>% 
  dplyr::mutate(cutoff=factor(cutoff,levels = unique(cutoff),ordered = T)) %>% 
  dplyr::mutate(type=str_remove(class,"[0-9]+$"))

plot_res_sature = plot_res %>% dplyr::filter(type=="TG_sature")
plot_res_unsature = plot_res %>% dplyr::filter(type=="TG_unsature")


pdf("../../plot/TG_cutoff_logistic_unsature_C1.pdf",width = 7,height = 4)
ggplot()+
  geom_segment(data = plot_res_sature,aes(x=as.numeric(cutoff)-0.15,xend=as.numeric(cutoff)-0.15,y=l_b,yend=u_b,col=type))+
  geom_segment(data = plot_res_sature,aes(x=as.numeric(cutoff)-0.25,xend=as.numeric(cutoff)-0.05,y=l_b,yend=l_b,col=type))+
  geom_segment(data = plot_res_sature,aes(x=as.numeric(cutoff)-0.25,xend=as.numeric(cutoff)-0.05,y=u_b,yend=u_b,col=type))+
  geom_point(data = plot_res_sature,aes(x=as.numeric(cutoff)-0.15,y=beita,col=type))+
  geom_line(data = plot_res_sature,aes(x=as.numeric(cutoff)-0.15,y=beita,col=type))+
  geom_segment(data = plot_res_unsature,aes(x=as.numeric(cutoff)+0.15,xend=as.numeric(cutoff)+0.15,y=l_b,yend=u_b,col=type))+
  geom_segment(data = plot_res_unsature,aes(x=as.numeric(cutoff)+0.05,xend=as.numeric(cutoff)+0.25,y=l_b,yend=l_b,col=type))+
  geom_segment(data = plot_res_unsature,aes(x=as.numeric(cutoff)+0.05,xend=as.numeric(cutoff)+0.25,y=u_b,yend=u_b,col=type))+
  geom_point(data = plot_res_unsature,aes(x=as.numeric(cutoff)+0.15,y=beita,col=type))+
  geom_line(data = plot_res_unsature,aes(x=as.numeric(cutoff)+0.15,y=beita,col=type))+
  scale_color_manual(values = c("#4169E1","#32CD32"))+
  geom_hline(yintercept = 0,lwd=0.5,lty=2,col="#4B0082")+
  xlab("cutoff")+
  ylab("beita")+
  scale_x_continuous(breaks = seq(1,14,1),labels = seq(0,13,1))+
  theme_bw()
dev.off()







##############################################################
mt_ls_need = TG_long_ls

plot_res_long = lapply(1:length(mt_ls_need),function(kk){
  
  norm_sel = mt_expr %>% dplyr::select(any_of(mt_ls_need[[kk]])) %>% rownames_to_column("SUBJID")
  
  library(RobustRankAggreg)
  library(survival)
  library(survminer)
  method_choose="mean"
  glist<-list()
  i=2
  for (i in 2:ncol(norm_sel)){
    data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
    glist[[i-1]]<-data$SUBJID
    names(glist)[i-1]<-colnames(norm_sel)[i]
  }
  library(RobustRankAggreg)
  r<-rankMatrix(glist,full = TRUE)
  RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
  
  RRA = RRA %>% dplyr::rename(SUBJID=Name) %>% inner_join(clinic_C1[,c("SUBJID","R","AVAL1","CNSR1")],by = "SUBJID") %>% 
    dplyr::mutate(class=names(mt_ls_need)[kk])
  
  # cox<-coxph(Surv(AVAL1,CNSR1)~Score, data = RRA) %>% summary()
  # 
  # cox$conf.int
  # return(data.frame(class=names(mt_ls_need)[kk],HR=cox$conf.int[1],l_HR=cox$conf.int[3],u_HR=cox$conf.int[4],
  #                   p=cox$coefficients[5]))
  
  model =glm(as.factor(R)~Score,data = RRA,family = binomial)
  s=summary(model)
  confint(model)
  return(data.frame(class=names(mt_ls_need)[kk],beita=s$coefficients["Score",1],l_b=confint(model)[2,1],u_b=confint(model)[2,2]))
  # return(RRA)
  
}) %>% do.call(rbind,.)


mt_ls_need = TG_short_ls

plot_res_short = lapply(1:length(mt_ls_need),function(kk){
  
  norm_sel = mt_expr %>% dplyr::select(any_of(mt_ls_need[[kk]])) %>% rownames_to_column("SUBJID")
  
  library(RobustRankAggreg)
  library(survival)
  library(survminer)
  method_choose="mean"
  glist<-list()
  i=2
  for (i in 2:ncol(norm_sel)){
    data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
    glist[[i-1]]<-data$SUBJID
    names(glist)[i-1]<-colnames(norm_sel)[i]
  }
  library(RobustRankAggreg)
  r<-rankMatrix(glist,full = TRUE)
  RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
  
  RRA = RRA %>% dplyr::rename(SUBJID=Name) %>% inner_join(clinic_C1[,c("SUBJID","R","AVAL1","CNSR1")],by = "SUBJID") %>% 
    dplyr::mutate(class=names(mt_ls_need)[kk])
  

  model =glm(as.factor(R)~Score,data = RRA,family = binomial)
  s=summary(model)
  confint(model)
  return(data.frame(class=names(mt_ls_need)[kk],beita=s$coefficients["Score",1],l_b=confint(model)[2,1],u_b=confint(model)[2,2]))
  # return(RRA)
  
}) %>% do.call(rbind,.)



plot_res = rbind(plot_res_short,plot_res_long) %>% 
  dplyr::mutate(cutoff=str_extract(class,"[0-9]+$")) %>% 
  dplyr::mutate(cutoff=as.numeric(cutoff)) %>% 
  arrange(cutoff) %>% 
  dplyr::mutate(cutoff=factor(cutoff,levels = unique(cutoff),ordered = T)) %>% 
  dplyr::mutate(type=str_remove(class,"[0-9]+$"))

plot_res_short = plot_res %>% dplyr::filter(type=="TG_short")
plot_res_long = plot_res %>% dplyr::filter(type=="TG_long")

length(unique(plot_res$cutoff))
pdf("../../plot/TG_cutoff_logistic_long_C1.pdf",width = 8,height = 4)
ggplot()+
  geom_segment(data = plot_res_short,aes(x=as.numeric(cutoff)-0.15,xend=as.numeric(cutoff)-0.15,y=l_b,yend=u_b,col=type))+
  geom_segment(data = plot_res_short,aes(x=as.numeric(cutoff)-0.25,xend=as.numeric(cutoff)-0.05,y=l_b,yend=l_b,col=type))+
  geom_segment(data = plot_res_short,aes(x=as.numeric(cutoff)-0.25,xend=as.numeric(cutoff)-0.05,y=u_b,yend=u_b,col=type))+
  geom_point(data = plot_res_short,aes(x=as.numeric(cutoff)-0.15,y=beita,col=type))+
  geom_line(data = plot_res_short,aes(x=as.numeric(cutoff)-0.15,y=beita,col=type))+
  geom_segment(data = plot_res_long,aes(x=as.numeric(cutoff)+0.15,xend=as.numeric(cutoff)+0.15,y=l_b,yend=u_b,col=type))+
  geom_segment(data = plot_res_long,aes(x=as.numeric(cutoff)+0.05,xend=as.numeric(cutoff)+0.25,y=l_b,yend=l_b,col=type))+
  geom_segment(data = plot_res_long,aes(x=as.numeric(cutoff)+0.05,xend=as.numeric(cutoff)+0.25,y=u_b,yend=u_b,col=type))+
  geom_point(data = plot_res_long,aes(x=as.numeric(cutoff)+0.15,y=beita,col=type))+
  geom_line(data = plot_res_long,aes(x=as.numeric(cutoff)+0.15,y=beita,col=type))+
  scale_color_manual(values = c("#4169E1","#32CD32"))+
  geom_hline(yintercept = 0,lwd=0.5,lty=2,col="#4B0082")+
  xlab("cutoff")+
  ylab("beita")+
  scale_x_continuous(breaks = seq(1,20,1),labels = unique(plot_res$cutoff))+
  theme_bw()
dev.off()






mt_expr = norm_2_C6 %>% column_to_rownames("SUBJID")

load("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/VD/mt_enricher_TG.rdata")

mt_ls_need = TG_unsature_ls

plot_res_unsature = lapply(1:length(mt_ls_need),function(kk){
  
  norm_sel = mt_expr %>% dplyr::select(any_of(mt_ls_need[[kk]])) %>% rownames_to_column("SUBJID")
  
  library(RobustRankAggreg)
  library(survival)
  library(survminer)
  method_choose="mean"
  glist<-list()
  i=2
  for (i in 2:ncol(norm_sel)){
    data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
    glist[[i-1]]<-data$SUBJID
    names(glist)[i-1]<-colnames(norm_sel)[i]
  }
  library(RobustRankAggreg)
  r<-rankMatrix(glist,full = TRUE)
  RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
  
  RRA = RRA %>% dplyr::rename(SUBJID=Name) %>% inner_join(clinic_C6[,c("SUBJID","R","AVAL1","CNSR1")],by = "SUBJID") %>% 
    dplyr::mutate(class=names(mt_ls_need)[kk])
  

  model =glm(as.factor(R)~Score,data = RRA,family = binomial)
  s=summary(model)
  confint(model)
  return(data.frame(class=names(mt_ls_need)[kk],beita=s$coefficients["Score",1],l_b=confint(model)[2,1],u_b=confint(model)[2,2]))
  # return(RRA)
  
}) %>% do.call(rbind,.)


mt_ls_need = TG_sature_ls

plot_res_sature = lapply(1:length(mt_ls_need),function(kk){
  
  norm_sel = mt_expr %>% dplyr::select(any_of(mt_ls_need[[kk]])) %>% rownames_to_column("SUBJID")
  
  library(RobustRankAggreg)
  library(survival)
  library(survminer)
  method_choose="mean"
  glist<-list()
  i=2
  for (i in 2:ncol(norm_sel)){
    data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
    glist[[i-1]]<-data$SUBJID
    names(glist)[i-1]<-colnames(norm_sel)[i]
  }
  library(RobustRankAggreg)
  r<-rankMatrix(glist,full = TRUE)
  RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
  
  RRA = RRA %>% dplyr::rename(SUBJID=Name) %>% inner_join(clinic_C6[,c("SUBJID","R","AVAL1","CNSR1")],by = "SUBJID") %>% 
    dplyr::mutate(class=names(mt_ls_need)[kk])
  

  model =glm(as.factor(R)~Score,data = RRA,family = binomial)
  s=summary(model)
  confint(model)
  return(data.frame(class=names(mt_ls_need)[kk],beita=s$coefficients["Score",1],l_b=confint(model)[2,1],u_b=confint(model)[2,2]))
  # return(RRA)
  
}) %>% do.call(rbind,.)



plot_res = rbind(plot_res_sature,plot_res_unsature) %>% 
  dplyr::mutate(cutoff=str_extract(class,"[0-9]+$")) %>% 
  dplyr::mutate(cutoff=as.numeric(cutoff)) %>% 
  arrange(cutoff) %>% 
  dplyr::mutate(cutoff=factor(cutoff,levels = unique(cutoff),ordered = T)) %>% 
  dplyr::mutate(type=str_remove(class,"[0-9]+$"))

plot_res_sature = plot_res %>% dplyr::filter(type=="TG_sature")
plot_res_unsature = plot_res %>% dplyr::filter(type=="TG_unsature")


pdf("../../plot/TG_cutoff_logistic_unsature_C6.pdf",width = 7,height = 4)
ggplot()+
  geom_segment(data = plot_res_sature,aes(x=as.numeric(cutoff)-0.15,xend=as.numeric(cutoff)-0.15,y=l_b,yend=u_b,col=type))+
  geom_segment(data = plot_res_sature,aes(x=as.numeric(cutoff)-0.25,xend=as.numeric(cutoff)-0.05,y=l_b,yend=l_b,col=type))+
  geom_segment(data = plot_res_sature,aes(x=as.numeric(cutoff)-0.25,xend=as.numeric(cutoff)-0.05,y=u_b,yend=u_b,col=type))+
  geom_point(data = plot_res_sature,aes(x=as.numeric(cutoff)-0.15,y=beita,col=type))+
  geom_line(data = plot_res_sature,aes(x=as.numeric(cutoff)-0.15,y=beita,col=type))+
  geom_segment(data = plot_res_unsature,aes(x=as.numeric(cutoff)+0.15,xend=as.numeric(cutoff)+0.15,y=l_b,yend=u_b,col=type))+
  geom_segment(data = plot_res_unsature,aes(x=as.numeric(cutoff)+0.05,xend=as.numeric(cutoff)+0.25,y=l_b,yend=l_b,col=type))+
  geom_segment(data = plot_res_unsature,aes(x=as.numeric(cutoff)+0.05,xend=as.numeric(cutoff)+0.25,y=u_b,yend=u_b,col=type))+
  geom_point(data = plot_res_unsature,aes(x=as.numeric(cutoff)+0.15,y=beita,col=type))+
  geom_line(data = plot_res_unsature,aes(x=as.numeric(cutoff)+0.15,y=beita,col=type))+
  scale_color_manual(values = c("#4169E1","#32CD32"))+
  geom_hline(yintercept = 0,lwd=0.5,lty=2,col="#4B0082")+
  xlab("cutoff")+
  ylab("beita")+
  scale_x_continuous(breaks = seq(1,14,1),labels = seq(0,13,1))+
  theme_bw()
dev.off()







##############################################################
mt_ls_need = TG_long_ls

plot_res_long = lapply(1:length(mt_ls_need),function(kk){
  
  norm_sel = mt_expr %>% dplyr::select(any_of(mt_ls_need[[kk]])) %>% rownames_to_column("SUBJID")
  
  library(RobustRankAggreg)
  library(survival)
  library(survminer)
  method_choose="mean"
  glist<-list()
  i=2
  for (i in 2:ncol(norm_sel)){
    data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
    glist[[i-1]]<-data$SUBJID
    names(glist)[i-1]<-colnames(norm_sel)[i]
  }
  library(RobustRankAggreg)
  r<-rankMatrix(glist,full = TRUE)
  RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
  
  RRA = RRA %>% dplyr::rename(SUBJID=Name) %>% inner_join(clinic_C6[,c("SUBJID","R","AVAL1","CNSR1")],by = "SUBJID") %>% 
    dplyr::mutate(class=names(mt_ls_need)[kk])
  

  model =glm(as.factor(R)~Score,data = RRA,family = binomial)
  s=summary(model)
  confint(model)
  return(data.frame(class=names(mt_ls_need)[kk],beita=s$coefficients["Score",1],l_b=confint(model)[2,1],u_b=confint(model)[2,2]))
  # return(RRA)
  
}) %>% do.call(rbind,.)


mt_ls_need = TG_short_ls

plot_res_short = lapply(1:length(mt_ls_need),function(kk){
  
  norm_sel = mt_expr %>% dplyr::select(any_of(mt_ls_need[[kk]])) %>% rownames_to_column("SUBJID")
  
  library(RobustRankAggreg)
  library(survival)
  library(survminer)
  method_choose="mean"
  glist<-list()
  i=2
  for (i in 2:ncol(norm_sel)){
    data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
    glist[[i-1]]<-data$SUBJID
    names(glist)[i-1]<-colnames(norm_sel)[i]
  }
  library(RobustRankAggreg)
  r<-rankMatrix(glist,full = TRUE)
  RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
  
  RRA = RRA %>% dplyr::rename(SUBJID=Name) %>% inner_join(clinic_C6[,c("SUBJID","R","AVAL1","CNSR1")],by = "SUBJID") %>% 
    dplyr::mutate(class=names(mt_ls_need)[kk])
  

  model =glm(as.factor(R)~Score,data = RRA,family = binomial)
  s=summary(model)
  confint(model)
  return(data.frame(class=names(mt_ls_need)[kk],beita=s$coefficients["Score",1],l_b=confint(model)[2,1],u_b=confint(model)[2,2]))
  # return(RRA)
  
}) %>% do.call(rbind,.)



plot_res = rbind(plot_res_short,plot_res_long) %>% 
  dplyr::mutate(cutoff=str_extract(class,"[0-9]+$")) %>% 
  dplyr::mutate(cutoff=as.numeric(cutoff)) %>% 
  arrange(cutoff) %>% 
  dplyr::mutate(cutoff=factor(cutoff,levels = unique(cutoff),ordered = T)) %>% 
  dplyr::mutate(type=str_remove(class,"[0-9]+$"))

plot_res_short = plot_res %>% dplyr::filter(type=="TG_short")
plot_res_long = plot_res %>% dplyr::filter(type=="TG_long")

length(unique(plot_res$cutoff))
pdf("../../plot/TG_cutoff_logistic_long_C6.pdf",width = 8,height = 4)
ggplot()+
  geom_segment(data = plot_res_short,aes(x=as.numeric(cutoff)-0.15,xend=as.numeric(cutoff)-0.15,y=l_b,yend=u_b,col=type))+
  geom_segment(data = plot_res_short,aes(x=as.numeric(cutoff)-0.25,xend=as.numeric(cutoff)-0.05,y=l_b,yend=l_b,col=type))+
  geom_segment(data = plot_res_short,aes(x=as.numeric(cutoff)-0.25,xend=as.numeric(cutoff)-0.05,y=u_b,yend=u_b,col=type))+
  geom_point(data = plot_res_short,aes(x=as.numeric(cutoff)-0.15,y=beita,col=type))+
  geom_line(data = plot_res_short,aes(x=as.numeric(cutoff)-0.15,y=beita,col=type))+
  geom_segment(data = plot_res_long,aes(x=as.numeric(cutoff)+0.15,xend=as.numeric(cutoff)+0.15,y=l_b,yend=u_b,col=type))+
  geom_segment(data = plot_res_long,aes(x=as.numeric(cutoff)+0.05,xend=as.numeric(cutoff)+0.25,y=l_b,yend=l_b,col=type))+
  geom_segment(data = plot_res_long,aes(x=as.numeric(cutoff)+0.05,xend=as.numeric(cutoff)+0.25,y=u_b,yend=u_b,col=type))+
  geom_point(data = plot_res_long,aes(x=as.numeric(cutoff)+0.15,y=beita,col=type))+
  geom_line(data = plot_res_long,aes(x=as.numeric(cutoff)+0.15,y=beita,col=type))+
  scale_color_manual(values = c("#4169E1","#32CD32"))+
  geom_hline(yintercept = 0,lwd=0.5,lty=2,col="#4B0082")+
  xlab("cutoff")+
  ylab("beita")+
  scale_x_continuous(breaks = seq(1,20,1),labels = unique(plot_res$cutoff))+
  theme_bw()
dev.off()






