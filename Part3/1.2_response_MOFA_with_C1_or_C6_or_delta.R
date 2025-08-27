#!/usr/bin/Rscript
# BiocManager::install("MOFA2")


rm(list = ls())
#### Load packages ####
library(reticulate)
use_condaenv("mofapy2",required = TRUE)
library(Seurat)
library(tidyverse)
library(scCustomize)
library(data.table)
library(MOFA2)


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


norm_3_C1 = norm_2_C1 %>% dplyr::filter(SUBJID %in% norm_3_process$SUBJID)
norm_3_C6 = norm_2_C6 %>% dplyr::filter(SUBJID %in% norm_3_process$SUBJID)

clinic_need = clinic_1[,c(2,6,7,9,11,12,13,14,18,19,20,21,25,26,35)] %>% distinct(.,.keep_all = T) %>% 
  dplyr::filter(SUBJID %in% norm_3_process$SUBJID)
clinic_info_1 = clinic_need %>% dplyr::select(1:12) %>% column_to_rownames("SUBJID")
clinic_info_3 = clinic_need %>% dplyr::select(1,15) %>% column_to_rownames("SUBJID")
meta.data <- cbind(clinic_info_1,clinic_info_3)





m = list()
m[[1]]=norm_3_C1 %>% column_to_rownames("SUBJID") %>% as.matrix() %>% t()
rownames(m[[1]])=str_c("C1",1:nrow(m[[1]]),sep = "_")
m[[2]]=norm_3_C6 %>% column_to_rownames("SUBJID") %>% as.matrix() %>% t()
rownames(m[[2]])=str_c("C6",1:nrow(m[[2]]),sep = "_")
m[[3]]=norm_3_process %>% column_to_rownames("SUBJID") %>% as.matrix() %>% t()
rownames(m[[3]])=str_c("delta",1:nrow(m[[3]]),sep = "_")
lapply(m,dim)

# names(m)=c("spliced","unspliced")
# groups = dt$group
MOFAobject <- create_mofa(m)

plot_data_overview(MOFAobject)
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)


MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

dir.create("./MOFA/")
outfile = file.path("./MOFA/all_MOFA.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile)



model <- load_model(outfile)


sum(model@dimensions$N)
dt=meta.data %>% rownames_to_column("SUBJID")
dim(dt)
dt = dt %>% dplyr::rename(sample=SUBJID)
dt = dt %>% dplyr::rename(condition=R)
samples_metadata(model) <- dt

head(model@samples_metadata, n=3)

head(model@cache$variance_explained$r2_total[[1]])

head(model@cache$variance_explained$r2_per_factor[[1]])

plot_variance_explained(model, x="view", y="factor")


print(p)




variance <- get_variance_explained(model,factors = "all",as.data.frame = T)
variance_total <- variance[[2]]
variance_total$view=c("C1","C6","delta")
variance_sub <- variance[[1]]
variance_sub$view = c(rep("C1",15),rep("C6",15),rep("delta",15))
variance_sub = variance_sub %>% dplyr::filter(factor %in% str_c("Factor",1:5,sep = ""))
variance_total = variance_total %>% dplyr::mutate(factor="total")
variance_all = rbind(variance_total,variance_sub)
variance_all$factor = factor(variance_all$factor,levels = c("total",str_c("Factor",1:5,sep = "")),ordered = T)
library(ggplot2)
variance_all = variance_all %>% dplyr::filter(!factor=="total")

pdf("./MOFA/variance_of_explanation.pdf",width = 4,height = 3)
ggplot(variance_all)+
  geom_bar(aes(x =factor,y=value,fill=view),stat = "identity",position = position_dodge(width = 0.8))+
  ylab("variance of explanation(%)")+
  scale_fill_manual(values = wes_palette("GrandBudapest1",3))+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  RotatedAxis()
dev.off()


