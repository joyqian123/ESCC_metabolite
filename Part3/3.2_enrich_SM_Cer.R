rm(list = ls())
library(tidyverse)
library(Seurat)
setwd("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/")

load("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/VD/mt_list_from_data.rdata")

load("../1_preclean/data/data_C1_preclean.rdata")
mt_info = data_in_C1 %>% dplyr::select(Compounds,`Class I`,`Class II`) %>% 
  dplyr::rename(mt=Compounds)
class_ID_1 = mt_info %>% dplyr::select(mt,`Class I`) %>% dplyr::rename(class=colnames(.)[2])
class_ID_2 = mt_info %>% dplyr::select(mt,`Class II`) %>% dplyr::rename(class=colnames(.)[2])

load("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/VD/only_response_VD_res.rdata")
plot1 = vd.res_def_C1 %>% dplyr::select(R) %>% rownames_to_column("mt")  %>% dplyr::mutate(data="C1")
plot2 = vd.res_def_C6 %>% dplyr::select(R) %>% rownames_to_column("mt")  %>% dplyr::mutate(data="C6")
plot3 = vd.res_def_delta %>% dplyr::select(R) %>% rownames_to_column("mt")  %>% dplyr::mutate(data="C6-C1")
sel_mt = rbind(plot1,plot2,plot3) %>% dplyr::filter(R>0) %>% pull(mt)

merge = rbind(plot1,plot2,plot3) %>% dplyr::filter(mt %in% sel_mt) %>% 
  pivot_wider(id_cols = mt,names_from = "data",values_from = "R")

colnames(FC_res_C1)[2:4]=str_c(colnames(FC_res_C1)[2:4],"C1",sep = "_")
colnames(FC_res_C6)[2:4]=str_c(colnames(FC_res_C6)[2:4],"C6",sep = "_")
colnames(FC_res_delta)[4:6]=str_c(colnames(FC_res_delta)[4:6],"delta",sep = "_")


merge = merge %>% inner_join(FC_res_C1,by = "mt") %>% inner_join(FC_res_C6,by = "mt") %>% 
  inner_join(FC_res_delta,by = "mt")



SM_ls = class_ID_2 %>% dplyr::filter(class=="SM") %>% 
  dplyr::mutate(tmp=str_remove(mt,"SM\\(")) %>% 
  separate(.,col="tmp",into=c("A","B"),sep="\\/") %>% 
  dplyr::mutate(B=str_remove(B,"\\)")) %>% inner_join(merge,by = "mt")



SP_ls = mt_info %>% dplyr::filter(`Class I`=="SP") %>% inner_join(merge,by = "mt") %>% 
  pivot_longer(cols=c("t_C6","t_C1"),names_to = "time",values_to = "t_value")

ggplot(SP_ls)+
  geom_boxplot(aes(x=`Class II`,y=t_value,col=time))+
  RotatedAxis()


kk=1
SP_mt_ls = lapply(1:length(unique(SP_ls$`Class II`)),function(kk){
  tmp = SP_ls %>% dplyr::filter(`Class II`==unique(SP_ls$`Class II`)[kk]) %>% pull(mt)
  return(tmp)
})
names(SP_mt_ls)=unique(SP_ls$`Class II`)

SP_mt_ls <- Filter(function(x) length(x) > 1, SP_mt_ls)



SP_ls_all = SP_ls %>% dplyr::select(mt,`Class I`) %>% dplyr::rename(class=colnames(.)[2])
#######################################################################
#########################################################################
SP_ls_only_SM_Cer = SP_ls %>% dplyr::select(mt,`Class II`) %>% 
  dplyr::mutate(class = str_extract(`Class II`,"Cer")) %>% 
  dplyr::mutate(class = ifelse(is.na(class),`Class II`,class)) %>% 
  dplyr::filter(!class=="SPH")

SP_ls_only_SM_Cer = rbind(SP_ls_all,SP_ls_only_SM_Cer[,c(1,3)])

geneList= merge$t_C1
names(geneList) = merge$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=SP_ls_only_SM_Cer[,c(2,1)], 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_SP = egmt@result
p1=gseaNb(object = egmt,
          geneSetID = c("SM","Cer","SP"),
          newGsea = T,
          addPval = T,
          pvalX = 0.5,pvalY = 0.75,
          pCol = 'black',
          pHjust = 1,
          subPlot = 2)

geneList= merge$t_C6
names(geneList) = merge$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=SP_ls_only_SM_Cer[,c(2,1)], 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_SP = egmt@result
p2=gseaNb(object = egmt,
          geneSetID = c("SM","Cer","SP"),
          newGsea = T,
          addPval = T,
          pvalX = 0.5,pvalY = 0.75,
          pCol = 'black',
          pHjust = 1,
          subPlot = 2)


library(patchwork)
pdf("./plot/SP/C1_C6_SM_Cer_SP_gsea.pdf",width = 10,height = 6)
wrap_plots(p1,p2,nrow = 1)
dev.off()





#######################################################################
#########################################################################
geneList= merge$t_C1
names(geneList) = merge$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=SP_ls[,c(3,1)], 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_SP_C1 = egmt@result %>% dplyr::mutate(time="C1")

geneList= merge$t_C6
names(geneList) = merge$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
set.seed(12345)
egmt <- GSEA(geneList, TERM2GENE=SP_ls[,c(3,1)], 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_SP_C6 = egmt@result %>% dplyr::mutate(time="C6")


final = rbind(egmt_res_SP_C1,egmt_res_SP_C6)

final_plot = final %>% dplyr::select(ID,NES,pvalue,time)
final_NES = final_plot %>% pivot_wider(id_cols = ID,names_from = time,values_from = NES) %>% 
  dplyr::mutate(C1=as.numeric(C1),C6=as.numeric(C6)) %>% arrange(desc(C6)) %>%  
  column_to_rownames("ID") %>% as.matrix()
rownames(final_NES) = factor(rownames(final_NES),levels = rownames(final_NES),ordered = T)
final_p = final_plot %>% pivot_wider(id_cols = ID,names_from = time,values_from = pvalue)%>% 
  dplyr::mutate(C1=as.numeric(C1),C6=as.numeric(C6)) %>% 
  dplyr::mutate(ID=factor(ID,levels = rownames(final_NES),ordered = T)) %>% 
  arrange(ID) %>% 
  dplyr::mutate(C1=ifelse(C1<0.001,"***",ifelse(C1<0.01,"**",ifelse(C1<0.05,"*","")))) %>% 
  dplyr::mutate(C6=ifelse(C6<0.001,"***",ifelse(C6<0.01,"**",ifelse(C6<0.05,"*","")))) %>% 
  column_to_rownames("ID") %>% as.matrix()


# m2 = matrix(rnorm(50*10), nrow = 50)
kk=1
library(enrichplot)
lt6 = lapply(1:length(rownames(final_NES)),function(kk){
  geneList= merge$t_C6
  names(geneList) = merge$mt
  geneList=sort(geneList,decreasing = T)
  head(geneList)
  library(clusterProfiler)
  tmp = SP_ls[,c(3,1)] %>% dplyr::filter(`Class II`==rownames(final_NES)[kk])
  egmt <- GSEA(geneList, TERM2GENE=tmp, 
               minGSSize = 1,
               pvalueCutoff = 0.99,
               verbose=FALSE)
  gsea_res <- gseaplot2(egmt,geneSetID = unique(tmp$`Class II`),color = "purple")[[1]]
  dt = gsea_res$data$runningScore
  
})
names(lt6)=rownames(final_NES)

lt1 = lapply(1:length(rownames(final_NES)),function(kk){
  geneList= merge$t_C1
  names(geneList) = merge$mt
  geneList=sort(geneList,decreasing = T)
  head(geneList)
  library(clusterProfiler)
  tmp = SP_ls[,c(3,1)] %>% dplyr::filter(`Class II`==rownames(final_NES)[kk])
  egmt <- GSEA(geneList, TERM2GENE=tmp, 
               minGSSize = 1,
               pvalueCutoff = 0.99,
               verbose=FALSE)
  gsea_res <- gseaplot2(egmt,geneSetID = unique(tmp$`Class II`),color = "purple")[[1]]
  dt = gsea_res$data$runningScore
  
})
names(lt1)=rownames(final_NES)
library(ComplexHeatmap)
row_anno <- rowAnnotation(GSEA_C1 = anno_horizon(lt1),
                          GSEA_C6 = anno_horizon(lt6))


color = c(colorRampPalette(colors = c("navy","white"))(100),colorRampPalette(colors = c("white","red"))(300))
pdf("./plot/SP/C1_C6_SP_subtype.pdf",width = 6,height = 6)
Heatmap(final_NES,cluster_columns = FALSE,cluster_rows = FALSE,border = TRUE,
        rect_gp = gpar(col = "white"),right_annotation = row_anno,
        col = color,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(final_p[i, j], x, y, gp = gpar(fontsize = 10))
        })
dev.off()
