######################################################################################
#T cell Analysis----
######################################################################################
rm(list = ls())
library(Seurat)
library(dplyr)
library(SAVER)
#library(ProjecTILs)
library(future)
library(patchwork)
library(harmony)
library(clusterProfiler)
library(pheatmap)
# library(CytoTRACE)
library(clustree)
library(monocle3)
library(ggsci)
library(VISION)
library(MEGENA)
library(GSVA)

# devtools::install_local("~/github/CytoTRACE_0.3.3.tar.gz")
# devtools::install_github("cole-trapnell-lab/monocle3")
# plan("multiprocess", workers = 20)
# plan()

# rm(list=ls())
setwd("~/help_for_others/DYQ/scRNA//")
dir.create("11_Tcell/")
scRNA<-
  qs::qread("4_Tcell/scRNA_CD8Tcell_SCT_AnnotatedWith0.8.qs")


DimPlot(scRNA,group.by = "Main_Celltype",reduction = "umap")  
FeaturePlot(scRNA,features = c("Cd3e","Cd4","Cd8a"))
DimPlot(scRNA,group.by = "subtype",reduction = "umap")  





pdf("./11_Tcell/Dimplot_subtype.pdf",width = 4.5,height = 3.2)
scRNA$subtype = factor(scRNA$subtype,levels = c("CD8-Tn","CD8-Tm","CD8-Temra","CD8-T-ISG","CD8-Teff","CD8-Tprolif","CD8-Tex"),ordered = T)
DimPlot(scRNA,reduction = "umap",group.by = "subtype",pt.size = 0.1,label = T,label.size = 2)+
  scale_color_manual(values = rev(c("#2d467b","#315098","#4077d0","#5591dc","#75aee5","#a2cbee","#c8dff5")))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()
pdf("./11_Tcell/Dimplot_subtype_split.pdf",width = 8,height = 5.5)
scRNA$group = factor(scRNA$group,levels = c("PP_tumor","SPP_tumor","IPP_tumor","PP_pbmc","SPP_pbmc","IPP_pbmc"),ordered = T)
DimPlot(scRNA,reduction = "umap",group.by = "subtype",pt.size = 0.15,split.by = "group",ncol=3)+
  scale_color_manual(values = rev(c("#2d467b","#315098","#4077d0","#5591dc","#75aee5","#a2cbee","#c8dff5")))+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()



CD8 = scRNA@meta.data
prop_tb = table(CD8$sample,CD8$subtype) %>% data.frame() %>% 
  group_by(Var1) %>% 
  dplyr::mutate(total=sum(Freq)) %>% dplyr::mutate(prop=Freq/total)

prop_tb$group = str_extract(prop_tb$Var1,"^[A-Z]+")
prop_tb$origin = str_extract(prop_tb$Var1,"[a-z]+$")
prop_tb$group = factor(prop_tb$group,levels = c("PP","SPP","IPP"),ordered = T)
pdf("./11_Tcell/sample_prop_change_CD8_subtype.pdf",width = 12,height = 6)
ggplot(prop_tb,aes(x=group,y=prop,col=group))+
  geom_boxplot()+
  geom_jitter()+
  # geom_line(aes(group=m_group))+
  scale_color_manual(values = c("#95295f","#1d3f7b","#88c4e8"))+
  stat_compare_means(comparisons = list(c("PP","SPP"),c("PP","IPP")),method = "t.test")+
  facet_wrap(origin~Var2,scales = "free",nrow = 2)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black",linewidth = 0.5), # 添加 x 和 y 轴线
        axis.ticks = element_line(color = "black"),
        # legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()


prop_tb_group = table(scRNA$group,scRNA$subtype) %>% data.frame() %>% 
  group_by(Var1) %>% 
  dplyr::mutate(total=sum(Freq)) %>% dplyr::mutate(prop=Freq/total)
prop_tb_group$Var1 = factor(prop_tb_group$Var1,levels = c("PP_tumor","SPP_tumor","IPP_tumor","PP_pbmc","SPP_pbmc","IPP_pbmc"),ordered = T)

library(ggalluvial)

pdf("./11_Tcell/all_prop_change.pdf",width = 8,height = 6)
ggplot(prop_tb_group, 
       aes(x = Var1, y = prop,
           stratum = Var2, 
           alluvium = Var2, 
           fill = Var2, 
           label = Var2)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(alpha = 0.8,lwd=0.3) +
  scale_fill_manual(values = rev(c("#2d467b","#315098","#4077d0","#5591dc","#75aee5","#a2cbee","#c8dff5")))+
  labs(x="",y="",fill="Main_Celltype")+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.border = element_blank(),
        axis.line = element_blank(), # 添加 x 和 y 轴线 
        axis.ticks = element_line(color = "black"),
        # legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  RotatedAxis()
dev.off()

# rev(c("#245c38","#215272","#349e9e","#56c497","#7ce394","#D0F4D0","#b2cdcc")


prop_tb_group = prop_tb_group %>% dplyr::filter(Var1 %in% c("PP_tumor","SPP_tumor","IPP_tumor"))
pdf("./11_Tcell/all_prop_change_only_tumor.pdf",width = 5,height = 4)
ggplot(prop_tb_group, 
       aes(x = Var1, y = prop,
           stratum = Var2, 
           alluvium = Var2, 
           fill = Var2, 
           label = Var2)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(alpha = 0.8,lwd=0.3) +
  scale_fill_manual(values = rev(c("#2d467b","#315098","#4077d0","#5591dc","#75aee5","#a2cbee","#c8dff5")))+
  labs(x="",y="",fill="Main_Celltype")+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.border = element_blank(),
        axis.line = element_blank(), # 添加 x 和 y 轴线 
        axis.ticks = element_line(color = "black"),
        # legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  RotatedAxis()
dev.off()


prop_tb_group = table(scRNA$group,scRNA$subtype) %>% data.frame() %>% 
  group_by(Var1) %>% 
  dplyr::mutate(total=sum(Freq)) %>% dplyr::mutate(prop=Freq/total)
prop_tb_group$Var1 = factor(prop_tb_group$Var1,levels = c("PP_tumor","SPP_tumor","IPP_tumor","PP_pbmc","SPP_pbmc","IPP_pbmc"),ordered = T)

prop_tb_group = prop_tb_group %>% dplyr::filter(Var1 %in% c("PP_pbmc","SPP_pbmc","IPP_pbmc"))
pdf("./11_Tcell/all_prop_change_only_pbmc.pdf",width = 5,height = 4)
ggplot(prop_tb_group, 
       aes(x = Var1, y = prop,
           stratum = Var2, 
           alluvium = Var2, 
           fill = Var2, 
           label = Var2)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(alpha = 0.8,lwd=0.3) +
  scale_fill_manual(values = rev(c("#2d467b","#315098","#4077d0","#5591dc","#75aee5","#a2cbee","#c8dff5")))+
  labs(x="",y="",fill="Main_Celltype")+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.border = element_blank(),
        axis.line = element_blank(), # 添加 x 和 y 轴线 
        axis.ticks = element_line(color = "black"),
        # legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  RotatedAxis()
dev.off()






################################
#####data imputation
expr.in.cells <- rowSums(scRNA[["RNA"]]$counts > 0)
genes <- names(expr.in.cells)[expr.in.cells >= 10]
expr.mat <- scRNA[["RNA"]]$data[genes, ] %>% as.matrix()# log normalized matrix
dim(expr.mat)
# RcppML::getRcppMLthreads()
RcppML::setRcppMLthreads(5) ## 使用5个线程
# k << number of genes
# Run NMF multiple times with different seeds and average the results
seeds <- c(8,16,32,1024, 2048)
imputed.mats <- list()

for(i in seq_along(seeds)) {
  model <- RcppML::nmf(expr.mat, k = 100, verbose = T, seed = seeds[i])
  model$w %>% dim
  model$h %>% dim
  model$d # a vector, length = k
  
  imputed.mats[[i]] <- model$w %*% diag(model$d) %*% model$h
}

# Calculate average imputed matrix
imputed.mat <- Reduce('+', imputed.mats) / length(seeds)

imputed.mat %>% dim
colnames(imputed.mat) <- colnames(expr.mat)
rownames(imputed.mat) <- rownames(expr.mat)
scRNA[["imputed"]] <- CreateAssay5Object(data = imputed.mat)

qs::qsave(scRNA,"./11_Tcell//CD8T_scRNA_imputation.qs")





scRNA = qs::qread("./11_Tcell/CD8T_scRNA_imputation.qs")
library(tidyverse)
DefaultAssay(scRNA)="imputed"
# VlnPlot(scRNA,group.by = "group",features = c("Ctla4","Pdcd1","Havcr2","Lag3","Tox","Tigit","Bach2","Layn","Entpd1","Dusp4","Itgae"))
tb = FetchData(scRNA,
               vars = c("group","Pdcd1","Havcr2","Lag3","Tigit","Layn"),
               layer = "data") %>% 
  pivot_longer(cols = !c("group"),names_to = "gene",values_to = "expr") %>% 
  dplyr::filter(grepl("tumor",group)) %>% 
  dplyr::mutate(group=as.character(group))
tb$group = factor(tb$group,levels = c("PP_tumor","SPP_tumor","IPP_tumor"),ordered = T)


median_tb <- tb %>%
  group_by(group,gene) %>%
  dplyr::summarise(median_expr = median(expr, na.rm = TRUE)) %>%
  ungroup()

library(ggpubr)
pdf("./11_Tcell/ex_marker.pdf",width = 10,height = 3)
ggplot(tb,aes(x=group,y=expr,col=group,fill=group))+
  geom_violin(alpha=0.5)+
  stat_compare_means(comparisons = list(c("PP_tumor","SPP_tumor"),c("PP_tumor",'IPP_tumor')))+
  scale_color_manual(values = c("#95295f","#1d3f7b","#88c4e8"))+
  scale_fill_manual(values = c("#95295f","#1d3f7b","#88c4e8"))+
  geom_segment(data = median_tb,  # 使用提前计算的中位数数据
               aes(x= as.numeric(as.factor(group)) - 0.1,  # 调整x的值来控制线段的起点
                   xend = as.numeric(as.factor(group)) + 0.1,  # 调整xend的值来控制线段的终点
                   y = median_expr,
                   yend = median_expr),
               linewidth=1,lty=1) +
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black",linewidth = 0.5), # 添加 x 和 y 轴线
        axis.ticks = element_line(color = "black"),
        # legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  facet_wrap(~gene,nrow=1,scale="free_y") +
  RotatedAxis()
dev.off()




tb = FetchData(scRNA,
               vars = c("group","Eomes","Cd28","Cd44","Cd69","Fgfbp2","Cx3cr1","Fcgr3a", "Klrg1","Cxcr5","Cxcr4","Cxcr3","Ccr4",
                        "Gzma","Cd39","Ifng","Il2","Tnf"),layer = "data") %>% 
  pivot_longer(cols = !c("group"),names_to = "gene",values_to = "expr")%>% 
  dplyr::filter(grepl("tumor",group))%>% 
  dplyr::mutate(group=as.character(group))

tb$group = factor(tb$group,levels = c("PP_tumor","SPP_tumor","IPP_tumor"),ordered = T)


median_tb <- tb %>%
  group_by(group,gene) %>%
  dplyr::summarise(median_expr = median(expr, na.rm = TRUE)) %>%
  ungroup()

library(ggpubr)
pdf("./11_Tcell/eff_marker.pdf",width = 14,height = 5)
ggplot(tb,aes(x=group,y=expr,col=group,fill=group))+
  geom_violin(alpha=0.5)+
  stat_compare_means(comparisons = list(c("PP_tumor","SPP_tumor"),c("PP_tumor",'IPP_tumor')))+
  scale_color_manual(values = c("#95295f","#1d3f7b","#88c4e8"))+
  scale_fill_manual(values = c("#95295f","#1d3f7b","#88c4e8"))+
  geom_segment(data = median_tb,  # 使用提前计算的中位数数据
               aes(x= as.numeric(as.factor(group)) - 0.1,  # 调整x的值来控制线段的起点
                   xend = as.numeric(as.factor(group)) + 0.1,  # 调整xend的值来控制线段的终点
                   y = median_expr,
                   yend = median_expr),
               linewidth=1,lty=1) +
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black",linewidth = 0.5), # 添加 x 和 y 轴线
        axis.ticks = element_line(color = "black"),
        # legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  facet_wrap(~gene,nrow=2,scale="free_y") +
  RotatedAxis()
dev.off()





