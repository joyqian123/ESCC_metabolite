######################################################################################
rm(list = ls())
library(Seurat)
library(SAVER)
# library(ProjecTILs)
library(future)
library(patchwork)
library(harmony)
library(cluster)
library(pheatmap)
library(tidyverse)
# source("Resource/function.R")
# plan("multiprocess", workers = 10)
# plan()
setwd("~/help_for_others/DYQ/scRNA/")
# rm(list=ls())
scRNA <- qs::qread("./9_merge/scRNA_final_annotated.qs")
table(scRNA$orig.ident)
table(scRNA$Sub_Celltype)
table(scRNA$subtype)

meta = scRNA@meta.data




#####################################################################密度图
########################################################################
#####################################################################密度图
########################################################################
umap_data = FetchData(scRNA,vars = c("umap_1","umap_2","group"))
library(scico)
pdf("./9_merge/2_dimplot_density.pdf",width = 9,height = 6)
ggplot(umap_data, aes(x = umap_1, y = umap_2)) +
  # stat_density2d(aes(fill = ..level..), geom = "polygon", bins = 20) +
  geom_bin2d(bins = 80, aes(fill = ..count..)) +
  # scale_fill_gradient2(low = ) + # 使用viridis颜色调色板
  scale_fill_scico(palette = 'vik')+
  # theme_minimal() +
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  labs(fill = "Density")+
  facet_wrap(~group,nrow=2)
dev.off()






key_marker = list(
  B = c("Ms4a1","Cd79a"),
  `T` = c("Ptprc","Cd3e","Cd3d","Cd4","Cd8a","Cd8b1"),
  NK = c("Nkg7","Ncr1","Klrd1"),
  cDC = c("Xcr1","Clec9a","Itgax","Fscn1","Ccr7"),
  mast = c("Ms4a2","Gata2"),
  Mono_Macro = c("Cd14","Cd68","Csf1r","Mafb"),
  Neutrophil = c("S100a9","S100a8"),
  Fibroblast = c("Col1a1","Acta2"),
  Endothelium =  c("Flt1","Pecam1","Ramp2")
)

pdf("./9_merge/3.marker_main_cluster.pdf",width = 8,height = 8)
DotPlot(scRNA,features = unique(unlist(key_marker)),group.by = "Main_Celltype")+
  scale_color_gradient2(low = "white",mid = "#ADD8E6",high = "#6A5ACD")+
  coord_flip()+
  RotatedAxis()
dev.off()



###################################################################
###################################################################



prop_tb = table(scRNA$sample,scRNA$group,scRNA$Main_Celltype) %>% data.frame() %>% 
  group_by(Var1) %>% dplyr::filter(Freq>0) %>% 
  dplyr::mutate(total=sum(Freq)) %>% dplyr::mutate(prop=Freq/total)
prop_tb$Var2 = factor(prop_tb$Var2,levels = c("PP_tumor","SPP_tumor","IPP_tumor","PP_pbmc","SPP_pbmc","IPP_pbmc"),ordered = T)
prop_tb$Var1 = as.character(prop_tb$Var1)
order = prop_tb %>% arrange(Var2) %>% distinct(Var1,.keep_all = T) 
prop_tb$Var1 = factor(prop_tb$Var1,levels = unique(prop_tb$Var1),ordered = T)

prop_tb_group = table(scRNA$group,scRNA$Main_Celltype) %>% data.frame() %>% 
  group_by(Var1) %>% 
  dplyr::mutate(total=sum(Freq)) %>% dplyr::mutate(prop=Freq/total)

library(ggalluvial)
pdf("./9_merge//4.all_prop_change.pdf",width = 8,height = 6)
ggplot(prop_tb_group, 
       aes(x = Var1, y = prop,
           stratum = Var2, 
           alluvium = Var2, 
           fill = Var2, 
           label = Var2)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(alpha = 0.8,lwd=0.3) +
  scale_fill_manual(values = c("#bdc3d2","#88c4e8","#0074b3","#aeb0d8",
                                "#e5ce81","#f47720","#f7bfa4","#e6966a",
                                
                                "#aedacb","#4eb69e","#7c6aa4"))+
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



pdf("./9_merge//5.all_prop_change_split.pdf",width = 8,height = 6)
ggplot(prop_tb, 
       aes(x = Var1, y = prop,
           fill = Var3, 
           label = Var3)) +
  geom_col(position = "stack",alpha = 0.8) +
  scale_fill_manual(values = c("#bdc3d2","#88c4e8","#0074b3","#aeb0d8",
                               "#e5ce81","#f47720","#f7bfa4","#e6966a",
                               
                               "#aedacb","#4eb69e","#7c6aa4"))+
  labs(x="",y="",fill="Sub_Celltype")+
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



library(ggalluvial)
pdf("./9_merge//4.all_prop_change_pbmc.pdf",width = 6,height = 6)
prop_tb_group_1 = prop_tb_group %>% dplyr::filter(grepl("pbmc",Var1))
ggplot(prop_tb_group_1, 
       aes(x = Var1, y = prop,
           stratum = Var2, 
           alluvium = Var2, 
           fill = Var2, 
           label = Var2)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(alpha = 0.8,lwd=0.3) +
  scale_fill_manual(values = c("#bdc3d2","#88c4e8","#0074b3","#aeb0d8",
                               "#e5ce81","#f47720","#f7bfa4","#e6966a",
                               
                               "#aedacb","#4eb69e","#7c6aa4"))+
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



library(ggalluvial)
pdf("./9_merge//4.all_prop_change_tumor.pdf",width = 6,height = 6)
prop_tb_group_1 = prop_tb_group %>% dplyr::filter(grepl("tumor",Var1))
ggplot(prop_tb_group_1, 
       aes(x = Var1, y = prop,
           stratum = Var2, 
           alluvium = Var2, 
           fill = Var2, 
           label = Var2)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(alpha = 0.8,lwd=0.3) +
  scale_fill_manual(values = c("#bdc3d2","#88c4e8","#0074b3","#aeb0d8",
                               "#e5ce81","#f47720","#f7bfa4","#e6966a",
                               
                               "#aedacb","#4eb69e","#7c6aa4"))+
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



