rm(list = ls())
# .libPaths(c("/home/rstudio_cu02/anaconda3/envs/R4.2/lib/R/library","/home/qianzy/R/4.0/"))
library(Seurat)
library("sf",lib.loc ="/home/rstudio_cu02/anaconda3/envs/R4.0/lib/R/library")
library(monocle3)
# find.package("monocle3")

# The tutorial shown below and on subsequent pages uses two additional packages:
library(ggplot2)
library(dplyr)
library(remotes)
# memory.limit()
# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)
library(dplyr)
library(patchwork)
library(ggplot2)
library(future)
setwd("~/help_for_others/DYQ/scRNA//")
options(future.globals.maxSize = 400 * 1024^3)

library(harmony)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(copykat)
library(Rtsne)
library(scCustomize)
library(igraph)
library(monocle)

scRNA<-
  qs::qread("4_Tcell/scRNA_CD8Tcell_SCT_AnnotatedWith0.8.qs")



######################################################
######################################################
expr_matrix <- as(as.matrix(scRNA@assays$RNA$counts),"sparseMatrix")
dim(expr_matrix)
p_data <- scRNA@meta.data
dim(p_data)
# p_data$Sub_Celltype <- scRNA_T@active.ident
f_data <- data.frame(gene_short_name=row.names(scRNA@assays$RNA$counts),
                     row.names = row.names(scRNA@assays$RNA$counts))
dim(f_data)

pd <- new("AnnotatedDataFrame", data = p_data)
fd <- new("AnnotatedDataFrame", data = f_data)
fd <- as(fd, "data.frame")
pd <- as(pd, "data.frame")

cds <- new_cell_data_set(expr_matrix, cell_metadata = pd, gene_metadata  = fd)

cds <- preprocess_cds(cds, num_dim = 50)
monocle3::plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,reduction_method="UMAP")
colData(cds)
plot_cells(cds, reduction_method="UMAP", color_cells_by="subtype")



cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="subtype",group_label_size = 4,cell_size = 1) + 
  ggtitle('int.umap')


cds <- cluster_cells(cds, resolution=0.005)
plot_cells(cds,reduction_method = "UMAP",color_cells_by = "subtype",group_label_size = 4,cell_size = 1)
# plot_cells(cds,reduction_method = "UMAP",color_cells_by = "partition",group_label_size = 4,cell_size = 1)
# plot_cells(cds,reduction_method = "UMAP",color_cells_by = "cluster",group_label_size = 4,cell_size = 1)

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE)+
  ggtitle("label by partitionID")


## 识别轨迹
# cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, color_cells_by = "subtype",
           label_branch_points = FALSE,group_label_size = 4,cell_size = 1)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, color_cells_by = "cluster",
           label_branch_points = FALSE,group_label_size = 4,cell_size = 1)
p=plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, color_cells_by = "cluster",
             label_branch_points = FALSE,group_label_size = 4,cell_size = 1)


##细胞按拟时排序
cds <- order_cells(cds) #存在bug，使用辅助线选择root细胞

cds <- order_cells(cds)


library(scico)
library(rcolors)
library(RColorBrewer)
# gradientColors <- brewer.pal("Hokusai3", n=6, type="continuous")
pdf("./11_Tcell/monocle/1.monocle_split.pdf",width = 23,height = 4)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE,  label_branch_points = FALSE,group_label_size = 4,cell_size = 1) +
  scale_color_scico(palette = 'vik', na.value = "transparent")+
  facet_grid(~group)
dev.off()
pdf("./11_Tcell/monocle/1.monocle.pdf",width = 5,height = 4)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE,  label_branch_points = FALSE,group_label_size = 4,cell_size = 1) +
  scale_color_scico(palette = 'vik', na.value = "transparent")#+
  # facet_grid(~group)
dev.off()

save_monocle_objects(cds=cds, directory_path='./11_Tcell/monocle/monocle_object')

cds <- load_monocle_objects(directory_path='./11_Tcell/monocle/monocle_object')


p=plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
             label_leaves = FALSE,  label_branch_points = FALSE,group_label_size = 4,cell_size = 1)



pseudotime = p$data$cell_color
dysfunction_gene = c("Havcr2","Rgs16","Layn","Srgap3","Dusp4","Csf1","Tnfrsf9",
                       "Lyst","Tnfrsf18","Ndfip2","Sqle","Id3","Sox4","Cd9",
                       "Phlda1","Ccl3","Ccl4","Klrc1","Klrd1","Klrb1","Cdk6",
                       "Pls3","Afap1l2","Ctsw","Il2ra","Ahl1","Rbpj",
                       "Gzmb","Gnly","Pdcd1","Tigit","Lag3")


dim(scRNA)
scRNA$pseudotime = pseudotime

library(UCell)
scRNA = AddModuleScore_UCell(scRNA,features = list(dysfunction_gene),name="dysfunction_score")

qs::qsave(scRNA,"./11_Tcell/scRNA_add_pseudotime_and_dysfunction_score.qs")






scRNA = qs::qread("./11_Tcell/scRNA_add_pseudotime_and_dysfunction_score.qs")
FeaturePlot(scRNA,features = "signature_1dysfunction_score")+
  scale_color_gradient2(low = "white",mid = "#f8d6e5",high = "#626e98",midpoint = 0.1)


colData(cds)$dysfunction_score = scRNA$signature_1dysfunction_score

pdf("./11_Tcell/monocle/1.monocle_dysfunction_score.pdf",width = 5,height = 4)
plot_cells(cds, color_cells_by = "dysfunction_score", label_cell_groups = FALSE,
           label_leaves = FALSE,  label_branch_points = FALSE,group_label_size = 4,cell_size = 1) +
  scale_color_gradient2(low = "white",mid = "#f8d6e5",high = "#626e98",midpoint = 0.05)+
  ggtitle("Dysfunction_score")
# facet_grid(~group)
dev.off()
pdf("./11_Tcell/monocle/1.monocle_dysfunction_score_split.pdf",width = 23,height = 4)
plot_cells(cds, color_cells_by = "dysfunction_score", label_cell_groups = FALSE,
           label_leaves = FALSE,  label_branch_points = FALSE,group_label_size = 4,cell_size = 1) +
  scale_color_gradient2(low = "white",mid = "#f8d6e5",high = "#626e98",midpoint = 0.05)+
  facet_grid(~group)
dev.off()


VlnPlot(scRNA,features = "signature_1dysfunction_score",group.by = "group")

scRNA_tumor = subset(scRNA,origin=="tumor")
tb = FetchData(scRNA_tumor,vars = c("group","signature_1dysfunction_score"))


median_tb <- tb %>%
  group_by(group) %>%
  dplyr::summarise(median_expr = median(signature_1dysfunction_score, na.rm = TRUE)) %>%
  ungroup()

library(ggpubr)
pdf("./11_Tcell/monocle/dysfunction_score.pdf",width = 4,height = 4)
ggplot(tb,aes(x=group,y=signature_1dysfunction_score,col=group,fill=group))+
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
  RotatedAxis()
dev.off()


tb = FetchData(scRNA_tumor,vars = c("group","pseudotime"))


median_tb <- tb %>%
  group_by(group) %>%
  dplyr::summarise(median_expr = median(pseudotime, na.rm = TRUE)) %>%
  ungroup()

library(ggpubr)
pdf("./11_Tcell/monocle/pseudotime.pdf",width = 4,height = 4)
ggplot(tb,aes(x=group,y=pseudotime,col=group,fill=group))+
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
  RotatedAxis()
dev.off()




