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
dir.create("15_NK/")
scRNA<-
  qs::qread("4_Tcell/NK_SCT_AnnotatedWith0.8.qs")
# table(scRNA$Main_Celltype)
# 
# scRNA = subset(scRNA,Main_Celltype %in% c("CD4T","NKT"))


# 
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures() %>%
  ScaleData(features = rownames(scRNA))

##Regress MT value and Cell Cycle----
# scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^R[bp][sl]")

# scRNA <- SCTransform(scRNA, return.only.var.genes = F,
#                            vars.to.regress = c("percent.mt",
#                                                "percent.rb"))
scRNA <- SCTransform(scRNA, return.only.var.genes = F,
                           vars.to.regress = c("nFeature_RNA",
                                               "nCount_RNA","S.Score", "G2M.Score"))

##Remove mt Gene

mt.genes <- grep("^mt-", rownames(scRNA), value=T, ignore.case=T)
rb.genes <- grep("^R[p][sl]", rownames(scRNA), value=T, ignore.case=T)
hsp.genes <- grep("^Hsp", rownames(scRNA), value=T, ignore.case=T)
library(readxl)
rm_gene = read_excel("~/CHP1_scRNA_merge/reference/Resource/remove gene.xlsx")
rm_gene_tb = rm_gene %>% unlist() %>% na.omit()
h_2_m = read.table("~/CHP1_scRNA_merge/reference/Resource/m2h.txt")
rm_gene_mouse = h_2_m %>% dplyr::filter(V2 %in% rm_gene_tb) %>% pull(V1)

FeatureGene<-
  rownames(scRNA)%>%
  .[!.%in%c(rm_gene_mouse,mt.genes,rb.genes,hsp.genes)]



hvg <- VariableFeatures(scRNA,nfeatures = 500)
scRNA <- RunPCA(scRNA, verbose = F,
                      features = c(intersect(FeatureGene,hvg)),      ####把所有Feature gene替换为hvg  ###markers?
                      assay = "SCT")

scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT")

ElbowPlot(scRNA, ndims = 50)
pc.num=1:20
scRNA <- scRNA %>% RunTSNE(dims=pc.num,reduction = "harmony") %>%
  RunUMAP(dims=pc.num, reduction = "harmony") %>%
  FindNeighbors(dims = pc.num,reduction = "harmony") %>% FindClusters(resolution=0.2)



DimPlot(scRNA,group.by = "Main_Celltype",reduction = "umap",split.by = "group")  
DimPlot(scRNA,reduction = "umap")  



diff.wilcox = FindAllMarkers(scRNA,group.by = "SCT_snn_res.0.2",features = FeatureGene)


library(tidyverse)
library(UCell)
scRNA = AddModuleScore_UCell(scRNA,features = list(cytotoxic,inflammatory,stress),name = c("cytotoxic","inflammatory","stress"))


VlnPlot(scRNA,features = c("signature_1cytotoxic","signature_2inflammatory","signature_3stress"),group.by = "subtype")



VlnPlot(scRNA,features = c("signature_1cytotoxic","signature_2inflammatory","signature_3stress"),group.by = "group")
tb = FetchData(scRNA,vars = c("subtype","signature_1cytotoxic","signature_2inflammatory","signature_3stress")) %>%
  pivot_longer(cols = !subtype,names_to = "sig",values_to = "score")


qs::qsave(scRNA,file = "./15_NK/scRNA_final_NK.qs")


scRNA = qs::qread("./15_NK/scRNA_final_NK.qs")
median_tb <- tb %>%
  group_by(subtype,sig) %>%
  dplyr::summarise(median_expr = median(score, na.rm = TRUE)) %>%
  ungroup()

compare_ls = combn(unique(scRNA$subtype),2) 
compare_list <- list()
for(i in 1:ncol(compare_ls)) {
  compare_list[[i]] <- as.character(compare_ls[, i])
}

library(ggpubr)
pdf("./15_NK/NK_function_score.pdf",width = 13,height = 6)
ggplot(tb,aes(x=subtype,y=score,col=subtype,fill=subtype))+
  geom_violin(alpha=0.5)+
  stat_compare_means(comparisons = compare_list)+
  scale_color_manual(values = c("#982b2d","#de5a69","#ee9d9f","#fccfc9"))+
  scale_fill_manual(values = c("#982b2d","#de5a69","#ee9d9f","#fccfc9"))+
  geom_segment(data = median_tb,  # 使用提前计算的中位数数据
               aes(x= as.numeric(as.factor(subtype)) - 0.1,  # 调整x的值来控制线段的起点
                   xend = as.numeric(as.factor(subtype)) + 0.1,  # 调整xend的值来控制线段的终点
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
  facet_wrap(~sig,nrow=1,scale="free_y")+
  RotatedAxis()
dev.off()



DimPlot(scRNA,group.by = "subtype",split.by = "group")






table(scRNA$subtype)

NK_gene=c("Cdk11b","Cd27","Cx3cr1","Ccl3","Klrc2","Nfkbia","Ngf","Plekhh2","Sirpb1a","Tlr7","Crem","Rgs1","Nr4a3",
          "Il7r","Gzmh","Dnajb1","Mki67")
top10 <- top_n(group_by(diff.wilcox,cluster),n=5,wt = avg_log2FC)

pdf("./15_NK/marker_main_cluster.pdf",width = 5,height = 6)
DotPlot(scRNA,features = unique(c(NK_gene)),group.by = "subtype")+
  scale_color_gradient2(low = "white",mid = "#ADD8E6",high = "#6A5ACD")+
  coord_flip()+
  RotatedAxis()
dev.off()


pdf("./15_NK/Dimplot_subtype.pdf",width = 5.3,height = 3.2)
DimPlot(scRNA,reduction = "umap",group.by = "subtype",pt.size = 0.3,label = T,label.size = 2)+
  scale_color_manual(values = c("#982b2d","#de5a69","#ee9d9f","#fccfc9"))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()
pdf("./15_NK/Dimplot_subtype_split.pdf",width = 9,height = 5.5)
scRNA$group = factor(scRNA$group,levels = c("PP_tumor","SPP_tumor","IPP_tumor","PP_pbmc","SPP_pbmc","IPP_pbmc"),ordered = T)
DimPlot(scRNA,reduction = "umap",group.by = "subtype",pt.size = 0.15,split.by = "group",ncol=3)+
  scale_color_manual(values = c("#982b2d","#de5a69","#ee9d9f","#fccfc9"))+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()



NK = scRNA@meta.data
prop_tb = table(NK$sample,NK$subtype) %>% data.frame() %>% 
  group_by(Var1) %>% 
  dplyr::mutate(total=sum(Freq)) %>% dplyr::mutate(prop=Freq/total)

prop_tb$group = str_extract(prop_tb$Var1,"^[A-Z]+")
prop_tb$origin = str_extract(prop_tb$Var1,"[a-z]+$")
prop_tb$group = factor(prop_tb$group,levels = c("PP","SPP","IPP"),ordered = T)
pdf("./15_NK/sample_prop_change_NK_subtype.pdf",width = 8,height = 6)
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

pdf("./15_NK/all_prop_change.pdf",width = 8,height = 6)
ggplot(prop_tb_group, 
       aes(x = Var1, y = prop,
           stratum = Var2, 
           alluvium = Var2, 
           fill = Var2, 
           label = Var2)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(alpha = 0.8,lwd=0.3) +
  scale_fill_manual(values = c("#982b2d","#de5a69","#ee9d9f","#fccfc9"))+
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




prop_tb_group = prop_tb_group %>% dplyr::filter(Var1 %in% c("PP_tumor","SPP_tumor","IPP_tumor"))
pdf("./15_NK/all_prop_change_only_tumor.pdf",width = 6,height = 4)
ggplot(prop_tb_group, 
       aes(x = Var1, y = prop,
           stratum = Var2, 
           alluvium = Var2, 
           fill = Var2, 
           label = Var2)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(alpha = 0.8,lwd=0.3) +
  scale_fill_manual(values = c("#982b2d","#de5a69","#ee9d9f","#fccfc9"))+
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
pdf("./15_NK/all_prop_change_only_pbmc.pdf",width = 6,height = 4)
ggplot(prop_tb_group, 
       aes(x = Var1, y = prop,
           stratum = Var2, 
           alluvium = Var2, 
           fill = Var2, 
           label = Var2)) +
  geom_flow(alpha = 0.5) +
  geom_stratum(alpha = 0.8,lwd=0.3) +
  scale_fill_manual(values = c("#982b2d","#de5a69","#ee9d9f","#fccfc9"))+
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









