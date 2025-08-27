rm(list = ls())
library(CellChat)
library(Seurat)
library(Matrix)
library(dplyr)
library(hdf5r)
library(tidyverse)
library(stringr)

setwd("~/help_for_others/DYQ/scRNA//")


scRNA = qs::qread("./9_merge/scRNA_final_annotated.qs")
meta = scRNA@meta.data

scRNA_NK = qs::qread("./15_NK/scRNA_final_NK.qs")
meta_NK = scRNA_NK@meta.data %>% dplyr::select(cellid,subtype)

meta_NK_add = meta %>% dplyr::filter(Main_Celltype=="NK") %>% dplyr::select(-subtype) %>% 
  inner_join(meta_NK,by = "cellid") %>% 
  dplyr::mutate(type_for_cellchat = subtype)

meta_not_NK = meta %>% dplyr::filter(!Main_Celltype=="NK") %>% 
  dplyr::mutate(type_for_cellchat = Main_Celltype)

meta_final = rbind(meta_NK_add,meta_not_NK)
rownames(meta_final)=meta_final$cellid

table(meta_final$type_for_cellchat)

scRNA = AddMetaData(scRNA,metadata = meta_final)
table(scRNA$type_for_cellchat)

scRNA = subset(scRNA,group %in% c("PP_tumor","IPP_tumor"))


co <- c("PP_tumor","IPP_tumor")

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
unique(CellChatDB$interaction$pathway_name)
cellchat_list <- list()
i=1
for (i in 1:2){
  object <- subset(scRNA,group==co[i])
  # meta <- object@assays$RNA$counts
  # idents = data.frame(row.names=rownames(object@meta.data), celltype=object$Sub_Celltype)
  Idents(object)=object$type_for_cellchat
  # head(idents)
  cellchat <- createCellChat(object, group.by='ident',assay = "RNA")
  
  # CellChatDB.ss <- subsetDB(CellChatDB, search = "Secreted Signaling", key='annotation')
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  # 识别在单个细胞类型中过表达配/受体
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  # 识别过表达互作对
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # 平滑表达值（目的是消除dropout影响，可选不用）
  # We also provide a function to project gene expression data onto protein-protein interaction (PPI) network. Specifically, a diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network. This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors. 
  # cellchat <- projectData(cellchat, PPI.human)
  
  
  # 互作可能性计算
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)    
  
  # 过滤表达细胞比例低的互作对
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # # 提取关注得细胞间通讯关系
  # df.net <- subsetCommunication(cellchat)
  # head(df.net)
  # 
  # # 提取特定细胞类型间通讯关系
  # df.net2 <- subsetCommunication(cellchat, sources.use=c('Naive CD4 T','Memory CD4 T','CD8 T') , targets.use= c('CD14+ Mono','B'))
  # head(df.net2)
  # 
  # # 提取特定信号通路中的通讯关系
  # df.net3 <- subsetCommunication(cellchat, signaling = c("IL16", "ANNEXIN"))
  # head(df.net3)
  
  
  # NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat@net
  
  
  # 整合通讯网络结果
  # USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
  cellchat <- aggregateNet(cellchat)
  
  cellchat_list[[i]]<- cellchat
  names(cellchat_list)[i]<-co[i]
}

dir.create("./15_NK/cellchat/IPP")
save(cellchat_list,file = "./15_NK/cellchat/IPP/cellchat_list.rdata")
