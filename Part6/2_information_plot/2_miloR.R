######################################################################################
rm(list = ls())

run=TRUE

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
options(future.globals.maxSize = 400 * 1024^3)
# rm(list=ls())
scRNA <- qs::qread("./9_merge/scRNA_final_annotated.qs")
table(scRNA$orig.ident)
table(scRNA$Sub_Celltype)
table(scRNA$subtype)
table(scRNA$group)
table(scRNA$sample)

meta = scRNA@meta.data


compare_ls = list(c("PP_tumor","SPP_tumor"),c("PP_tumor","IPP_tumor"),c("PP_pbmc","SPP_pbmc"),c("PP_pbmc","IPP_pbmc"))
scRNA_all = scRNA

i=1
for (i in 1:length(compare_ls)){
  
  sel_group = compare_ls[[i]]
  scRNA = subset(scRNA_all,group %in% sel_group)
  
  if(run){
    dir.create(paste0("./9_merge/miloR/",sel_group[1],"_",sel_group[2]))
    ##########################################################################
    ###########4miloR
    # BiocManager::install("miloR")
    library(miloR)
    library(SingleCellExperiment)
    library(SingleCellExperiment)
    library(scater)
    # library(scran)
    library(dplyr)
    library(patchwork)
    # library(reticulate)
    # 增设行名
    meta.data <- scRNA@meta.data
    pca.embeddings <- Embeddings(scRNA, reduction = "harmony")
    umap.embeddings <- Embeddings(scRNA, reduction = "umap")
    sce <- SingleCellExperiment(assays=list(counts=scRNA@assays$RNA$counts),
                                reducedDims=SimpleList(PCA=pca.embeddings,  
                                                       UMAP=umap.embeddings))
    
    milo <- Milo(sce)
    colData(milo) <- DataFrame(meta.data)
    #### 2.Construct KNN graph
    reducedDimNames(milo)
    milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA")
    
    #### 3.Defining representative neighbourhoods on the KNN graph
    milo <- makeNhoods(milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "PCA")
    
    if(i <=2){
      #### 4.Counting cells in neighbourhoods
      milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="orig.ident")
      
      #### 5.Defining experimental design
      milo_design <- as.data.frame(colData(milo))
      milo_design <- milo_design %>% dplyr::select(group,orig.ident) %>% distinct(.,.keep_all = T)
      milo_design$group = factor(milo_design$group,levels = sel_group)
      rownames(milo_design) <- milo_design$orig.ident
      milo_design
    }else{
      #### 4.Counting cells in neighbourhoods
      milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="sample")
      
      #### 5.Defining experimental design
      milo_design <- as.data.frame(colData(milo))
      
      milo_design <- milo_design %>% dplyr::select(group,sample) %>% distinct(.,.keep_all = T) 
      milo_design$group = factor(milo_design$group,levels = sel_group)
      rownames(milo_design) <- milo_design$sample
      milo_design
    }
    
    #### 6.Computing neighbourhood connectivity
    #### 需要时间
    milo <- calcNhoodDistance(milo, d=30, reduced.dim = "PCA")
    # dir.create("./11.2_miloR")
    saveRDS(milo,file = paste0("./9_merge/miloR/",sel_group[1],"_",sel_group[2],"/milo_spatial_fdr.rds"))
  
  }
  milo = readRDS( paste0("./9_merge/miloR/",sel_group[1],"_",sel_group[2],"/milo_spatial_fdr.rds"))
  
  
  
  if(run){
    #### 7.差异分析
    da_results <- testNhoods(milo, design = ~ group, design.df = milo_design, reduced.dim = "PCA")
    
    da_results <- annotateNhoods(milo, da_results, coldata_col = "Main_Celltype")
    head(da_results)
    saveRDS(da_results,file = paste0("./9_merge/miloR/",sel_group[1],"_",sel_group[2],"/data_da_results.rds"))
  }
  da_results = readRDS(paste0("./9_merge/miloR/",sel_group[1],"_",sel_group[2],"/data_da_results.rds"))
  
  
  ### 标记细胞占比
  milo <- buildNhoodGraph(milo)
  
  pdf(paste0("./9_merge/miloR/",sel_group[1],"_",sel_group[2],"/nodeplot.pdf"),width = 6,height = 5)
  print(plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=1,res_column = "logFC")) 
  dev.off()

  order=da_results %>% group_by(Main_Celltype) %>% 
    dplyr::summarise(mean=mean(logFC)) %>% 
    arrange(desc(mean))
  da_results$Main_Celltype = factor(da_results$Main_Celltype,levels = rev(order$Main_Celltype),ordered = T)
  
  pdf(paste0("./9_merge/miloR/",sel_group[1],"_",sel_group[2],"/milo_plot.pdf"),width = 7,height = 7)
  print(
    plotDAbeeswarm(da_results, group.by = "Main_Celltype",alpha = 1)+
      geom_hline(yintercept=0, linetype="dashed")+
      scale_color_gradient2(midpoint=0, low="#1d3f7b", mid="white",high="#95295f", space ="Lab")
  )
  dev.off()
  
}
