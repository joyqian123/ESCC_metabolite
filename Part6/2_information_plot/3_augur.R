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
library(Augur)
library(miloR)
library(SingleCellExperiment)
library(SingleCellExperiment)
library(scater)
# library(scran)
library(dplyr)
library(patchwork)

# source("Resource/function.R")
# plan("multiprocess", workers = 10)
# plan()
setwd("~/help_for_others/DYQ/scRNA/")
options(future.globals.maxSize = 400 * 1024^3)
# rm(list=ls())


compare_ls = list(c("PP_tumor","SPP_tumor"),c("PP_tumor","IPP_tumor"),c("PP_pbmc","SPP_pbmc"),c("PP_pbmc","IPP_pbmc"))


# table(scRNA$response)


i=1


i=4
for(i in 1:2){
  
  sel_group = compare_ls[[i]]
  
  scRNA = qs::qread(paste0("./9_merge/augur/",sel_group[1],"_",sel_group[2],"/scRNA.qs"))
  augur = calculate_auc(scRNA,
                        cell_type_col = "Main_Celltype",
                        label_col = "group",
                        n_threads = 8)

  saveRDS(augur,file = paste0("./9_merge/augur/",sel_group[1],"_",sel_group[2],"/augur.rds"))
  
  
  
  
  
  ### 可视化
  augur = readRDS( paste0("./9_merge/augur/",sel_group[1],"_",sel_group[2],"/augur.rds"))
  data.plot <- subset(augur$results, metric == "roc_auc")
  order_1 = data.plot$cell_type %>% unique %>% as.character() %>% sort()
  data.plot$cell_type <- factor(data.plot$cell_type, levels = rev(augur$AUC$cell_type))
  
  col_order = match(levels(data.plot$cell_type),order_1)
  
  pdf(paste0("./9_merge/augur/",sel_group[1],"_",sel_group[2],"/augur.pdf"),width = 4,height = 5)
  print(
    ggplot(data.plot, aes(cell_type, estimate, color = cell_type)) +
      geom_violin() +
      geom_boxplot(width=0.3)+
      geom_jitter(size = 0.5) +
      geom_hline(yintercept = 0.5,lty=2)+
      scale_color_manual(values = c("#bdc3d2","#88c4e8","#0074b3","#aeb0d8",
                                    "#e5ce81","#f47720","#f7bfa4","#e6966a",
                                    
                                    "#aedacb","#4eb69e","#7c6aa4")[col_order])+
      coord_flip() +
      theme_bw()+
      theme(legend.position = "none")
  )
  dev.off()
  
  
  
  
  data2 <- subset(augur$results, metric == "roc_auc")
  data2 <- data2 %>% group_by(cell_type) %>% summarise(AUC = mean(estimate),yerr=sd(estimate))
  colnames(data2) <- c("Main_Celltype", "auc", "yerr")
  
  col_dt = data.frame(Main_Celltype=c("B","CD4T","CD8T","NK",
                                      "DC","mast","Mono_Macro","Neutrophil",
                                      "Fibroblast",
                                      "Endothelium","Malignant"),
                      color = c("#bdc3d2","#88c4e8","#0074b3","#aeb0d8",
                                "#e5ce81","#f47720","#f7bfa4","#e6966a",
                                "#aedacb","#4eb69e","#7c6aa4"))
  data2 = data2 %>% inner_join(col_dt,by="Main_Celltype") %>% arrange(auc) %>% 
    dplyr::mutate(Main_Celltype=factor(Main_Celltype,levels = Main_Celltype,ordered = T))
  
  pdf(paste0("./9_merge/augur/",sel_group[1],"_",sel_group[2],"/augur_erroplot.pdf"),width = 5,height = 5)
  print(
    ggplot(data2, aes(x = auc, y = Main_Celltype, color = color, label = color)) +
      geom_point(aes(size = auc)) +
      geom_errorbar(aes(xmin = auc - yerr, xmax = auc + yerr), width = 0.2,alpha=1) +
      scale_color_identity()+

      # geom_errorbarh(aes(xmin = logFC- xerr, xmax = logFC + xerr), height = 0.1,alpha=0.5) +
      geom_vline(xintercept = 0.6, linetype="dashed", color='black') +
      # geom_vline(xintercept = c(-1,1), linetype="dashed", color='black') +
      # ggrepel::geom_text_repel(show.legend = F,color="black") +
      theme_bw(base_size = 15)+
      ggtitle(paste0(sel_group[1]," vs ",sel_group[2]))+
      labs(x = "AUC score by augur") +
      theme(legend.position = "none")
  )
  dev.off()
  

  
  
  library(dplyr)
  da_results <- readRDS(file = paste0("./9_merge/miloR/",sel_group[1],"_",sel_group[2],"/data_da_results.rds"))
  data1 <- da_results %>% group_by(Main_Celltype) %>% summarise(logFC_mean = mean(logFC),xerr=sd(logFC))
  colnames(data1) <- c("Main_Celltype", "logFC", "xerr")

  
  # augur.res <- readRDS(file = "./11.3_augur/augur.rds")
  data2 <- subset(augur$results, metric == "roc_auc")
  data2 <- data2 %>% group_by(cell_type) %>% summarise(AUC = mean(estimate),yerr=sd(estimate))
  colnames(data2) <- c("Main_Celltype", "auc", "yerr")
  
  data_plot <- inner_join(data1, data2, by = "Main_Celltype")
  
  
  pdf(paste0("./9_merge/augur/",sel_group[1],"_",sel_group[2],"/combine_milo_augur_plot_type.pdf"),width = 5,height = 5)
  print(
    ggplot(data_plot, aes(x = logFC, y = auc, color = Main_Celltype, label = Main_Celltype)) +
      geom_point(size = 4) +
      geom_errorbar(aes(ymin = auc - yerr, ymax = auc + yerr), width = 0.1,alpha=1) +
      scale_color_manual(values = c("#bdc3d2","#88c4e8","#0074b3","#aeb0d8",
                                    "#e5ce81","#f47720","#f7bfa4","#e6966a",
                                    
                                    "#aedacb","#4eb69e","#7c6aa4"))+
      # geom_errorbarh(aes(xmin = logFC- xerr, xmax = logFC + xerr), height = 0.1,alpha=0.5) +
      geom_hline(yintercept = 0.5, linetype="dashed", color='black') +
      # geom_vline(xintercept = c(-1,1), linetype="dashed", color='black') +
      ggrepel::geom_text_repel(show.legend = F,color="black") +
      theme_bw(base_size = 15)+
      labs(x = "logFC by miloR",y = "AUC score by augur") +
      theme(legend.position = "none")
  )
  dev.off()
  
}
