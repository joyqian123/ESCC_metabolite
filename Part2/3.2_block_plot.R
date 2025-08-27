#!/usr/bin/Rscript
rm(list = ls())

setwd("~/help_for_others/DYQ/output/3_baseline_model/mt_block//")
load("./final_clone_0.8.rdata")

# 初始化相关系数矩阵为单位矩阵
n <- length(unique(final$mt))
cor_matrix <- diag(n)

# 遍历每一对代谢物
for (i in 1:n) {
  for (j in 1:n) {
    if (i == j) {
      cor_matrix[i, j] <- 1
    } else {
      # 获取第i个和第j个代谢物的名称
      metabolite_i <- unique(final$mt)[i]
      metabolite_j <- unique(final$mt)[j]
      
      # 获取它们的组别
      group_i <- final$clone.id[final$mt == metabolite_i]
      group_j <- final$clone.id[final$mt == metabolite_j]
      
      # 如果同一组，相关系数设为1
      if (group_i == group_j) {
        cor_matrix[i, j] <- 1
      } else {
        cor_matrix[i, j] <- 0
      }
    }
  }
}

# 查看相关系数矩阵
print(cor_matrix)
colnames(cor_matrix)=final$mt
rownames(cor_matrix)=final$mt

cor_matrix[upper.tri(cor_matrix,diag = TRUE)]=NA

check.correlation = cor_matrix  %>% data.frame(.,check.names = FALSE) %>% 
  rownames_to_column("A") %>% pivot_longer(cols = !A,names_to = "B",values_to = "cor") %>% 
  dplyr::filter(cor>0)
check.correlation$B = str_remove(check.correlation$B,"X")

load("~/help_for_others/DYQ/output/1_preclean/data/data_C1_preclean.rdata")
mt_info = data_in_C1 %>% dplyr::select(1:3) %>% dplyr::rename(mt=Compounds)

check.correlation = check.correlation %>% dplyr::filter(A %in% mt_info$mt,B %in% mt_info$mt)

library(igraph)
library(data.table)
graph.cor <- graph_from_data_frame(check.correlation, directed = FALSE)


# save(mt_info,file = "~/help_for_others/DYQ/data/mt_info.rdata")

graph_mt_ls =  names(V(graph.cor)) %>% data.frame()
colnames(graph_mt_ls)[1]="mt"
graph_mt_ls$mt <- factor(graph_mt_ls$mt,levels = graph_mt_ls$mt,ordered = T)



library(RColorBrewer)
random_col=c(brewer.pal(8,"Dark2"),brewer.pal(12,"Paired")[c(1:2)]) %>% rev()

col_distribution = data.frame(clone.id=unique(mt_info$`Class I`),col=rep(random_col,2)[-20],shape=c(rep("circle",10),rep("square",9)))
graph_mt_ls_col = mt_info %>% dplyr::filter(mt %in% graph_mt_ls$mt) %>% dplyr::rename(clone.id=`Class I`) %>% 
  inner_join(col_distribution,by = "clone.id") %>% dplyr::select(mt,col,shape)
colnames(graph_mt_ls_col)[1]="mt"
graph_mt_ls_col$mt <- factor(graph_mt_ls_col$mt,levels = graph_mt_ls$mt,ordered = T)
graph_mt_ls_col <- graph_mt_ls_col %>% arrange(mt)


V(graph.cor)$color <- graph_mt_ls_col$col
V(graph.cor)$shape <- graph_mt_ls_col$shape
E(graph.cor)$width <- 0.05
E(graph.cor)$opacity <- 0.05
E(graph.cor)$color <- "#E0EEEE"
edge_attr(graph.cor)

length(V(graph.cor)) # 查看节点数
length(E(graph.cor)) # 查看边数

set.seed(42)
coords <- layout_with_fr(graph.cor, niter=9999,grid="nogrid")
pdf("./graph_class_I_new.pdf",width=20,height=15)
plot(graph.cor,vertex.label="",vertex.size=1,vertex.frame.color = "transparent",
     layout=coords)
legend("bottomright",legend = col_distribution$clone.id,col = col_distribution$col,
       pch=c(rep(16,10),rep(15,10)))
dev.off()
