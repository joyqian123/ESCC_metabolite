#!/usr/bin/Rscript
rm(list = ls())
load("~/help_for_others/DYQ/output/3_baseline_model/csv/df_training_norm.rdata")

dir.create("~/help_for_others/DYQ/output/3_baseline_model/mt_block")
setwd("~/help_for_others/DYQ/output/3_baseline_model//")





########用全部C1、C6的数据，原始且不归一化的数据来做
cor_matrix = cor(training[,-c(1)])
mt_ls = colnames(cor_matrix)


flattenCorrMatrix <- function(cormat){
  ut <- upper.tri(cormat) 
  return(data.frame(row = rownames(cormat)[row(cormat)[ut]], 
                    column = rownames(cormat)[col(cormat)[ut]], cor =(cormat)[ut]))
}
#举个栗子
sparse_mat = flattenCorrMatrix(cor_matrix)

correlation.cutoff=0.8

check.correlation = sparse_mat %>% dplyr::filter(cor>correlation.cutoff)
colnames(check.correlation) <- c("row", "col","weight")

check.correlation <- as.matrix(check.correlation)

library(igraph)
library(data.table)
graph.cor <- graph_from_data_frame(check.correlation, directed = FALSE)


load("~/help_for_others/DYQ/output/1_preclean/data/data_C1_preclean.rdata")
mt_info = data_in_C1 %>% dplyr::select(1:3) %>% dplyr::rename(mt=Compounds)

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
pdf("./mt_block/graph_class_I.pdf",width=30,height=20)
plot(graph.cor,vertex.label="",vertex.size=1,vertex.frame.color = "transparent",
     layout=coords)
legend("bottomright",legend = col_distribution$clone.id,col = col_distribution$col,
       pch=c(rep(16,10),rep(15,10)))
dev.off()



#######################################################################
########################################################################
groups.cor <- split(unique(as.vector(check.correlation[,c(1,2)])),
                    clusters(graph.cor)$membership)

l <- seq(1, length(groups.cor))
df.conv <- apply(as.matrix(l), 1, function(x) {
  data.table(clone.id = x, mt = groups.cor[[x]])
}) %>% do.call(rbind,.)




clone_ls = df.conv$clone.id %>% unique()

res_each_simi <- lapply(1:length(clone_ls),function(tt){
  
  clone_mt_ls = df.conv %>% dplyr::filter(clone.id==clone_ls[tt]) %>% pull(mt)
  tmp = training %>% dplyr::select(any_of(clone_mt_ls))
  
  cor_tb = cor(tmp)
  
  return(data.frame(clone_No=clone_ls[tt],right_cor=table(cor_tb>=0.8)["TRUE"]-length(clone_mt_ls),
                    total=length(clone_mt_ls)*(length(clone_mt_ls)-1)))
  
}) %>% do.call(rbind,.)
res_each_simi$high_cor_prob = res_each_simi$right_cor/res_each_simi$total

not_1_clone_ls_need_split <- res_each_simi %>% dplyr::filter(high_cor_prob<1) %>% pull(clone_No)





##################################################
res_each_deep_split <- lapply(1:length(not_1_clone_ls_need_split),function(tt){
  print(tt)
  clone_mt_ls = df.conv  %>% dplyr::filter(clone.id==not_1_clone_ls_need_split[tt]) %>% pull(mt)
  
  cutoff_seq=seq(0.8,0.99,0.01)
  
  tmp = training %>% dplyr::select(any_of(clone_mt_ls))
  cor_tb = cor(tmp)
  
  for(ss in cutoff_seq){
    
    print(ss)
    
    correlation.cutoff=ss
    
    sparse_mat = flattenCorrMatrix(cor_tb)
    
    check.correlation = sparse_mat %>% dplyr::filter(cor>correlation.cutoff)
    colnames(check.correlation) <- c("row", "col","weight")
    
    if(nrow(check.correlation)==0){
      break
    }
    
    check.correlation <- as.matrix(check.correlation)
    graph.cor <- graph_from_data_frame(check.correlation, directed = FALSE)
    
    groups.cor <- split(unique(as.vector(check.correlation[,c(1,2)])),
                        clusters(graph.cor)$membership)
    
    l <- seq(1, length(groups.cor))
    df.conv <- apply(as.matrix(l), 1, function(x) {
      data.table(clone.id = x, mt = groups.cor[[x]])
    }) %>% do.call(rbind,.)
    
    
    
    df.conv$origin_clone_id = not_1_clone_ls_need_split[tt]
    
    colnames(df.conv)[1]=str_c("cutoff",ss,sep = "_")
    
    if(ss==0.8){
      tmp_2 = df.conv
    }  else{
      tmp_2 = tmp_2 %>% full_join(df.conv,by = c("mt","origin_clone_id"))
    }
  }
  
  tmp_2 = tmp_2 %>% dplyr::select(mt,origin_clone_id,everything())
  return(tmp_2)
}) 





########################################
unmatch_sub_clone <- lapply(1:length(res_each_deep_split),function(uu){
  
  print(uu)
  clone_tb = unique(res_each_deep_split[[uu]])
  
  pp=6
  res_each_simi <- lapply(3:ncol(clone_tb),function(pp){
    
    print(pp)
    clone_ls <- clone_tb %>% dplyr::select(1:2,any_of(pp)) %>% dplyr::rename(clone.id=colnames(.)[3]) %>% 
      na.omit()
    clone_No <- unique(unlist(clone_ls[,3]))
    
    tt=1
    res_each_simi_sub <- lapply(1:length(clone_No),function(tt){
      
      print(tt)
      clone_mt_ls = clone_ls  %>% dplyr::filter(clone.id==clone_No[tt]) %>% pull(mt)
      
      tmp = training %>% dplyr::select(any_of(clone_mt_ls))
      
      cor_tb = cor(tmp)
      
      return(data.frame(clone_split=colnames(clone_tb)[pp],clone_No=clone_No[tt],right_cor=table(cor_tb>=0.8)["TRUE"]-length(clone_mt_ls),
                        total=length(clone_mt_ls)*(length(clone_mt_ls)-1)))
      
      
    }) %>% do.call(rbind,.)
    
  }) %>% do.call(rbind,.)
  
  res_each_simi$high_cor_prob = res_each_simi$right_cor/res_each_simi$total  
  
  return(res_each_simi)
}) 





i=1

final_clone_deep_split = lapply(1:length(res_each_deep_split),function(i){
  
  tmp = res_each_deep_split[[i]]
  black_ls = unmatch_sub_clone[[i]]
  
  tmp$final_cutoff = sapply(1:nrow(tmp),function(rr){
    print(rr)
    tmp_line = tmp[rr,]
    for (qq in 3:ncol(tmp_line)){
      print(qq)
      exist_cutoff=colnames(tmp_line)[qq]
      exist_clone= tmp_line %>% dplyr::select(any_of(qq)) %>% as.numeric()
      if(is.na(exist_clone)){
        final_cutoff=exist_cutoff
        final_clone=NA
        break
      }
      is_in_black_ls = black_ls %>% dplyr::filter(clone_split==exist_cutoff,clone_No==exist_clone) %>% pull(high_cor_prob)
      if(is_in_black_ls<1){
        next
      }else if(is_in_black_ls==1){
        final_cutoff=exist_cutoff
        final_clone=exist_clone
        break
      }
    }
    if(exists("final_cutoff")){
      return(final_cutoff)
    }else{
      return(exist_cutoff)
    }
  }) 
  
  tmp$final_clone = sapply(1:nrow(tmp),function(rr){
    tmp_line = tmp[rr,]
    for (qq in 3:ncol(tmp_line)){
      exist_cutoff=colnames(tmp_line)[qq]
      exist_clone= tmp_line %>% dplyr::select(any_of(qq)) %>% as.numeric()
      if(is.na(exist_clone)){
        final_cutoff=exist_cutoff
        final_clone=NA
        break
      }
      is_in_black_ls = black_ls %>% dplyr::filter(clone_split==exist_cutoff,clone_No==exist_clone) %>% pull(high_cor_prob)
      if(is_in_black_ls<1){
        next
      }else if(is_in_black_ls==1){
        final_cutoff=exist_cutoff
        final_clone=exist_clone
        break
      }
    }
    if(exists("final_clone")){
      return(final_clone)
    }else{
      return(exist_clone)
    }
  }) 
  
  return(tmp)
})






tt=1
final_clone_unpair = lapply(1:length(final_clone_deep_split),function(tt){
  tmp = final_clone_deep_split[[tt]] %>% dplyr::select(mt,origin_clone_id,final_cutoff,final_clone) %>% 
    dplyr::mutate(final_clone_name=str_c(final_cutoff,final_clone,sep = "_"))
  origin = unique(tmp$origin_clone_id)
  tmp_clone_size = table(tmp$final_clone_name) %>% as.data.frame() %>% 
    dplyr::mutate(name=str_c(origin,"sub",1:nrow(.),sep = "_")) %>% 
    dplyr::rename(final_clone_name=Var1)
  tmp_match = tmp %>% dplyr::filter(!is.na(final_clone_name)) %>% full_join(tmp_clone_size,by = "final_clone_name") %>% 
    dplyr::select(mt,origin_clone_id,name) %>% dplyr::rename(clone_id=name)
  tmp_unmatch =  tmp %>% dplyr::filter(is.na(final_clone_name))
  if(nrow(tmp_unmatch)>0){
    tmp_unmatch = tmp_unmatch %>% dplyr::mutate(name=str_c(origin,"sub",nrow(tmp_clone_size)+(1:nrow(tmp_unmatch)),sep = "_")) %>% 
      dplyr::select(mt,origin_clone_id,name) %>% dplyr::rename(clone_id=name)
    tmp_final = rbind(tmp_match,tmp_unmatch)
  }else{
    tmp_final = tmp_match
  }
  return(tmp_final)
}) %>% do.call(rbind,.)


final_pair = df.conv %>% dplyr::filter(!clone.id %in% unique(final_clone_unpair$origin_clone_id))
final_clone_unpair = final_clone_unpair %>% dplyr::select(clone_id,mt) %>% dplyr::rename(clone.id=clone_id)

final=rbind(final_pair,final_clone_unpair)

all_unpair = data.frame(mt=setdiff(mt_ls,final$mt)) %>% dplyr::mutate(clone.id=seq(168,168+nrow(.)-1,1))

final = rbind(final,all_unpair)
length(unique(final$clone.id))


save(final,file = "./mt_block/final_clone_0.8.rdata")
