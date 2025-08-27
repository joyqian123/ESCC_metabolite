rm(list = ls())
setwd("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/")

load("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/VD/mt_list_from_data.rdata")

load("../1_preclean/data/data_C1_preclean.rdata")
mt_info = data_in_C1 %>% dplyr::select(Compounds,`Class I`,`Class II`) %>% 
  dplyr::rename(mt=Compounds)
class_ID_1 = mt_info %>% dplyr::select(mt,`Class I`) %>% dplyr::rename(class=colnames(.)[2])
class_ID_2 = mt_info %>% dplyr::select(mt,`Class II`) %>% dplyr::rename(class=colnames(.)[2])



kk=1
library(clusterProfiler)
en_res_1 = lapply(1:length(mt_list),function(kk){
  tmp = enricher(mt_list[[kk]],TERM2GENE = class_ID_2[,c(2,1)])
  tmp_1 = tmp@result %>% dplyr::mutate(class=names(mt_list)[kk])
  return(tmp_1)
}) %>% do.call(rbind,.)


en_res_1_plot = en_res_1 %>% dplyr::filter(pvalue<0.1) %>% dplyr::select(ID,GeneRatio,pvalue,class)
en_res_1_plot = separate(en_res_1_plot,col = "GeneRatio",into = c("A","B"),sep = "/") %>% 
  dplyr::mutate(A=as.numeric(A),B=as.numeric(B)) %>% 
  dplyr::mutate(r=A/B) %>% dplyr::select(ID,r,pvalue,class) 

# en_res_1_plot = en_res_1_plot %>% dplyr::filter(!ID=="TG")
r_tb = en_res_1_plot %>% dplyr::select(ID,r,class) %>%
  pivot_wider(id_cols = ID,names_from = class,values_from = r) %>% column_to_rownames("ID")
p_tb = en_res_1_plot %>% dplyr::select(ID,pvalue,class) %>%
  pivot_wider(id_cols = ID,names_from = class,values_from = pvalue) %>% column_to_rownames("ID")

library(ggplot2)
library(Seurat)
pdf("./plot/enricher_of_class_mt_type.pdf",width = 5.5,height = 5)
ggplot(en_res_1_plot)+
  geom_point(aes(x=class,y=ID,col=r,size=-log(pvalue)))+
  scale_color_gradient2(low = "#FFB6C1",mid = "#DDA0DD",high = "#C71585")+
  geom_vline(xintercept = seq(1.5,5.5,1),lwd=0.75,lty=3,col="lightgrey")+
  geom_hline(yintercept = seq(1.5,13.5,1),lwd=0.75,lty=3,col="lightgrey")+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()) +
  RotatedAxis()
dev.off()  



load("../../reference/kegg_ID_sum.rdata")
library(clusterProfiler)
en_res_2 = lapply(1:length(mt_list),function(kk){
  tmp = enricher(mt_list[[kk]],TERM2GENE = ID_sum[,c(2,1)])
  tmp_1 = tmp@result %>% dplyr::mutate(class=names(mt_list)[kk])
  return(tmp_1)
}) %>% do.call(rbind,.)

en_res_2_plot = en_res_2 %>% dplyr::filter(pvalue<0.1) %>% dplyr::select(ID,GeneRatio,pvalue,class,geneID) 
rm_pathway = en_res_2_plot %>% dplyr::select(ID,geneID) %>% 
  distinct(geneID,.keep_all = T) %>% pull(ID)
en_res_2_plot = en_res_2_plot %>% dplyr::filter(ID %in% rm_pathway)
en_res_2_plot = separate(en_res_2_plot,col = "GeneRatio",into = c("A","B"),sep = "/") %>% 
  dplyr::mutate(A=as.numeric(A),B=as.numeric(B)) %>% 
  dplyr::mutate(r=A/B) %>% dplyr::select(ID,r,pvalue,class) 

library(ggplot2)
library(Seurat)
pdf("./plot/enricher_of_class_mt_pathway.pdf",width =6.5,height = 7)
ggplot(en_res_2_plot)+
  geom_point(aes(x=class,y=ID,col=r,size=-log(pvalue)))+
  scale_color_gradient2(low = "#FFB6C1",mid = "#DDA0DD",high = "#C71585")+
  geom_vline(xintercept = seq(1.5,5.5,1),lwd=0.75,lty=3,col="lightgrey")+
  geom_hline(yintercept = seq(1.5,24.5,1),lwd=0.75,lty=3,col="lightgrey")+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()) +
  RotatedAxis()
dev.off()  






########################################################################################
#########################################################################################
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

TG = mt_info %>% dplyr::filter(`Class II`=="TG")

merge_TG = merge %>% dplyr::filter(mt %in% TG$mt)



###################################################################
###################################################################
continuous_R_TG = mt_list[[1]] %>% .[grep("TG",.)] %>% data.frame() %>% dplyr::rename(mt=colnames(.)[1]) %>% 
  dplyr::mutate(mt=str_remove(mt,"TG\\(")) %>% dplyr::mutate(mt=str_remove(mt,"\\)")) %>% 
  separate(.,col=mt,into=c("A","B","C"),sep="_")
continuous_R_sum = continuous_R_TG %>% pivot_longer(cols = 1:3,names_to = "sn",values_to = "type") %>% 
  dplyr::mutate(count=1) %>% 
  group_by(type) %>% dplyr::summarize(sum1=sum(count)) %>% 
  dplyr::mutate(prob=sum1/sum(sum1)) %>% 
  dplyr::mutate(class="continuous_R")

baseline_NR_TG = mt_list[[4]] %>% .[grep("TG",.)] %>% data.frame() %>% dplyr::rename(mt=colnames(.)[1]) %>% 
  dplyr::mutate(mt=str_remove(mt,"TG\\(")) %>% dplyr::mutate(mt=str_remove(mt,"\\)")) %>% 
  separate(.,col=mt,into=c("A","B","C"),sep="_")
baseline_NR_sum = baseline_NR_TG %>% pivot_longer(cols = 1:3,names_to = "sn",values_to = "type") %>% 
  dplyr::mutate(count=1) %>% 
  group_by(type) %>% dplyr::summarize(sum1=sum(count)) %>% 
  dplyr::mutate(prob=sum1/sum(sum1)) %>% 
  dplyr::mutate(class="baseline_NR")


response_R_TG = mt_list[[6]] %>% .[grep("TG",.)] %>% data.frame() %>% dplyr::rename(mt=colnames(.)[1]) %>% 
  dplyr::mutate(mt=str_remove(mt,"TG\\(")) %>% dplyr::mutate(mt=str_remove(mt,"\\)")) %>% 
  separate(.,col=mt,into=c("A","B","C"),sep="_")
response_R_sum = response_R_TG %>% pivot_longer(cols = 1:3,names_to = "sn",values_to = "type") %>% 
  dplyr::mutate(count=1) %>% 
  group_by(type) %>% dplyr::summarize(sum1=sum(count)) %>% 
  dplyr::mutate(prob=sum1/sum(sum1)) %>% 
  dplyr::mutate(class="response_R")



TG_final = rbind(continuous_R_sum,baseline_NR_sum,response_R_sum) %>% 
  dplyr::select(type,class,prob) %>% 
  pivot_wider(id_cols = c("type"),names_from = "class",values_from = "prob")

TG_final = TG_final %>% arrange(type) %>% column_to_rownames("type")
library(ComplexHeatmap)

col_fun = colorRamp2(breaks = c(0,0.1,0.3),colors = c("lightgrey","#DDA0DD","#C71585"))

pdf("./plot/TG_compare_1.pdf",width = 3,height = 7)
Heatmap(as.matrix(TG_final),na_col = "lightgrey",cluster_rows = FALSE,col = col_fun,
        cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd = 0.75),border = TRUE)
dev.off()





###################################################################
###################################################################
TG_ls = merge_TG$mt %>% .[grep("TG",.)] %>% data.frame() %>% dplyr::rename(mt=colnames(.)[1]) %>% 
  dplyr::mutate(mt1=str_remove(mt,"TG\\(")) %>% dplyr::mutate(mt1=str_remove(mt1,"\\)")) %>% 
  separate(.,col=mt1,into=c("A","B","C"),sep="_") %>% 
  separate(.,col=A,into=c("A1","A2"),remove = FALSE,sep=":")%>% 
  separate(.,col=B,into=c("B1","B2"),remove = FALSE,sep=":")%>% 
  separate(.,col=C,into=c("C1","C2"),remove = FALSE,sep=":") %>% 
  dplyr::mutate(A1=as.numeric(A1),A2=as.numeric(A2),B1=as.numeric(B1),B2=as.numeric(B2),C1=as.numeric(C1),C2=as.numeric(C2)) %>% 
  dplyr::mutate(total_C = A1+B1+C1,total_unsature=A2+B2+C2)

type_ls = c(TG_ls$A,TG_ls$B,TG_ls$C) %>% unique()

chain_type = lapply(1:length(type_ls),function(kk){
  tmp = TG_ls %>% dplyr::filter(A %in% type_ls[kk] | B %in% type_ls[kk] | C %in% type_ls[kk])
  res_tmp = merge_TG %>% dplyr::select(mt,FC_C1,p_C1,FC_C6,p_C6) %>% 
    dplyr::filter(mt %in% tmp$mt) %>% dplyr::mutate(type=type_ls[kk]) %>% 
    dplyr::select(-mt)
  return(res_tmp)
}) %>% do.call(rbind,.)

ggplot(chain_type)+
  geom_boxplot(aes(x=type,y=FC_C1))
ggplot(chain_type)+
  geom_boxplot(aes(x=type,y=FC_C6))+
  ylim(0.8,1.2)+
  RotatedAxis()


###############################################################
##############################################################
carbon_ls = unique(TG_ls$total_C)
carbon_N = lapply(1:length(carbon_ls),function(kk){
  tmp = TG_ls %>% dplyr::filter(total_C==carbon_ls[kk])
  res_tmp = merge_TG %>% dplyr::select(mt,FC_C1,p_C1,FC_C6,p_C6) %>% 
    dplyr::filter(mt %in% tmp$mt) %>% dplyr::mutate(type=carbon_ls[kk]) %>% 
    dplyr::select(-mt)
  return(res_tmp)
}) %>% do.call(rbind,.) %>% arrange(type)
median(carbon_N$type)

carbon_N$type = factor(carbon_N$type,levels = unique(carbon_N$type),ordered = TRUE)
carbon_N_tb = table(carbon_N$type) %>% data.frame() %>% dplyr::rename(type=Var1) 

table(carbon_N$type)
# class(carbon_N$type[1])
carbon_N = carbon_N %>% dplyr::mutate(type_ls = ifelse(type<=53,"<=53",">53")) %>% 
  pivot_longer(cols = c("FC_C1","FC_C6"),names_to = "time",values_to = "FC")


pdf("./plot/TG_compare_carbon.pdf",width = 7,height = 4)
ggplot(carbon_N)+
  geom_bar(data = carbon_N_tb,aes(x=type,y=Freq/165),stat = "identity",alpha=0.5,fill="#B0C4DE")+
  geom_boxplot(aes(x=type,y=FC-0.92,col=time),fill="transparent")+
  scale_color_manual(values = wes_palette("GrandBudapest1",2))+
  geom_hline(yintercept = 1-0.92,lty=2,lwd=0.5,col="purple")+
  ggtitle("carbon")+
  xlab("carbon")+
  scale_y_continuous(limits = c(0,32/165),
                     breaks = seq(0,32/165,8/165),
                     labels = seq(0,32,8),
                     expand = c(0,0),
                     name = "Freq",
                     sec.axis = sec_axis(~.+0.92,breaks = seq(0.95,1.1,0.05),name = "FC"))+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()






##################################################################
##################################################################
unsa_ls = unique(TG_ls$total_unsature)
unsa_N = lapply(1:length(unsa_ls),function(kk){
  tmp = TG_ls %>% dplyr::filter(total_unsature==unsa_ls[kk])
  res_tmp = merge_TG %>% dplyr::select(mt,FC_C1,p_C1,FC_C6,p_C6) %>% 
    dplyr::filter(mt %in% tmp$mt) %>% dplyr::mutate(type=unsa_ls[kk]) %>% 
    dplyr::select(-mt)
  return(res_tmp)
}) %>% do.call(rbind,.)
unsa_N$type = factor(unsa_N$type,levels = unique(unsa_N$type),ordered = TRUE)
unsa_N_tb = table(unsa_N$type) %>% data.frame() %>% dplyr::rename(type=Var1) 


# class(unsa_N$type[1])
unsa_N = unsa_N %>% dplyr::mutate(type_ls = ifelse(type<=4,"<=4",">4")) %>% 
  pivot_longer(cols = c("FC_C1","FC_C6"),names_to = "time",values_to = "FC")


pdf("./plot/TG_compare_unsaturation.pdf",width = 6,height = 4)
ggplot(unsa_N)+
  geom_bar(data = unsa_N_tb,aes(x=type,y=Freq/165),stat = "identity",alpha=0.5,fill="#B0C4DE")+
  geom_boxplot(aes(x=type,y=FC-0.94,col=time),fill="transparent")+
  scale_color_manual(values = wes_palette("GrandBudapest1",2))+
  geom_hline(yintercept = 1-0.94,lty=2,lwd=0.5,col="purple")+
  ggtitle("unsaturation")+
  xlab("unsaturation")+
  scale_y_continuous(limits = c(0,32/165),
                     breaks = seq(0,32/165,8/165),
                     labels = seq(0,32,8),
                     expand = c(0,0),
                     name = "Freq",
                     sec.axis = sec_axis(~.+0.94,breaks = seq(0.95,1.1,0.05),name = "FC"))+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()

median(unsa_N$type)




######################################总结需要验证的mt_ls
SM_all = class_ID_2 %>% dplyr::filter(class=="SM") %>% pull(mt)
Cer_all = class_ID_2 %>% dplyr::filter(grepl("Cer",class)) %>% pull(mt)
SM_sel = mt_list[[5]] %>% .[grepl("SM",.)]
Cer_sel = mt_list[[5]] %>% .[grepl("Cer",.)]

PI_all = class_ID_2 %>% dplyr::filter(class=="PI") %>% pull(mt)
PI_sel = mt_list[[1]] %>% .[grepl("PI",.)]

PEP_all = class_ID_2 %>% dplyr::filter(class=="PE-P") %>% pull(mt)
PEP_sel = mt_list[[3]] %>% .[grepl("PE\\(P",.)]
PEO_all = class_ID_2 %>% dplyr::filter(class=="PE-O") %>% pull(mt)
PEO_sel = mt_list[[1]] %>% .[grepl("PE\\(O",.)]
LPE_all = class_ID_2 %>% dplyr::filter(class=="LPE") %>% pull(mt)
LPE_sel = mt_list[[6]] %>% .[grepl("LPE\\(",.)]

sugar_sel = class_ID_2 %>% dplyr::filter(class=="Sugars") %>% pull(mt) %>% intersect(.,mt_list[[1]])
phenoic_sel = class_ID_2 %>% dplyr::filter(class=="Phenolic acids") %>% pull(mt) %>% intersect(.,mt_list[[6]])

mt_ls_need = list(SM_all,Cer_all,SM_sel,Cer_sel,
                  PI_all,PI_sel,sugar_sel,phenoic_sel,
                  PEP_all,PEP_sel,PEO_all,PEO_sel,LPE_all,LPE_sel)
names(mt_ls_need) = c("SM_all", "Cer_all", "SM_sel", "Cer_sel",
                      "PI_all", "PI_sel", "sugar_sel", "phenoic_sel",
                      "PEP_all","PEP_sel","PEO_all","PEO_sel","LPE_all","LPE_sel")
save(mt_ls_need,file = "./VD/mt_enricher.rdata")




range(TG_ls$total_C)
tt=1
TG_long_ls = lapply(1:(length(unique(TG_ls$total_C))),function(tt){
  tmp = TG_ls %>% dplyr::filter(total_C>=unique(TG_ls$total_C)[tt]) %>% pull(mt)
  return(tmp)
})
names(TG_long_ls)=str_c("TG_long",unique(TG_ls$total_C))
TG_short_ls = lapply(1:(length(unique(TG_ls$total_C))),function(tt){
  tmp = TG_ls %>% dplyr::filter(total_C<=unique(TG_ls$total_C)[tt]) %>% pull(mt)
  return(tmp)
})
names(TG_short_ls)=str_c("TG_short",unique(TG_ls$total_C))
TG_unsature_ls = lapply(1:(length(unique(TG_ls$total_unsature))),function(tt){
  tmp = TG_ls %>% dplyr::filter(total_unsature>=unique(TG_ls$total_unsature)[tt]) %>% pull(mt)
  return(tmp)
})
names(TG_unsature_ls)=str_c("TG_unsature",unique(TG_ls$total_unsature))
TG_sature_ls = lapply(1:(length(unique(TG_ls$total_unsature))),function(tt){
  tmp = TG_ls %>% dplyr::filter(total_unsature<=unique(TG_ls$total_unsature)[tt]) %>% pull(mt)
  return(tmp)
})
names(TG_sature_ls)=str_c("TG_sature",unique(TG_ls$total_unsature))

save(TG_long_ls,TG_short_ls,TG_sature_ls,TG_unsature_ls,file = "./VD/mt_enricher_TG.rdata")
