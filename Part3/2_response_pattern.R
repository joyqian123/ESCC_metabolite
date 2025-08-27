rm(list = ls())
setwd("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/")
library(tidyverse)
dir.create("metabo_process")
setwd("./metabo_process/")
dir.create("./all_C1C6")
setwd("./all_C1C6/")

load("~/help_for_others/DYQ/output/1_preclean/data/data_preclean.rdata")
table(clinic_in$cycle)



clinic_1 = clinic_in %>% dplyr::filter(cycle %in% c("C1D1","C6D1"))


# .libPaths()
library(MetaboAnalystR)
library(tidyverse)

mSet<-InitDataObjects("conc", "stat", FALSE)
mSet<-Read.TextData(mSet, "../../csv/data_all_C1C6.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet)  ####替换0值或缺失值
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE,ratioNum = 20)

mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###按代谢物作图
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###按样本作图

norm = mSet$dataSet$norm
info = mSet$dataSet$meta.info %>% unlist()

norm_1 = norm %>% rownames_to_column("Sample") %>% 
  inner_join(clinic_1[,c("Sample","SUBJID","cycle")],by = "Sample") %>% 
  dplyr::select(Sample,SUBJID,cycle,everything())
# norm_2 = norm_1 %>% pivot_longer(cols = !c(1:3),names_to = "mt",values_to="value")

norm_2_C1 = norm_1%>% dplyr::filter(cycle=="C1D1")  %>% dplyr::select(-c(Sample,cycle)) 
norm_2_C6 = norm_1 %>% dplyr::filter(cycle=="C6D1")  %>% dplyr::select(-c(Sample,cycle)) 


clinic_C1 = clinic_1[,c(2,6,7,9,11,12,13,14,18,19,20,21,25,26,35)] %>% distinct(.,.keep_all = T) %>% 
  dplyr::filter(SUBJID %in% norm_2_C1$SUBJID)

clinic_C6 = clinic_1[,c(2,6,7,9,11,12,13,14,18,19,20,21,25,26,35)] %>% distinct(.,.keep_all = T) %>% 
  dplyr::filter(SUBJID %in% norm_2_C6$SUBJID)

  mt_expr = norm_2_C1 %>% column_to_rownames("SUBJID")

  clinic_info_1 = clinic_C1 %>% dplyr::select(1:12) %>% column_to_rownames("SUBJID")
  clinic_info_3 = clinic_C1 %>% dplyr::select(1,15) %>% column_to_rownames("SUBJID")

  source("~/help_for_others/DYQ/code/8_compare_with_mt_change_in_R_and_NR/VD.R")
  
  vd.vars <- c(colnames(clinic_info_1),colnames(clinic_info_3))
  meta.data <- cbind(clinic_info_1,clinic_info_3)
  table(meta.data$R)
  
  # unknown = meta.data %>% dplyr::filter(PFS_larger=="unknown") %>% rownames(.)
  
  ras.data <- mt_expr #%>% rownames_to_column("ID") %>% dplyr::filter(!ID %in% unknown) %>% column_to_rownames("ID")
  meta.data = meta.data #%>% rownames_to_column("ID") %>% dplyr::filter(!ID %in% unknown) %>% column_to_rownames("ID")
  
  vd.res_def_C1 <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = 5)
  
  
  FC_tb = norm_2_C1 %>% column_to_rownames("SUBJID") %>% cbind(clinic_info_3,.) 
  
  # tt=1057
  FC_res_C1 = lapply(2:ncol(FC_tb),function(tt){
    print(tt)
    tmp = FC_tb %>% dplyr::select(R,any_of(tt))
    if(length(unique(tmp[,2]))>1){
      t_res = t.test(tmp[,2]~tmp$R)
      return(data.frame(mt=colnames(FC_tb)[tt],t=t_res$statistic,
                        FC=t_res$estimate["mean in group R"]/t_res$estimate["mean in group NR"],p=t_res$p.value))      
    }
  }) %>% do.call(rbind,.)
  

  
#################################################################################
#################################################################################
#################################################################################
  mt_expr = norm_2_C6 %>% column_to_rownames("SUBJID")

  clinic_info_1 = clinic_C6 %>% dplyr::select(1:12) %>% column_to_rownames("SUBJID")
  clinic_info_3 = clinic_C6 %>% dplyr::select(1,15) %>% column_to_rownames("SUBJID")

  source("~/help_for_others/DYQ/code/8_compare_with_mt_change_in_R_and_NR/VD.R")
  
  vd.vars <- c(colnames(clinic_info_1),colnames(clinic_info_3))
  meta.data <- cbind(clinic_info_1,clinic_info_3)
  table(meta.data$R)
  

  ras.data <- mt_expr 
  meta.data = meta.data
  
  vd.res_def_C6 <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = 5)
  
  
  
  FC_tb = norm_2_C6 %>% column_to_rownames("SUBJID") %>% cbind(clinic_info_3,.) 
  
  # tt=1057
  FC_res_C6 = lapply(2:ncol(FC_tb),function(tt){
    print(tt)
    tmp = FC_tb %>% dplyr::select(R,any_of(tt))
    if(length(unique(tmp[,2]))>1){
      t_res = t.test(tmp[,2]~tmp$R)
      return(data.frame(mt=colnames(FC_tb)[tt],t=t_res$statistic,
                        FC=t_res$estimate["mean in group R"]/t_res$estimate["mean in group NR"],p=t_res$p.value))      
    }
  }) %>% do.call(rbind,.)
  
  
###################################################################################
###################################################################################
###################################################################################
  norm_1 = norm %>% rownames_to_column("Sample") %>% 
    inner_join(clinic_1[,c("Sample","SUBJID","cycle")],by = "Sample") %>% 
    dplyr::select(Sample,SUBJID,cycle,everything())
  norm_2 = norm_1 %>% pivot_longer(cols = !c(1:3),names_to = "mt",values_to="value")
  
  norm_3_C1 = norm_2 %>% dplyr::select(-Sample) %>% dplyr::filter(cycle=="C1D1") %>% 
    arrange(SUBJID,mt) %>% dplyr::rename(C1D1=value)
  norm_3_C6 = norm_2 %>% dplyr::select(-Sample) %>% dplyr::filter(cycle=="C6D1") %>% 
    arrange(SUBJID,mt) %>% dplyr::rename(C6D1=value)
  
  norm_3 = norm_3_C1[,c(1,3,4)] %>% inner_join(norm_3_C6[,c(1,3,4)],by = c("SUBJID","mt"))
  norm_3$delta = norm_3$C6D1-norm_3$C1D1
  norm_3_process = norm_3 %>% dplyr::select(1,2,5) %>% pivot_wider(id_cols = c(1),names_from = "mt",values_from = "delta")
  
  
  
  clinic_need = clinic_1[,c(2,6,7,9,11,12,13,14,18,19,20,21,25,26,35)] %>% distinct(.,.keep_all = T) %>% 
    dplyr::filter(SUBJID %in% norm_3_process$SUBJID)
  
  mt_expr = norm_3_process %>% column_to_rownames("SUBJID")
  clinic_info_1 = clinic_need %>% dplyr::select(1:12) %>% column_to_rownames("SUBJID")
  clinic_info_3 = clinic_need %>% dplyr::select(1,15) %>% column_to_rownames("SUBJID")
  
  source("~/help_for_others/DYQ/code/8_compare_with_mt_change_in_R_and_NR/VD.R")
  
  vd.vars <- c(colnames(clinic_info_1),colnames(clinic_info_3))
  meta.data <- cbind(clinic_info_1,clinic_info_3)
  table(meta.data$R)
  

  ras.data <- mt_expr 
  meta.data = meta.data 
  
  vd.res_def_delta <- VarDecompose(data = ras.data, meta.data = meta.data, vd.vars = vd.vars, cores = 5)

  
  FC_tb = norm_3_process %>% column_to_rownames("SUBJID") %>% cbind(clinic_info_3,.) 
  
  # tt=1057
  tt=2
  FC_res_delta = lapply(2:ncol(FC_tb),function(tt){
    print(tt)
    tmp = FC_tb %>% dplyr::select(R,any_of(tt))
    if(length(unique(tmp[,2]))>1){
      t_res = t.test(tmp[,2]~tmp$R)
      return(data.frame(mt=colnames(FC_tb)[tt],t=t_res$statistic,
                        R_delta = t_res$estimate["mean in group R"],
                        NR_delta = t_res$estimate["mean in group NR"],
                        FC=t_res$estimate["mean in group R"]-t_res$estimate["mean in group NR"],p=t_res$p.value))      
    }
  }) %>% do.call(rbind,.)
  
  
  
  save(vd.res_def_C1,vd.res_def_C6,vd.res_def_delta,
       FC_res_C1,FC_res_C6,FC_res_delta,
       file = "~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/VD/only_response_VD_res.rdata")
  
  
  
  load("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/VD/only_response_VD_res.rdata")
########################################################################################################
  ######################################################################################################
  plot1 = vd.res_def_C1 %>% dplyr::select(R) %>% rownames_to_column("mt")  %>% dplyr::mutate(data="C1")
  plot2 = vd.res_def_C6 %>% dplyr::select(R) %>% rownames_to_column("mt")  %>% dplyr::mutate(data="C6")
  plot3 = vd.res_def_delta %>% dplyr::select(R) %>% rownames_to_column("mt")  %>% dplyr::mutate(data="C6-C1")

  
  plot = rbind(plot1,plot2,plot3) %>% dplyr::filter(R>0.05)  
  
  library(ggridges)
  library(wesanderson)
  pdf("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/plot/response_VD_with_C1_C6_delta.pdf",width = 5,height = 4)
  ggplot(plot, aes(x=R, y=data,fill=data))+
    geom_density_ridges()+
    scale_fill_manual(values = wes_palette(name = "GrandBudapest1",3))+
    scale_x_log10()+
    labs(x="variance of explanation")+
    theme_ridges(font_size = 13, grid = FALSE)+
    theme(axis.title.y = element_blank())#+
    # theme_bw()
  dev.off()
  
  
  sel_mt = rbind(plot1,plot2,plot3) %>% dplyr::filter(R>0) %>% pull(mt)
  
  merge = rbind(plot1,plot2,plot3) %>% dplyr::filter(mt %in% sel_mt) %>% 
    pivot_wider(id_cols = mt,names_from = "data",values_from = "R")

  colnames(FC_res_C1)[2:4]=str_c(colnames(FC_res_C1)[2:4],"C1",sep = "_")
  colnames(FC_res_C6)[2:4]=str_c(colnames(FC_res_C6)[2:4],"C6",sep = "_")
  colnames(FC_res_delta)[4:6]=str_c(colnames(FC_res_delta)[4:6],"delta",sep = "_")

  
  merge = merge %>% inner_join(FC_res_C1,by = "mt") %>% inner_join(FC_res_C6,by = "mt") %>% 
    inner_join(FC_res_delta,by = "mt")
  
  
  
##############################################################################
  sig_at_least_one = merge %>% dplyr::filter(p_C1<0.1|p_C6<0.1)
  
  mt_continuous_R = sig_at_least_one %>% dplyr::filter(p_C1<0.1,p_C6<0.1,log(FC_C1)*log(FC_C6)>0,log(FC_C1)>0) %>% pull(mt)
  mt_continuous_NR = sig_at_least_one %>% dplyr::filter(p_C1<0.1,p_C6<0.1,log(FC_C1)*log(FC_C6)>0,log(FC_C1)<0) %>% pull(mt)
  mt_baseline_R = sig_at_least_one %>% dplyr::filter(p_C1<0.1,p_C6>=0.1,log(FC_C1)>0) %>% pull(mt)
  mt_baseline_NR = sig_at_least_one %>% dplyr::filter(p_C1<0.1,p_C6>=0.1,log(FC_C1)<0) %>% pull(mt)
  mt_response_R = sig_at_least_one %>% dplyr::filter(p_C1>=0.1,p_C6<0.1,log(FC_C6)>0) %>% pull(mt)
  mt_response_NR = sig_at_least_one %>% dplyr::filter(p_C1>=0.1,p_C6<0.1,log(FC_C6)<0) %>% pull(mt)
  
  
  
  mt_list = list(mt_continuous_R,mt_continuous_NR,mt_baseline_R,mt_baseline_NR,mt_response_R,mt_response_NR)
  names(mt_list) = c("continuous_R","continuous_NR","baseline_R","baseline_NR","response_R","response_NR")
  
save(mt_list,file = "../../VD/mt_list_from_data.rdata")  
  
  

mt_continuous_R = sig_at_least_one %>% dplyr::filter(p_C1<0.1,p_C6<0.1,log(FC_C1)*log(FC_C6)>0,log(FC_C1)>0) %>% dplyr::mutate(class="continuous_R")
mt_continuous_NR = sig_at_least_one %>% dplyr::filter(p_C1<0.1,p_C6<0.1,log(FC_C1)*log(FC_C6)>0,log(FC_C1)<0) %>% dplyr::mutate(class="continuous_NR")
mt_baseline_R = sig_at_least_one %>% dplyr::filter(p_C1<0.1,p_C6>=0.1,log(FC_C1)>0) %>% dplyr::mutate(class="baseline_R")
mt_baseline_NR = sig_at_least_one %>% dplyr::filter(p_C1<0.1,p_C6>=0.1,log(FC_C1)<0) %>% dplyr::mutate(class="baseline_NR")
mt_response_R = sig_at_least_one %>% dplyr::filter(p_C1>=0.1,p_C6<0.1,log(FC_C6)>0) %>% dplyr::mutate(class="response_R")
mt_response_NR = sig_at_least_one %>% dplyr::filter(p_C1>=0.1,p_C6<0.1,log(FC_C6)<0)%>% dplyr::mutate(class="response_NR")
all_mt = rbind(mt_continuous_R,mt_continuous_NR,mt_baseline_R,mt_baseline_NR,mt_response_R,mt_response_NR)
# others = sig_at_least_one %>% dplyr::filter(!mt %in% all_mt$mt) %>% dplyr::mutate(class="others")

# all_mt = rbind(all_mt,others)

p_mt = all_mt %>% dplyr::select(mt,p_C1,p_C6) %>% column_to_rownames("mt") %>% log()*(-1)
C1_FC = all_mt %>% dplyr::select(mt,FC_C1,p_C1) %>% column_to_rownames("mt") 
C1_FC$FC_C1 = log(C1_FC$FC_C1)
C1_FC$FC_C1 = sapply(1:nrow(C1_FC),function(kk){
  if(C1_FC$FC_C1[kk]>0.2){
    return(0.2)
  }else if(C1_FC$FC_C1[kk]<(-0.1)){
    return(-0.1)
  }else{
    return(C1_FC$FC_C1[kk])
  }
})
C6_FC = all_mt %>% dplyr::select(mt,FC_C6,p_C6) %>% column_to_rownames("mt") 
C6_FC$FC_C6 = log(C6_FC$FC_C6)
C6_FC$FC_C6 = sapply(1:nrow(C6_FC),function(kk){
  if(C6_FC$FC_C6[kk]>0.2){
    return(0.2)
  }else if(C6_FC$FC_C6[kk]<(-0.1)){
    return(-0.1)
  }else{
    return(C6_FC$FC_C6[kk])
  }
})

library(circlize)
library(circlize)
library(RColorBrewer)
library(stringr)
library(ComplexHeatmap)
library(grImport2)
library(gridBase)
library(ggplot2)
col_fun1 = colorRamp2(c(0, -log(0.1),40), c("#6A5ACD", "#FFF0F5", "#FF69B4"))
# 
col_anno_gene = structure(c("#C71585","#483D8B","#FF69B4","#7B68EE","#FFB6C1","#B0C4DE"), names = unique(all_mt$class))
# col_anno_gene = structure(c("#e5965a","#746fb1","#f0c2a2","#b5a1e3","#FFB6C1","#B0C4DE"), names = unique(all_mt$class))

circlize_plot<-function(){
  circos.clear()
  circos.par(gap.degree = c(rep(2,5),30), cell.padding = c(0.01, 0, 0.01, 0),
             track.margin = c(0.01,0.01),start.degree=90)
  circos.heatmap.initialize(p_mt, split = all_mt$class)
  
  
  circos.heatmap(all_mt$class,split = all_mt$class, col = col_anno_gene, track.height = 0.05) 
  circos.heatmap(p_mt, split = all_mt$class,col = col_fun1,track.height = 0.1,show.sector.labels = TRUE)
 
  # CELL_META$row_dend
  # CELL_META$row_order
  # CELL_META$subset
  
  circos.track(track.height = 0.3,ylim=c(-0.1,0.2),
               panel.fun = function(x, y) {
                 y = C1_FC$FC_C1[CELL_META$subset]
                 y = y[CELL_META$row_order]
                 circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
                 # circos.points(seq_along(y) - 0.5, y, col = ifelse(y > 0, "red", "blue"))
                 circos.lines(seq_along(y)-0.5, y, baseline = 0,col=ifelse(C1_FC$p_C1[CELL_META$subset]>0.1,"grey",
                                                                           ifelse(C1_FC$FC_C1[CELL_META$subset]>0,"#FF69B4","#6A5ACD")),
                              type = "h",lwd = 2)
               })
  # "#FF69B4","#6A5ACD"
  circos.track(track.height = 0.3,ylim=c(-0.1,0.2),
               panel.fun = function(x, y) {
                 y = C6_FC$FC_C6[CELL_META$subset]
                 y = y[CELL_META$row_order]
                 circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
                 # circos.points(seq_along(y) - 0.5, y, col = ifelse(y > 0, "red", "blue"))
                 circos.lines(seq_along(y)-0.5, y, baseline = 0,col=ifelse(C6_FC$p_C6[CELL_META$subset]>0.1,"grey",
                                                                           ifelse(C6_FC$FC_C6[CELL_META$subset]>0,"#FF69B4","#6A5ACD")),
                              type = "h",lwd = 2)
               })
  # "#FF69B4","#6A5ACD"
  circos.yaxis("left", at=c(-0.1,0,0.1,0.2),labels = c("<-0.1","0","0.1",">0.2"), track.index =  3, 
               sector.index = "continuous_R",
               labels.cex = 0.5,tick=T,labels.niceFacing=F)
  circos.yaxis("left", at=c(-0.1,0,0.1,0.2),labels = c("<-0.1","0","0.1",">0.2"), track.index =  4, 
               sector.index = "continuous_R",
               labels.cex = 0.5,tick=T,labels.niceFacing=F)
  
}


lgd_NES = Legend(title = "logFC", at = c("up_in_R","up_in_NR","unsig"), 
                 legend_gp = gpar(fill = c("#FF69B4","#6A5ACD","grey")))
# 

lgd_NES_2 = Legend(title = "-log(P value)", col_fun = col_fun1)


pdf("../../plot/classify_of_mt_into_6_types_new.pdf",width = 8,height = 5)
library(gridBase)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region

pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circlize_plot()  ####把前面的circlize定义为函数
upViewport()

h = dev.size()[2]
lgd_list = packLegend(lgd_NES,lgd_NES_2, max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = circle_size, just = "left")
dev.off()
  
  
  
