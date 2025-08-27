rm(list = ls())
setwd("~/help_for_others/DYQ/output/")
load("./8_compare_with_mt_change_in_R_and_NR/VD/mt_list_from_data.rdata")

load("./7_JM/joint_res_1/joint_tb_and_res_1.rdata")
res_1 = res %>% dplyr::select(metabo,coef,p)
colnames(res_1)[2:3]=str_c(colnames(res_1)[2:3],"R1",sep = "_")
res_1$coef_R1 = log(res_1$coef_R1)
load("./7_JM/joint_res_2/joint_tb_and_res_2.rdata")
res_2 = res %>% dplyr::select(metabo,coef,p)
colnames(res_2)[2:3]=str_c(colnames(res_2)[2:3],"R2",sep = "_")

res = res_1 %>% inner_join(res_2,by = "metabo") %>% na.omit()





##########################################################
contin_R = c(mt_list[[1]])

res_contin_R = res %>% dplyr::filter(metabo %in% contin_R)%>% dplyr::filter(coef_R1<0,coef_R2<0)

response_R = c(mt_list[[6]])

res_response_R = res %>% dplyr::filter(metabo %in% response_R)%>% dplyr::filter(coef_R1<0,coef_R2<0)

time_R = res %>% dplyr::filter(coef_R1<0,coef_R2<0) %>% pull(metabo)


library(eulerr)

fit <- euler(list(
  continuous_R_mt = contin_R,
  response_R_mt = response_R,
  time_dependent_R_mt = time_R
))


pdf("./8_compare_with_mt_change_in_R_and_NR/good_mt/good_mt_intersection.pdf",width = 4,height = 4)
plot(fit,
     fills = list(fill = c("#B0C4DE", "#778899", "#4682B4"), alpha = 0.4),
     labels = list(font = 2),
     edges = list(
       col = "black",     # 边框颜色
       lwd = 1,           # 边框线宽
       lty = 2            # 线型：1 = 实线，2 = 虚线，3 = 点线等
     ),
     quantities = TRUE)
dev.off()



load("~/help_for_others/DYQ/code/8_compare_with_mt_change_in_R_and_NR/HMDB/metabolite_type.rdata")

load("./1_preclean/data/data_preclean.rdata")
mt_info = data_in %>% dplyr::select(1:4)

R_with_cont_R_ID = mt_info %>% dplyr::filter(Compounds %in% res_contin_R$metabo)

R_with_cont_R_ID$HMDB[5]="HMDB0034323"
R_with_cont_R_ID$HMDB[11]="HMDB0011359"
R_with_cont_R_ID$HMDB[12]="HMDB0011392"
R_with_cont_R_ID$HMDB[24]="HMDB0011360"

R_with_cont_R_ID = R_with_cont_R_ID %>% dplyr::filter(!HMDB=="-") %>% 
  dplyr::mutate(exogenous=ifelse(HMDB %in% exo$HMDB_ID,TRUE,FALSE)) %>% 
  dplyr::mutate(food=ifelse(HMDB %in% food$HMDB_ID,TRUE,FALSE)) %>% 
  dplyr::mutate(plant=ifelse(HMDB %in% plant$HMDB_ID,TRUE,FALSE)) %>%
  dplyr::mutate(drug=ifelse(HMDB %in% drug$HMDB_ID,TRUE,FALSE)) 
R_with_cont_R_ID$exogenous[1]=TRUE  ###官网修正  
R_with_cont_R_ID$food[1]=TRUE
R_with_cont_R_ID$exogenous[2]=TRUE  ###官网修正  
R_with_cont_R_ID$food[2]=TRUE
R_with_cont_R_ID$exogenous[3]=TRUE  ###官网修正  
R_with_cont_R_ID$food[3]=TRUE
R_with_cont_R_ID$exogenous[4]=TRUE  ###官网修正  
R_with_cont_R_ID$food[4]=TRUE

R_with_cont_R_ID$plant[5]=TRUE

R_with_cont_R_ID$endogenous = c(TRUE,TRUE,TRUE,TRUE,FALSE,rep(TRUE,13))



R_with_cont_R_ID = R_with_cont_R_ID %>% dplyr::select(Compounds,endogenous,exogenous,food,plant,drug) %>% 
  column_to_rownames("Compounds")
library(ComplexHeatmap)
R_with_cont_R_ID[R_with_cont_R_ID==TRUE]=1
R_with_cont_R_ID[R_with_cont_R_ID==FALSE]=0
library(circlize)
col_fun = colorRamp2(c(0, 0.5, 1), c("lightgrey", "white", "#6495ED"))
pdf("./8_compare_with_mt_change_in_R_and_NR/good_mt/R_with_cont_R.pdf",width = 5,height = 7)
Heatmap(as.matrix(R_with_cont_R_ID),border = TRUE,col = col_fun,cluster_columns = FALSE,
        width = unit(50, "mm"),  # 每个单元格宽度
        height = unit(100, "mm"), # 每个单元格高度
        cluster_rows = FALSE,
        rect_gp = gpar(lwd = 2, col = "white"))
dev.off()








##############################################################
##############################################################
R_with_resp_R_ID = mt_info %>% dplyr::filter(Compounds %in% res_response_R$metabo)

R_with_resp_R_ID$HMDB[9]="HMDB0114742"
R_with_resp_R_ID$HMDB[10]="HMDB0010383"
R_with_resp_R_ID$HMDB[11]="HMDB0010387"
R_with_resp_R_ID$HMDB[12]="HMDB0011479"
R_with_resp_R_ID$HMDB[14]="HMDB0059745"
R_with_resp_R_ID$HMDB[18]="HMDB0031189"
R_with_resp_R_ID$HMDB[20]="HMDB0000722"
R_with_resp_R_ID$HMDB[21]="HMDB0028829"
# R_with_resp_R_ID$HMDB[27]="HMDB0341549"

R_with_resp_R_ID = R_with_resp_R_ID %>% dplyr::filter(!HMDB=="-") %>% 
  dplyr::mutate(exogenous=ifelse(HMDB %in% exo$HMDB_ID,TRUE,FALSE)) %>% 
  dplyr::mutate(food=ifelse(HMDB %in% food$HMDB_ID,TRUE,FALSE)) %>% 
  dplyr::mutate(plant=ifelse(HMDB %in% plant$HMDB_ID,TRUE,FALSE)) %>%
  dplyr::mutate(drug=ifelse(HMDB %in% drug$HMDB_ID,TRUE,FALSE)) 


# R_with_resp_R_ID$exogenous[1]=TRUE  ###官网修正  
R_with_resp_R_ID$food[1]=TRUE  ###1同时来自微生物（肠道菌群）
R_with_resp_R_ID$food[2]=TRUE  
R_with_resp_R_ID$food[3]=TRUE  ###3同时来自微生物（肠道菌群）
R_with_resp_R_ID$food[4]=TRUE  

R_with_resp_R_ID$exogenous[9]=TRUE  
R_with_resp_R_ID$food[9]=TRUE  
R_with_resp_R_ID$plant[11]=TRUE  ###11同时来自微生物（肠道菌群）
R_with_resp_R_ID$exogenous[12]=TRUE  
R_with_resp_R_ID$plant[13]=TRUE
R_with_resp_R_ID$exogenous[15]=TRUE  
R_with_resp_R_ID$food[15]=TRUE 

R_with_resp_R_ID$endogenous = c(TRUE,TRUE,TRUE,FALSE,rep(TRUE,4),TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,rep(TRUE,24))



R_with_resp_R_ID = R_with_resp_R_ID %>% dplyr::select(Compounds,endogenous,exogenous,food,plant,drug) %>% 
  column_to_rownames("Compounds")
library(ComplexHeatmap)
R_with_resp_R_ID[R_with_resp_R_ID==TRUE]=1
R_with_resp_R_ID[R_with_resp_R_ID==FALSE]=0
library(circlize)
col_fun = colorRamp2(c(0, 0.5, 1), c("lightgrey", "white", "#6495ED"))
pdf("./8_compare_with_mt_change_in_R_and_NR/good_mt/R_with_resp_R.pdf",width = 5.5,height = 7)
Heatmap(as.matrix(R_with_resp_R_ID)[1:20,],border = TRUE,col = col_fun,cluster_columns = FALSE,
        width = unit(50, "mm"),  # 每个单元格宽度
        height = unit(100, "mm"), # 每个单元格高度
        cluster_rows = FALSE,
        rect_gp = gpar(lwd = 2, col = "white"))
Heatmap(as.matrix(R_with_resp_R_ID)[21:40,],border = TRUE,col = col_fun,cluster_columns = FALSE,
        width = unit(50, "mm"),  # 每个单元格宽度
        height = unit(100, "mm"), # 每个单元格高度
        cluster_rows = FALSE,
        rect_gp = gpar(lwd = 2, col = "white"))
dev.off()







#############################################################
############################################################
# rm(list = ls())
load("~/help_for_others/DYQ/output/1_preclean/data/data_preclean.rdata")
table(clinic_in$cycle)



clinic_1 = clinic_in %>% dplyr::filter(cycle %in% c("C1D1","C6D1"))
data_1 = data_in %>% dplyr::select(1:4,any_of(clinic_1$Sample))

clinic_need = clinic_1 %>% dplyr::select(SUBJID,R) %>% distinct(.,.keep_all = T)

setwd("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/")
dir.create("metabo_process")
setwd("./metabo_process/")
dir.create("./all_C1C6")
setwd("./all_C1C6/")


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
norm_2 = norm_1 %>% pivot_longer(cols = !c(1:3),names_to = "mt",values_to="value")

norm_3_C1 = norm_2 %>% dplyr::select(-Sample) %>% dplyr::filter(cycle=="C1D1") %>% 
  arrange(SUBJID,mt) %>% dplyr::rename(C1D1=value)
norm_3_C6 = norm_2 %>% dplyr::select(-Sample) %>% dplyr::filter(cycle=="C6D1") %>% 
  arrange(SUBJID,mt) %>% dplyr::rename(C6D1=value)

norm_3 = norm_3_C1[,c(1,3,4)] %>% inner_join(norm_3_C6[,c(1,3,4)],by = c("SUBJID","mt")) %>% 
  inner_join(clinic_need,by="SUBJID")




sel_mt_1 = R_with_cont_R_ID %>% dplyr::filter(endogenous==FALSE) %>% rownames_to_column("mt") %>% pull(mt)
sel_mt_2 = R_with_resp_R_ID %>% dplyr::filter(endogenous==FALSE) %>% rownames_to_column("mt") %>% pull(mt)
sel_mt = c(sel_mt_1,sel_mt_2[1])
library(ggpubr)
library(ggpmisc)
kk=1
plot_ls = lapply(1:length(sel_mt),function(kk){
  print(kk)
  tmp = norm_3 %>% dplyr::filter(mt==sel_mt[kk]) %>% pivot_longer(cols = c("C1D1","C6D1"),names_to = "time",values_to = "value") %>% 
    dplyr::mutate(group=str_c(R,time,sep = "_"))
  fl = aov(value ~ time * R + Error(SUBJID/(time)), data = tmp)
  a = summary(fl)
  median_tb <- tmp %>%
    group_by(group) %>%
    dplyr::summarise(median_expr = median(value, na.rm = TRUE)) %>%
    ungroup()
  
  p=ggplot(tmp,aes(x=group,y=value))+
    geom_violin(aes(fill=time),alpha=1,scale = "width",trim = TRUE)+
    # geom_boxplot(aes(color=response),alpha=0.5,width=0.1,outliers=FALSE)+
    stat_compare_means(comparisons = list(c("NR_C1D1","NR_C6D1"),c("R_C1D1","R_C6D1"),
                                          c("NR_C1D1","R_C1D1"),c("NR_C6D1","R_C6D1")),method = "t.test")+
    geom_segment(data = median_tb,  # 使用提前计算的中位数数据
                 aes(x= as.numeric(as.factor(group)) - 0.15,  # 调整x的值来控制线段的起点
                     xend = as.numeric(as.factor(group)) + 0.15,  # 调整xend的值来控制线段的终点
                     y = median_expr,
                     yend = median_expr),
                 linewidth=1,lty=1) +
    scale_fill_manual(values = c("#b5a1e3","#f0c2a2"))+
    scale_color_manual(values = c("#b5a1e3","#f0c2a2"))+
    xlab("abundance")+
    # facet_wrap(~mt,ncol=5,scale="free_y")+
    theme_bw() +
    theme(panel.grid.major=element_line(colour=NA),
          strip.background = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.miNRr = element_blank())+
    RotatedAxis()+
    ggtitle(sel_mt[kk])
  
  return(p)
})

library(patchwork)
pdf("../../good_mt//final_sel_mt.pdf",width = 10,height = 5)
wrap_plots(plot_ls,nrow = 1)
dev.off()









load("../../../1_preclean/data/data_preclean.rdata")
data_use = data_in %>% dplyr::select(-c(2:4)) %>% dplyr::filter(Compounds %in% sel_mt) %>% 
  column_to_rownames("Compounds") %>% t() %>% data.frame(.,check.names = FALSE) %>% 
  rownames_to_column("Sample") %>% 
  inner_join(clinic_in[,c("Sample","SUBJID","cycle","AVAL1","CNSR1","R")],by = "Sample")

data_use$CNSR1 = ifelse(data_use$CNSR1==0,1,0)

load("../../../../data/blood_time.rdata")
blood_time = blood_time %>% dplyr::mutate(cycle=ifelse(Cycle=="筛选期(D-7~D-1)","C1D1",Cycle))

data_use = data_use %>% inner_join(blood_time[,c("SUBJID","cycle","Time")],by = c("SUBJID","cycle")) %>% 
  dplyr::mutate(measure_time=ifelse(Time<0,0,Time)) %>% dplyr::mutate(measure_time=measure_time/30.44)


RRA_PL = data_use %>% dplyr::mutate(Score=log(`Indole 3-carbinol`))
RRA_PL$cycle = factor(RRA_PL$cycle,levels = c("C1D1","C6D1","C13D1","C25D1","C33D1"),ordered = TRUE)

res = lapply(4:24,function(tt){
  
  RRA_PL$type = sapply(1:nrow(RRA_PL),function(kk){
    if(RRA_PL$AVAL1[kk]<tt & RRA_PL$CNSR1[kk]==1){
      return("short")
    }else if(RRA_PL$AVAL1[kk]<tt & RRA_PL$CNSR1[kk]==0){
      return("unknown")
    }else if(RRA_PL$AVAL1[kk]>=tt){
      #   return("long no progression(>=6m)")
      # }else{
      return("long")
    }
  })
  
  table(RRA_PL$type)
  RRA_PL_plot = RRA_PL %>% dplyr::filter(!type=="unknown") %>% dplyr::filter(measure_time<tt)
  RRA_PL_plot$cycle = factor(RRA_PL_plot$cycle,levels = c("C1D1","C6D1","C13D1","C25D1","C33D1"),ordered = TRUE)
  RRA_PL_plot = RRA_PL_plot %>% group_by(SUBJID) %>% dplyr::mutate(mean_s=mean(Score)) %>% 
    distinct(SUBJID,.keep_all = T) %>% dplyr::mutate(pred_time=tt)
  
  return(RRA_PL_plot)

}) %>% do.call(rbind,.)


median_tb <- res %>%
  group_by(type,pred_time) %>%
  dplyr::summarise(median_expr = median(mean_s, na.rm = TRUE)) %>%
  ungroup()

p=ggplot(res,aes(x=type,y=mean_s))+
  geom_violin(aes(fill=type),alpha=1,scale = "width",trim = TRUE)+
  # geom_boxplot(aes(color=response),alpha=0.5,width=0.1,outliers=FALSE)+
  stat_compare_means(comparisons = list(c("long","short")),method = "t.test")+
  geom_segment(data = median_tb,  # 使用提前计算的中位数数据
               aes(x= as.numeric(as.factor(type)) - 0.15,  # 调整x的值来控制线段的起点
                   xend = as.numeric(as.factor(type)) + 0.15,  # 调整xend的值来控制线段的终点
                   y = median_expr,
                   yend = median_expr),
               linewidth=1,lty=1) +
  scale_fill_manual(values = c("#b5a1e3","#f0c2a2"))+
  scale_color_manual(values = c("#b5a1e3","#f0c2a2"))+
  xlab("abundance")+
  # ggtitle(paste(tt,"m",sep = ""))+
  facet_wrap(~pred_time,nrow=3,scale="free_y")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        legend.position = "none",
        strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

pdf("../../good_mt/indo_follow.pdf",width = 12,height = 8)
# wrap_plots(plot_ls,nrow = 3)
p
dev.off()






RRA_PL = data_use %>% dplyr::mutate(Score=log(`S-Allyl-L-cysteine`))
RRA_PL$cycle = factor(RRA_PL$cycle,levels = c("C1D1","C6D1","C13D1","C25D1","C33D1"),ordered = TRUE)

res = lapply(4:24,function(tt){
  
  RRA_PL$type = sapply(1:nrow(RRA_PL),function(kk){
    if(RRA_PL$AVAL1[kk]<tt & RRA_PL$CNSR1[kk]==1){
      return("short")
    }else if(RRA_PL$AVAL1[kk]<tt & RRA_PL$CNSR1[kk]==0){
      return("unknown")
    }else if(RRA_PL$AVAL1[kk]>=tt){
      #   return("long no progression(>=6m)")
      # }else{
      return("long")
    }
  })
  
  table(RRA_PL$type)
  RRA_PL_plot = RRA_PL %>% dplyr::filter(!type=="unknown") %>% dplyr::filter(measure_time<tt)
  RRA_PL_plot$cycle = factor(RRA_PL_plot$cycle,levels = c("C1D1","C6D1","C13D1","C25D1","C33D1"),ordered = TRUE)
  RRA_PL_plot = RRA_PL_plot %>% group_by(SUBJID) %>% dplyr::mutate(mean_s=mean(Score)) %>% 
    distinct(SUBJID,.keep_all = T) %>% dplyr::mutate(pred_time=tt)
  
  return(RRA_PL_plot)
  
}) %>% do.call(rbind,.)


median_tb <- res %>%
  group_by(type,pred_time) %>%
  dplyr::summarise(median_expr = median(mean_s, na.rm = TRUE)) %>%
  ungroup()

p=ggplot(res,aes(x=type,y=mean_s))+
  geom_violin(aes(fill=type),alpha=1,scale = "width",trim = TRUE)+
  # geom_boxplot(aes(color=response),alpha=0.5,width=0.1,outliers=FALSE)+
  stat_compare_means(comparisons = list(c("long","short")),method = "t.test")+
  geom_segment(data = median_tb,  # 使用提前计算的中位数数据
               aes(x= as.numeric(as.factor(type)) - 0.15,  # 调整x的值来控制线段的起点
                   xend = as.numeric(as.factor(type)) + 0.15,  # 调整xend的值来控制线段的终点
                   y = median_expr,
                   yend = median_expr),
               linewidth=1,lty=1) +
  scale_fill_manual(values = c("#b5a1e3","#f0c2a2"))+
  scale_color_manual(values = c("#b5a1e3","#f0c2a2"))+
  xlab("abundance")+
  # ggtitle(paste(tt,"m",sep = ""))+
  facet_wrap(~pred_time,nrow=3,scale="free_y")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        legend.position = "none",
        strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())


pdf("../../good_mt/SAC_follow.pdf",width = 12,height = 8)
# wrap_plots(plot_ls,nrow = 3)
p
dev.off()






#################################################################
#################################################################
library(survival)
library(survminer)
library(parallel)
library(JMbayes2) #,lib.loc = "/hwdata/home/fengzq7/R/4.0/JMbayes2" 
a = data_use %>% dplyr::mutate(Score=log(`Indole 3-carbinol`))
a$measure_time <- as.numeric(a$measure_time)
colnames(a)

lmefit1 <- lme(Score ~ ns(measure_time, 2), data = a,
               random = ~ ns(measure_time, 2) | SUBJID,
               control = lmeControl(opt = 'optim'))
a.id <- a[!duplicated(a$SUBJID), ]
coxfit <- coxph(Surv(AVAL1,CNSR1) ~ 1,
                data = a.id, x = TRUE)

jmod<-jm(coxfit,list(lmefit1),
         time_var ="measure_time")#,


fForms1 <- list("Score" = ~ value(Score))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s1 = data.frame(
  type="value",
  coef=s$Survival["value(Score)","Mean"],p=s$Survival["value(Score)","P"],
  l=s$Survival["value(Score)","2.5%"],h=s$Survival["value(Score)","97.5%"]
)


fForms1 <- list("Score" = ~ area(Score))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s2 = data.frame(
  type="area",
  coef=s$Survival["area(Score)","Mean"],p=s$Survival["area(Score)","P"],
  l=s$Survival["area(Score)","2.5%"],h=s$Survival["area(Score)","97.5%"]
)


fForms1 <- list("Score" = ~ slope(Score))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s3 = data.frame(
  type="slope",
  coef=s$Survival["slope(Score)","Mean"],p=s$Survival["slope(Score)","P"],
  l=s$Survival["slope(Score)","2.5%"],h=s$Survival["slope(Score)","97.5%"]
)

