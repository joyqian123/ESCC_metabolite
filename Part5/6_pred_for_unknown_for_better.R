#!/usr/bin/Rscript
rm(list = ls())
setwd("~/help_for_others/DYQ/output/")
library(stringr)
library(gridExtra)
library(future)
# library(sva)
# library(e1071)
library(pROC)
library(ROCit)
library(caret)
library(doParallel)
library(cancerclass)
library(tidyverse)
library(JMbayes2)
setwd("~/help_for_others/DYQ/output/7_JM_3/")

load("./ratio_pred/all_mt_use_for_model.rdata")



###################
setwd("~/help_for_others/DYQ/output/7_JM_3/")
dir.create("pre_metabo")
setwd("./pre_metabo/")

load("norm.rdata")
ref = norm



load("../../1_preclean/validation_all_follow_up/data_preclean.rdata")
sample_all = clinic_in$Sample
data_in_for_csv = data_in %>% dplyr::select(Compounds,any_of(sample_all)) %>% 
  column_to_rownames("Compounds") %>% t() %>% data.frame(.,check.names = FALSE) %>% 
  rownames_to_column("Sample") %>% inner_join(clinic_in[,c("Sample","R")],by = "Sample") %>% 
  dplyr::select(Sample,R,everything())
table(duplicated(colnames(data_in)))
table(duplicated(clinic_in$Sample))
table(duplicated(data_in_for_csv$Sample))

write.csv(data_in_for_csv,file = "~/help_for_others/DYQ/output/7_JM_3/csv/unknown_df.csv",row.names = FALSE)



##########################################
setwd("~/help_for_others/DYQ/output/7_JM_3/")
dir.create("pre_metabo_unknown")
setwd("./pre_metabo_unknown//")
# write.csv(N_df,"all_response_df.csv",row.names = FALSE)

if(T){
  library(MetaboAnalystR)
  mSet<-InitDataObjects("conc", "stat", FALSE)
  mSet<-Read.TextData(mSet, "../csv/unknown_df.csv", "rowu", "disc")
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet)  ####
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "NULL", "LogNorm", "AutoNorm", "S10T0", ratio=FALSE, ratioNum=20)
  norm=mSet$dataSet$norm %>% rownames_to_column("Sample")  
  save(norm,file = "norm.rdata")
}


load("norm.rdata")
unknown = norm


tt=15

setwd("../")
load("./ratio_pred/area_slope_value_model_for_better.rdata")

mm_nn=3
for(mm_nn in c(3)){
  

  all_SUBJID_pred_res = lapply(1:length(unique(clinic_in$SUBJID)),function(tt){
    
    print(tt)
    load("../1_preclean/validation_all_follow_up/data_preclean.rdata")
    subj_id_in = unique(clinic_in$SUBJID)[tt]
    sample_in = clinic_in %>% dplyr::filter(SUBJID==subj_id_in) %>% pull(Sample)
    norm_for_RRA = unknown %>% dplyr::filter(Sample %in% sample_in) %>% rbind(ref,.)
    clinic_for_RRA = clinic_in %>% dplyr::filter(Sample %in% sample_in) %>% dplyr::select(Sample,SUBJID,PFS,status,response,blood_time) %>% 
      dplyr::rename(AVAL1=PFS,CNSR1=status,R=response) %>% dplyr::mutate(measure_time=ifelse(blood_time<0,0,blood_time)) %>% 
      dplyr::select(-blood_time) %>% dplyr::mutate(cycle=str_c("c",1:nrow(.))) %>% 
      dplyr::mutate(measure_time=measure_time/30.44)
    
    #######################################################################################
    SP_norm = norm_for_RRA %>% dplyr::select(1,any_of(Cer),any_of(SM))
    
    library(RobustRankAggreg)
    method_choose="median"
    glist<-list()
    i=2
    for (i in 2:ncol(SP_norm)){
      data<-SP_norm %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
      glist[[i-1]]<-data$Sample
      names(glist)[i-1]<-colnames(SP_norm)[i]
    }
    library(RobustRankAggreg)
    r<-rankMatrix(glist,full = TRUE)
    # r<-ScoreMatrix(glist,full = TRUE)
    RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
    
    load("..//1_preclean/data/data_preclean.rdata")
    clinic_in$CNSR1 = ifelse(clinic_in$CNSR1==0,1,0)
    clinic_in = rbind(clinic_in[,c("Sample","SUBJID","AVAL1","CNSR1","R","cycle")],
                      clinic_for_RRA[,c("Sample","SUBJID","AVAL1","CNSR1","R","cycle")])
    RRA_SP = RRA %>% dplyr::rename(Sample=Name) %>% 
      inner_join(clinic_in,by = "Sample")
    
    
    load("../../data/blood_time.rdata")
    blood_time = blood_time %>% dplyr::mutate(cycle=ifelse(Cycle=="筛选期(D-7~D-1)","C1D1",Cycle)) %>% 
      dplyr::select(SUBJID,cycle,Time) %>% 
      dplyr::mutate(measure_time=ifelse(Time<0,0,Time)) %>% dplyr::mutate(measure_time=measure_time/30.44) %>% 
      dplyr::select(-Time) %>% distinct(.keep_all = T)
    blood_time_unknown = clinic_for_RRA %>% dplyr::select(SUBJID,cycle,measure_time)
    
    blood_time = rbind(blood_time,blood_time_unknown)
    
    RRA_SP = RRA_SP %>% inner_join(blood_time[,c("SUBJID","cycle","measure_time")],by = c("SUBJID","cycle")) 
    RRA_SP = RRA_SP %>% dplyr::rename(SP=Score)
    
    
    
    
    #######################################################################################
    PL_norm = norm_for_RRA %>% dplyr::select(1,any_of(PL))
    
    library(RobustRankAggreg)
    method_choose="median"
    glist<-list()
    i=2
    for (i in 2:ncol(PL_norm)){
      data<-PL_norm %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
      glist[[i-1]]<-data$Sample
      names(glist)[i-1]<-colnames(PL_norm)[i]
    }
    library(RobustRankAggreg)
    r<-rankMatrix(glist,full = TRUE)
    # r<-ScoreMatrix(glist,full = TRUE)
    RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
    
    load("..//1_preclean/data/data_preclean.rdata")
    clinic_in$CNSR1 = ifelse(clinic_in$CNSR1==0,1,0)
    clinic_in = rbind(clinic_in[,c("Sample","SUBJID","AVAL1","CNSR1","R","cycle")],
                      clinic_for_RRA[,c("Sample","SUBJID","AVAL1","CNSR1","R","cycle")])
    RRA_PL = RRA %>% dplyr::rename(Sample=Name) %>% 
      inner_join(clinic_in,by = "Sample")
    
    
    load("../../data/blood_time.rdata")
    blood_time = blood_time %>% dplyr::mutate(cycle=ifelse(Cycle=="筛选期(D-7~D-1)","C1D1",Cycle)) %>% 
      dplyr::select(SUBJID,cycle,Time) %>% 
      dplyr::mutate(measure_time=ifelse(Time<0,0,Time)) %>% dplyr::mutate(measure_time=measure_time/30.44) %>% 
      dplyr::select(-Time) %>% distinct(.keep_all = T)
    blood_time_unknown = clinic_for_RRA %>% dplyr::select(SUBJID,cycle,measure_time)
    
    blood_time = rbind(blood_time,blood_time_unknown)
    
    RRA_PL = RRA_PL %>% inner_join(blood_time[,c("SUBJID","cycle","measure_time")],by = c("SUBJID","cycle")) 
    RRA_PL = RRA_PL %>% dplyr::rename(PL=Score)
    
    RRA = RRA_PL %>% inner_join(RRA_SP,by = c("SUBJID","Sample","cycle","measure_time","AVAL1","CNSR1","R"))
    
    
    
    
    print(mm_nn)
    jmod1 = model_ls[[mm_nn]]

    pred_tb = RRA %>% dplyr::select(SUBJID,PL,SP,measure_time) %>% dplyr::filter(SUBJID==subj_id_in)
    class(pred_tb$measure_time)
    
    time_ls = pred_tb %>% arrange(measure_time) %>% pull(measure_time) %>% unique() %>% .[-1]+0.01 %>% as.numeric() %>% sort()
    
    source("../../code/7_JM_3/source_script/pred_combined.R")
    library(JMbayes2)
    # library(JM)
    pred_res_all = lapply(time_ls,function(tt){
      t0=tt
      pred_res <- predict_combined_original(jmod1, newdata2 = pred_tb, 
                                            process = "event",
                                            Tstart = t0,
                                            times = seq(t0,36,0.2),return_newdata = TRUE)
      return(pred_res)
    })
    
    res_pred2 = RRA %>% dplyr::select(SUBJID,AVAL1,CNSR1) %>% distinct(SUBJID,.keep_all = T)
    
    
    
    ########################################################  
    library(parallel)
    mc = getOption("mc.cores",5)
    tune_matrix_1 = expand.grid(target_cif=seq(0.1,0.9,0.1),cut_off_delta=seq(1,24,1))
    
    pred_tb_ref = RRA %>% dplyr::select(SUBJID,PL,SP,measure_time) %>% dplyr::filter(!SUBJID==subj_id_in)
    pred_res_ref = lapply(time_ls,function(tt){
      t0=tt
      pred_res <- predict_combined_original(jmod1, newdata2 = pred_tb_ref, 
                                            process = "event",
                                            Tstart = t0,
                                            times = seq(t0,36,0.2),return_newdata = TRUE)
      return(pred_res)
    })
    
    plot_large_res = lapply(1:length(pred_res_ref),function(uu){
      
      data_on1 = pred_res_ref[[uu]]
      res_pred1 = RRA %>% dplyr::select(SUBJID,AVAL1,CNSR1) %>% distinct(SUBJID,.keep_all = T) %>% 
        dplyr::filter(SUBJID %in% data_on1$SUBJID)
      t0=time_ls[uu]
      
      CVdats <- create_folds(data_on1, V = 5, id_var = "SUBJID",seed = 123)
      
      #########################
      CV_res =lapply(1:5,function(vv){
        
        print(str_c("vv=",vv,sep=""))
        
        data_on = CVdats$training[[vv]]
        data_on1 = data_on %>% dplyr::select(SUBJID,AVAL1,pred_CIF) %>% 
          pivot_wider(id_cols = AVAL1,names_from = SUBJID,values_from = pred_CIF) %>% 
          column_to_rownames("AVAL1")
        
        id = unique(data_on$SUBJID)
        
        
        #################
        c_res_1 = mclapply(1:nrow(tune_matrix_1),mc.cores = mc,function(kk){
          
          print(str_c("kk=",kk,sep=""))
          target_cif = tune_matrix_1$target_cif[kk]
          cut_off_delta = tune_matrix_1$cut_off_delta[kk]
          
          
          res_0.5 = lapply(1:length(id),function(i){
            tmp = data_on %>% dplyr::filter(SUBJID %in% id[i])
            
            time_at_cif_0.5 <- approx(x = tmp$pred_CIF, y = tmp$AVAL1, xout = target_cif)$y
            return(data.frame(SUBJID=id[i],time_at_CIF_0.5=time_at_cif_0.5))
          }) %>% do.call(rbind,.)
          
          
          res_pred = res_pred1 %>% 
            inner_join(res_0.5,by = "SUBJID") %>% 
            dplyr::filter(AVAL1 >= t0) %>% 
            dplyr::mutate(
              pred_group = ifelse(is.na(time_at_CIF_0.5) | time_at_CIF_0.5>cut_off_delta+t0,"low_risk","high_risk"))
          
          if(length(table(res_pred$pred_group))==1){
            return(data.frame(pred_time=t0,co=paste0("training",vv),c_index=NA,target=target_cif,cut_off=cut_off_delta))
          }else{
            cox<-coxph(Surv(AVAL1, CNSR1) ~ pred_group, data = res_pred)
            c=summary(cox)$concordance[1]
            return(data.frame(pred_time=t0,co=paste0("training",vv),c_index=c,target=target_cif,cut_off=cut_off_delta))
          }
        }) %>% do.call(rbind,.)
        
      }) %>% do.call(rbind,.)
      
    }) %>% do.call(rbind,.)
    
    
    
    
    library(survival)
    library(survminer)
    CV_ana = plot_large_res %>% #dplyr::filter(pred_time>=4) %>% 
      pivot_wider(id_cols = c(pred_time,target,cut_off),names_from = "co",values_from="c_index")  %>% 
      dplyr::mutate(mean_c = rowMeans(.[,c(4:8)]))
    
    try_1 = CV_ana %>% #dplyr::filter(pred_time>=4) %>%
      group_by(target,cut_off) %>% dplyr::summarize(mean_mean_c=mean(mean_c,na.rm=T))
    
    tune_res = CV_ana %>% group_by(pred_time) %>% dplyr::reframe(mean_c=max(mean_c,na.rm = T)) %>% ungroup() %>% 
      inner_join(CV_ana,by = c("pred_time","mean_c")) %>% distinct(pred_time,.keep_all = T)
    
    
    
    
    i=1
    res = lapply(1:nrow(tune_res),function(i){
      
      target_cif <- tune_res$target[i]
      cut_off_delta  = tune_res$cut_off[i]
      
      data_on = pred_res_all[[i]]
      t0 = tune_res$pred_time[i]
      
      data_on1 = data_on %>% dplyr::select(SUBJID,AVAL1,pred_CIF) %>%
        pivot_wider(id_cols = AVAL1,names_from = SUBJID,values_from = pred_CIF) %>%
        column_to_rownames("AVAL1")
      
      id = unique(data_on$SUBJID)
      
      i=1
      res_0.5 = lapply(1:length(id),function(i){
        tmp = data_on %>% dplyr::filter(SUBJID %in% id[i])
        
        time_at_cif_0.5 <- approx(x = tmp$pred_CIF, y = tmp$AVAL1, xout = target_cif)$y
        return(data.frame(SUBJID=id[i],time_at_CIF_0.5=time_at_cif_0.5))
      }) %>% do.call(rbind,.)
      
      res_pred = res_pred2 %>% 
        inner_join(res_0.5,by = "SUBJID") %>%
        # dplyr::filter(AVAL1 >= t0) %>%
        dplyr::mutate(
          pred_group = ifelse(is.na(time_at_CIF_0.5) | time_at_CIF_0.5>cut_off_delta+t0,"low_risk","high_risk")) %>% 
        dplyr::mutate(pred_time = t0,cut_delta=cut_off_delta,cut=cut_off_delta+t0)
      
    }) %>% do.call(rbind,.)
    
    return(res)
    
  }) %>% do.call(rbind,.)
  
  
  
  dir.create(paste0("./ratio_pred/Model_",mm_nn,"/pred_unknown/"))
  save(all_SUBJID_pred_res,file = paste0("./ratio_pred/Model_",mm_nn,"/pred_unknown/all_SUBJID_pred_res_specific_time.rdata"))
  
  

  load(paste0("./ratio_pred/Model_",mm_nn,"/pred_unknown/all_SUBJID_pred_res_specific_time.rdata"))
  


}







#########################################
plot_ls = lapply(1:12,function(dd){
  res_pred_tmp = all_SUBJID_pred_res %>% dplyr::mutate(if_dead_in_delta=ifelse(AVAL1<pred_time+dd & CNSR1==1,1,
                                                                          ifelse(AVAL1>pred_time+dd,0,NA))) %>%
    na.omit()
  
  res_pred_tmp$if_dead_in_delta = factor(res_pred_tmp$if_dead_in_delta,levels = c(0,1))
  res_pred_tmp$pred_group = ifelse(res_pred_tmp$pred_group=="high_risk",1,0)
  
  roc = pROC::roc(res_pred_tmp$if_dead_in_delta,res_pred_tmp$pred_group,levels = c(0,1),direction="<")
  
  pred_data <- res_pred_tmp %>%
    count(pred_group, if_dead_in_delta) %>%
    group_by(pred_group) %>%
    mutate(
      total = sum(n),
      proportion = n / total,
      percentage = round(proportion * 100, 1)
    )
  
  # 计算总体CR比例
  overall_cr_rate <- sum(pred_data$if_dead_in_delta == 1) / nrow(pred_data)
  
  # 计算每个预测组的CR比例
  cr_rates <- pred_data %>%
    group_by(pred_group) %>%
    summarise(
      cr_rate = mean(if_dead_in_delta == 1),
      cr_percentage = round(mean(if_dead_in_delta == 1) * 100, 1),
      n = n(),
      .groups = 'drop'
    )
  
  # 绘制堆砌柱状图
  p <- ggplot(pred_data, aes(x = pred_group, y = proportion, fill = if_dead_in_delta)) +
    geom_col(width = 0.7, color = "white", size = 0.8) +
    geom_text(aes(label = paste0(percentage, "%")),
              position = position_stack(vjust = 0.5),
              color = "black", fontweight = "bold", size = 4) +
    scale_fill_manual(values = c("#b6cada","#f4edd3")) +
    # scale_y_continuous(labels = scales::percent_format(accuracy = 1),
    #                    expand = c(0, 0, 0.05, 0)) +
    # labs(
    #   x = "tVAT area",
    #   y = "Proportion") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  return(p)
  
})

pdf(paste0("./ratio_pred/Model_",mm_nn,"/pred_unknown/pred_data_progression_risk.pdf"),width = 24,height = 6)
wrap_plots(plot_ls,nrow=2)
dev.off()
