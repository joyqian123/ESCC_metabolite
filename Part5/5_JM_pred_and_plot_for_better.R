#!/usr/bin/Rscript
rm(list = ls())

not_run=FALSE

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

load("./ratio_pred/data_RRA_for_better_2.rdata")


load("./ratio_pred/area_slope_value_model_for_better.rdata")

pred_tb = RRA[,c(4,2,3,10)]


source("../../code/7_JM_3/source_script/pred_combined.R")


mm_nn=3
for(mm_nn in 1:length(model_ls)){
  
  print(mm_nn)
  jmod1 = model_ls[[mm_nn]]
  
  if(not_run){
    pred_res_all = lapply(1:25,function(tt){
      t0=tt
      pred_res <- predict_combined_original(jmod1, newdata2 = pred_tb, 
                                            process = "event",
                                            Tstart = t0,
                                            times = seq(t0,36,0.2),return_newdata = TRUE)
      return(pred_res)
    })
    
    dir.create(paste0("./ratio_pred/Model_",mm_nn))
    save(pred_res_all,file = paste0("./ratio_pred/Model_",mm_nn,"/pred_res_all.rdata"))
  }
  
  
  load(paste0("./ratio_pred/Model_",mm_nn,"/pred_res_all.rdata"))
  
  
  # kk=2
  
  if(T){
    all_res = lapply(1:25,function(kk){
    
      data_on1 = pred_res_all[[kk]]
      res_pred = RRA %>% dplyr::select(SUBJID,AVAL1,CNSR1) %>% distinct(SUBJID,.keep_all = T)
      t0=kk
      # vv=2
    
      res = lapply(seq(t0+1,t0+34,1),function(vv){
        pred_time = vv
        tmp = data_on1 %>% dplyr::filter(AVAL1==pred_time)
        res_pred_tmp = res_pred %>% dplyr::mutate(if_dead_in_delta=ifelse(AVAL1<pred_time & CNSR1==1,1,
                                                                          ifelse(AVAL1>pred_time,0,NA))) %>%
          dplyr::filter(SUBJID %in% tmp$SUBJID) %>% na.omit()
    
        if(length(table(res_pred_tmp$if_dead_in_delta))>1){
          #   return(NULL)
          # }else{
          auc_tb = res_pred_tmp %>% inner_join(tmp[,c(1,4)],by = "SUBJID")
          auc_tb$pred_CIF = as.numeric(auc_tb$pred_CIF)
          auc_tb$if_dead_in_delta = factor(auc_tb$if_dead_in_delta,levels = c(0,1))
          roc = pROC::roc(auc_tb$if_dead_in_delta,auc_tb$pred_CIF,levels = c(0,1),direction="<")
    
          roc_result <- coords(roc, "best")
          auc <- roc %>% auc()
          auc_ci <- ci(roc)
          
          df_validation <- data.frame(t0=t0,
                                      pred_time=pred_time,
                                      ROC=auc_ci[2], ROC_l=auc_ci[1],ROC_u=auc_ci[3])
          return(df_validation)
        }
      }) %>% do.call(rbind,.)
    
    }) %>% do.call(rbind,.)
    
    

    
    

    ################################################################
    ##################################################################
    tune_matrix_1 = expand.grid(target_cif=seq(0.1,0.9,0.1),cut_off_delta=seq(1,24,1))
    # tune_matrix_2 = expand.grid(target_delta=seq(1,12,1),cut_off_cif=seq(0.1,0.8,0.1))
    
    library(parallel)
    mc<-getOption("mc.cores",12)
    
    
    plot_large_res = lapply(1:25,function(uu){
      print(str_c("uu=",uu,sep=""))
      # load(paste("./plot/JM/tune/output_seq12_wo_chemo_rm_add_coef_rm_0.3_hr/pred_res/",model_ls[uu],sep = ""))
      
      
      data_on1 = pred_res_all[[uu]]
      res_pred1 = RRA %>% dplyr::select(SUBJID,AVAL1,CNSR1) %>% distinct(SUBJID,.keep_all = T) %>% 
        dplyr::filter(SUBJID %in% data_on1$SUBJID)
      t0=uu
    
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
            return(data.frame(pred_time=uu,co=paste0("training",vv),c_index=NA,target=target_cif,cut_off=cut_off_delta,
                              hr_l = NA ,p_l = NA))
          }else{
            cox<-coxph(Surv(AVAL1, CNSR1) ~ pred_group, data = res_pred)
            c=summary(cox)$concordance[1]
            hr_l = summary(cox)$coef[2]
            p_l = summary(cox)$coef[5]
            return(data.frame(pred_time=uu,co=paste0("training",vv),c_index=c,target=target_cif,cut_off=cut_off_delta,
                              hr_l = hr_l ,p_l = p_l))
          }
        }) %>% do.call(rbind,.)
        
      }) %>% do.call(rbind,.)
      
    }) %>% do.call(rbind,.) 
    
    
    save(plot_large_res,file = paste0("./ratio_pred/Model_",mm_nn,"/CV5_pred_cutoff.rdata"))
    
    
  }
  
  
  
  library(survival)
  library(survminer)
  load(paste0("./ratio_pred/Model_",mm_nn,"/CV5_pred_cutoff.rdata"))
  CV_ana = plot_large_res %>% #dplyr::filter(hr_l<1) %>% 
    pivot_wider(id_cols = c(pred_time,target,cut_off),names_from = "co",values_from="c_index")  %>% 
    dplyr::mutate(mean_c = rowMeans(.[,c(4:8)]))
  
  
  try_1 = CV_ana %>% #dplyr::filter(pred_time>=4) %>%
    group_by(target,cut_off) %>% dplyr::summarize(mean_mean_c=mean(mean_c,na.rm=T))
  
  
  
  tune_res = CV_ana %>% group_by(pred_time) %>% dplyr::reframe(mean_c=max(mean_c,na.rm = T)) %>% ungroup() %>% 
    inner_join(CV_ana,by = c("pred_time","mean_c")) %>% distinct(pred_time,.keep_all = T)

  
  res_pred1 = RRA %>% dplyr::select(SUBJID,AVAL1,CNSR1) %>% distinct(SUBJID,.keep_all = T)
  
  if(not_run){
    
    plot_large_res = lapply(1:nrow(tune_res),function(i){
      
      target_cif <- tune_res$target[i]
      cut_off_delta  = tune_res$cut_off[i]
      
      data_on = pred_res_all[[i]]
      t0 = tune_res$pred_time[i]
      
      data_on1 = data_on %>% dplyr::select(SUBJID,AVAL1,pred_CIF) %>%
        pivot_wider(id_cols = AVAL1,names_from = SUBJID,values_from = pred_CIF) %>%
        column_to_rownames("AVAL1")
      
      
      library(ComplexHeatmap)
      library(circlize)
      library(dendextend)
      library(ggsci)
      
      id = unique(data_on$SUBJID)
      
      i=1
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
        return(NULL)
      }else{
        fit<-survfit(Surv(AVAL1, CNSR1) ~ pred_group, data = res_pred)
        summary(fit)
        
        cox<-coxph(Surv(AVAL1, CNSR1) ~ pred_group, data = res_pred)
        c=summary(cox)$concordance[1]
        
        p=ggsurvplot(fit, data = res_pred,
                     surv.median.line = "hv", 
                     conf.int = FALSE, 
                     risk.table = TRUE, 
                     pval = TRUE,
                     palette = "jco",  
                     xlab = "Follow up time(m)", 
                     legend = c(0.8,0.9), 
                     title = paste("pred time:",round(t0,2),sep = ""), 
                     legend.labs = c("High risk", "Low risk"), 
                     xlim=c(0,40),
                     break.x.by = 10)  
        p1 = p$plot + geom_vline(xintercept = t0,lty=2)+
          geom_text(x=5,y=0.1,label=paste0("C_index:",round(c,3),sep=""))
        return(p1)
      }
    })
    
    plot_large_res <- Filter(Negate(is.null), plot_large_res)
    
    
    
    library(patchwork)
    pdf(paste0("./ratio_pred/Model_",mm_nn,"/pred_res_plot_best_1.pdf"),
        width = 6,height = 4)
    lapply(1:length(plot_large_res),function(kk){
      print(wrap_plots(plot_large_res[[kk]],ncol = 1))
    })
    dev.off()
    
    
    
  }
  
  
  
  
  plot_large_res = lapply(1:nrow(tune_res),function(i){

    target_cif <- tune_res$target[i]
    cut_off_delta  = tune_res$cut_off[i]

    data_on = pred_res_all[[i]]
    t0 = tune_res$pred_time[i]

    data_on1 = data_on %>% dplyr::select(SUBJID,AVAL1,pred_CIF) %>%
      pivot_wider(id_cols = AVAL1,names_from = SUBJID,values_from = pred_CIF) %>%
      column_to_rownames("AVAL1")


    library(ComplexHeatmap)
    library(circlize)
    library(dendextend)
    library(ggsci)

    id = unique(data_on$SUBJID)

    i=1
    res_0.5 = lapply(1:length(id),function(i){
      tmp = data_on %>% dplyr::filter(SUBJID %in% id[i])

      time_at_cif_0.5 <- approx(x = tmp$pred_CIF, y = tmp$AVAL1, xout = target_cif)$y
      return(data.frame(SUBJID=id[i],time_at_CIF_0.5=time_at_cif_0.5))
    }) %>% do.call(rbind,.)


    res_pred = res_pred1 %>%
      inner_join(res_0.5,by = "SUBJID") %>%
      dplyr::filter(AVAL1 >= t0) %>%
      dplyr::mutate(
        pred_group = ifelse(is.na(time_at_CIF_0.5) | time_at_CIF_0.5>cut_off_delta+t0,"low_risk","high_risk")) %>%
      dplyr::mutate(pred_time=t0)

    return(res_pred)
  }) %>% do.call(rbind,.)


  plot_large_res = plot_large_res %>% dplyr::filter(pred_time<AVAL1)

  all_res = lapply(1:25,function(dd){

    res_pred_tmp = plot_large_res %>% dplyr::mutate(if_dead_in_delta=ifelse(AVAL1<pred_time+dd & CNSR1==1,1,
                                                                                 ifelse(AVAL1>pred_time+dd,0,NA))) %>%
      na.omit()

    res_pred_tmp$if_dead_in_delta = factor(res_pred_tmp$if_dead_in_delta,levels = c(0,1))
    res_pred_tmp$pred_group = ifelse(res_pred_tmp$pred_group=="high_risk",1,0)

    roc = pROC::roc(res_pred_tmp$if_dead_in_delta,res_pred_tmp$pred_group,levels = c(0,1),direction="<")


    roc_result <- pROC::coords(roc, "best")
    auc <- roc %>% auc()
    auc_ci <- ci(roc)
    a = table(res_pred_tmp$if_dead_in_delta==res_pred_tmp$pred_group)
    accur = a["TRUE"]/(a["FALSE"]+a["TRUE"])
    df_validation <- data.frame(delta=dd,accur=accur,
                                # pred_time=pred_time,
                                ROC=auc_ci[2], ROC_l=auc_ci[1],ROC_u=auc_ci[3])
    return(df_validation)
    # }
  }) %>% do.call(rbind,.)


  save(all_res,file = paste0("./ratio_pred/Model_",mm_nn,"/pred_data_delta_res.rdata"))


  load(paste0("./ratio_pred/Model_",mm_nn,"/pred_data_delta_res.rdata"))




  #########################################
  plot_ls = lapply(1:12,function(dd){
    res_pred_tmp = plot_large_res %>% dplyr::mutate(if_dead_in_delta=ifelse(AVAL1<pred_time+dd & CNSR1==1,1,
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

  pdf(paste0("./ratio_pred/Model_",mm_nn,"/pred_data_progression_risk.pdf"),width = 24,height = 6)
  wrap_plots(plot_ls,nrow=2)
  dev.off()

  
  all_res = lapply(1:12,function(kk){

    data_on1 = plot_large_res %>% dplyr::filter(pred_time==kk)

    # res_pred = RRA %>% dplyr::select(SUBJID,AVAL1,CNSR1) %>% distinct(SUBJID,.keep_all = T)
    t0=kk
    # vv=2

    res = lapply(seq(t0+1,t0+12,1),function(vv){
      end_time = vv
      tmp = data_on1 %>% dplyr::filter(AVAL1>=t0)
      res_pred_tmp = tmp %>% dplyr::mutate(if_dead_in_delta=ifelse(AVAL1<end_time & CNSR1==1,1,
                                                                        ifelse(AVAL1>end_time,0,NA))) %>%
        na.omit()

      if(length(table(res_pred_tmp$if_dead_in_delta))>1){
        res_pred_tmp$if_dead_in_delta = factor(res_pred_tmp$if_dead_in_delta,levels = c(0,1))
        res_pred_tmp$pred_group = ifelse(res_pred_tmp$pred_group=="high_risk",1,0)

        roc = pROC::roc(res_pred_tmp$if_dead_in_delta,res_pred_tmp$pred_group,levels = c(0,1),direction="<")

        roc_result <- coords(roc, "best")
        auc <- roc %>% auc()
        auc_ci <- ci(roc)

        df_validation <- data.frame(t0=t0,
                                    pred_time=end_time,
                                    ROC=auc_ci[2], ROC_l=auc_ci[1],ROC_u=auc_ci[3])
        return(df_validation)
      }
    }) %>% do.call(rbind,.)

  }) %>% do.call(rbind,.)


  all_res = all_res %>% dplyr::mutate(delta=pred_time-t0) %>% dplyr::filter(delta<=12)



  
  plot_large_res = lapply(1:nrow(tune_res),function(i){
    
    target_cif <- tune_res$target[i]
    cut_off_delta  = tune_res$cut_off[i]
    
    data_on = pred_res_all[[i]]
    t0 = tune_res$pred_time[i]
    
    data_on1 = data_on %>% dplyr::select(SUBJID,AVAL1,pred_CIF) %>%
      pivot_wider(id_cols = AVAL1,names_from = SUBJID,values_from = pred_CIF) %>%
      column_to_rownames("AVAL1")
    
    
    library(ComplexHeatmap)
    library(circlize)
    library(dendextend)
    library(ggsci)
    
    id = unique(data_on$SUBJID)
    
    i=1
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
      return(NULL)
    }else{
      fit<-survfit(Surv(AVAL1, CNSR1) ~ pred_group, data = res_pred)
      summary(fit)
      
      cox<-coxph(Surv(AVAL1, CNSR1) ~ pred_group, data = res_pred)
      c=summary(cox)$concordance[1]
      hr=summary(cox)$coefficient[1]
      # se=summary(cox)$coefficient[3]
      l_hr = summary(cox)$conf.int[3]
      h_hr = summary(cox)$conf.int[4]
      z=summary(cox)$coefficient[4]
      p=summary(cox)$coefficient[5]
      
      return(data.frame(model=mm_nn,pred_time=t0,cindex=c,hr=hr,l_hr=l_hr,h_hr=h_hr,z_hr=z,p=p))
    }
  }) %>% do.call(rbind,.)
  
  
  save(plot_large_res,file = paste0("./ratio_pred/Model_",mm_nn,"/pred_data_cindex.rdata"))
  

}






res_compare = lapply(1:9,function(mm_nn){
  
  load(paste0("./ratio_pred/Model_",mm_nn,"/pred_data_cindex.rdata"))
  all_res = plot_large_res %>% dplyr::mutate(model=as.character(model)) 
  return(all_res)  
}) %>% do.call(rbind,.)

res_compare = res_compare %>% dplyr::filter(pred_time<=12)


cindex_matrix = res_compare %>% dplyr::select(model,pred_time,cindex) %>% 
  pivot_wider(id_cols = "model",names_from = "pred_time",values_from = "cindex")

library(RobustRankAggreg)
method_choose="median"
glist<-list()
i=2
for (i in 2:ncol(cindex_matrix)){
  data<-cindex_matrix %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(pred_time=colnames(.)[2]) %>% arrange(desc(pred_time))
  glist[[i-1]]<-data$model
  names(glist)[i-1]<-colnames(cindex_matrix)[i]
}
library(RobustRankAggreg)
library(ComplexHeatmap)
library(circlize)
r<-rankMatrix(glist,full = TRUE)

RRA<-aggregateRanks(glist,rmat=r,
                    method=method_choose,full=TRUE,exact=F,topCutoff=NA)

RRA
# row_ha = rowAnnotation(mean_cindex = anno_barplot(rowMeans(cindex_matrix[,-1])),gp = gpar(fill = c("#C71585","#7B68EE","#B0E0E6","#E1FFFF")))
col_fun = colorRamp2(c(0.50, 0.56, 0.66), c("navy", "white", "red"))
# cindex_matrix = cindex_matrix 

pdf("./ratio_pred/cindex_model_heatmap.pdf",width = 6,height = 4)
Heatmap(cindex_matrix[c(1,2,3,7,8,9,4,5,6),-1],cluster_rows = FALSE,cluster_columns = FALSE,
        border = TRUE,col = col_fun,
        rect_gp = gpar(col = "white", lwd = 2),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", cindex_matrix[c(1,2,3,7,8,9,4,5,6),-1][i, j]), x, y, gp = gpar(fontsize = 8))
        }
)
dev.off()




mean_cindex = data.frame(SP=rep(c("value","area","slope"),each=3),GP=rep(c("value","slope","area"),3),mean_index = rowMeans(cindex_matrix[,-1])) %>% 
  pivot_wider(id_cols = SP,names_from = "GP",values_from = "mean_index")
mean_cindex = mean_cindex[c(1,3,2),]
col_fun = colorRamp2(c(0.53, 0.59, 0.6), c("navy", "white", "red"))
pdf("./ratio_pred/cindex_model_heatmap_compare.pdf",width = 2,height =2)
Heatmap(mean_cindex[,-1],cluster_rows = FALSE,cluster_columns = FALSE,
        border = TRUE,col = col_fun,
        rect_gp = gpar(col = "white", lwd = 2),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", mean_cindex[,-1][i, j]), x, y, gp = gpar(fontsize = 8))
        }
)
dev.off()







