#!/usr/bin/Rscript
rm(list = ls())
library(stringr)
library(gridExtra)
library(future)
library(sva)
# library(e1071)
library(pROC)
library(ROCit)
library(caret)
library(doParallel)
library(cancerclass)

file.copy("~/help_for_others/DYQ/output/2_response_mt_enricher/csv/df_chemo_norm.rdata",
          to = "~/help_for_others/DYQ/output/3_baseline_model/csv/df_chemo_norm.rdata",overwrite = TRUE)
file.copy("~/help_for_others/DYQ/output/2_response_mt_enricher/csv/df_surg_norm.rdata",
          to = "~/help_for_others/DYQ/output/3_baseline_model/csv/df_surg_norm.rdata",overwrite = TRUE)
file.copy("~/help_for_others/DYQ/output/2_response_mt_enricher/csv/df_rad_norm.rdata",
          to = "~/help_for_others/DYQ/output/3_baseline_model/csv/df_rad_norm.rdata",overwrite = TRUE)


dir.create("~/help_for_others/DYQ/output/3_baseline_model/mt_block_model_res/")  #######在最开始设置工作目录
setwd("~/help_for_others/DYQ/output/3_baseline_model/mt_block_model_res/")
load("../mt_block/final_clone_0.8.rdata")


load("../csv/df_training_norm.rdata")

kk=1
block_dt = lapply(1:length(unique(final$clone.id)),function(kk){
  sel_mt = final %>% dplyr::filter(clone.id==unique(final$clone.id)[kk]) %>% pull(mt)
  tmp = training %>% dplyr::select(any_of(sel_mt)) %>% dplyr::mutate(value=rowMeans(.)) %>% 
    dplyr::select(value)
  colnames(tmp)=str_c("bl",unique(final$clone.id)[kk],sep = "_")
  return(tmp)
}) %>% do.call(cbind,.)

block_dt_train = data.frame(response=training$response,block_dt) %>% rownames_to_column("Sample")

write.csv(block_dt_train,file = "../csv/block_training.csv",row.names = FALSE)



load("../csv/df_validation_norm.rdata")

kk=1
block_dt = lapply(1:length(unique(final$clone.id)),function(kk){
  sel_mt = final %>% dplyr::filter(clone.id==unique(final$clone.id)[kk]) %>% pull(mt)
  tmp = validation %>% dplyr::select(any_of(sel_mt)) %>% dplyr::mutate(value=rowMeans(.)) %>% 
    dplyr::select(value)
  colnames(tmp)=str_c("bl",unique(final$clone.id)[kk],sep = "_")
  return(tmp)
}) %>% do.call(cbind,.)

block_dt_validation = data.frame(response=validation$response,block_dt) %>% rownames_to_column("Sample")

write.csv(block_dt_validation,file = "../csv/block_validation.csv",row.names = FALSE)




load("../csv/df_chemo_norm.rdata")


kk=1
block_dt = lapply(1:length(unique(final$clone.id)),function(kk){
  sel_mt = final %>% dplyr::filter(clone.id==unique(final$clone.id)[kk]) %>% pull(mt)
  tmp = chemo %>% dplyr::select(any_of(sel_mt)) %>% dplyr::mutate(value=rowMeans(.)) %>% 
    dplyr::select(value)
  colnames(tmp)=str_c("bl",unique(final$clone.id)[kk],sep = "_")
  return(tmp)
}) %>% do.call(cbind,.)

block_dt_chemo = data.frame(response=chemo$response,block_dt) %>% rownames_to_column("Sample")

write.csv(block_dt_chemo,file = "../csv/block_chemo.csv",row.names = FALSE)





load("../csv/df_surg_norm.rdata")

kk=1
block_dt = lapply(1:length(unique(final$clone.id)),function(kk){
  sel_mt = final %>% dplyr::filter(clone.id==unique(final$clone.id)[kk]) %>% pull(mt)
  tmp = surg %>% dplyr::select(any_of(sel_mt)) %>% dplyr::mutate(value=rowMeans(.)) %>% 
    dplyr::select(value)
  colnames(tmp)=str_c("bl",unique(final$clone.id)[kk],sep = "_")
  return(tmp)
}) %>% do.call(cbind,.)

block_dt_surg = data.frame(response=surg$response,block_dt) %>% rownames_to_column("Sample")

write.csv(block_dt_surg,file = "../csv/block_surg.csv",row.names = FALSE)




load("../csv/df_rad_norm.rdata")

kk=1
block_dt = lapply(1:length(unique(final$clone.id)),function(kk){
  sel_mt = final %>% dplyr::filter(clone.id==unique(final$clone.id)[kk]) %>% pull(mt)
  tmp = rad %>% dplyr::select(any_of(sel_mt)) %>% dplyr::mutate(value=rowMeans(.)) %>% 
    dplyr::select(value)
  colnames(tmp)=str_c("bl",unique(final$clone.id)[kk],sep = "_")
  return(tmp)
}) %>% do.call(cbind,.)

block_dt_rad = data.frame(response=rad$response,block_dt) %>% rownames_to_column("Sample")

write.csv(block_dt_rad,file = "../csv/block_rad.csv",row.names = FALSE)

