#!/usr/bin/Rscript
rm(list = ls())
setwd("~/help_for_others/DYQ/output/3_baseline_model//")

load("~/help_for_others/DYQ/output/3_baseline_model/csv/df_validation_norm.rdata")

load("./mt_block/final_clone_0.8.rdata")


clone_ls = unique(final$clone.id)

kk=1
res_each_simi_1 = lapply(1:length(clone_ls),function(kk){
  sel_mt = final %>% dplyr::filter(clone.id==clone_ls[kk]) %>% pull(mt)
  tmp = validation %>% dplyr::select(any_of(sel_mt))
  cor_tb = cor(tmp)
  return(data.frame(clone_No=clone_ls[kk],right_cor=table(cor_tb>=0.7)["TRUE"]-length(sel_mt),
                    total=length(sel_mt)*(length(sel_mt)-1)))
  
}) %>% do.call(rbind,.)

res_each_simi_1$high_cor_prob = res_each_simi_1$right_cor/res_each_simi_1$total




load("~/help_for_others/DYQ/output/2_response_mt_enricher//csv/df_chemo_norm.rdata")
kk=1
res_each_simi_2 = lapply(1:length(clone_ls),function(kk){
  sel_mt = final %>% dplyr::filter(clone.id==clone_ls[kk]) %>% pull(mt)
  tmp = chemo %>% dplyr::select(any_of(sel_mt))
  cor_tb = cor(tmp)
  return(data.frame(clone_No=clone_ls[kk],right_cor=table(cor_tb>=0.7)["TRUE"]-length(sel_mt),
                    total=length(sel_mt)*(length(sel_mt)-1)))
  
}) %>% do.call(rbind,.)

res_each_simi_2$high_cor_prob = res_each_simi_2$right_cor/res_each_simi_2$total




load("~/help_for_others/DYQ/output/2_response_mt_enricher//csv/df_surg_norm.rdata")
kk=1
res_each_simi_3 = lapply(1:length(clone_ls),function(kk){
  sel_mt = final %>% dplyr::filter(clone.id==clone_ls[kk]) %>% pull(mt)
  tmp = surg %>% dplyr::select(any_of(sel_mt))
  cor_tb = cor(tmp)
  return(data.frame(clone_No=clone_ls[kk],right_cor=table(cor_tb>=0.7)["TRUE"]-length(sel_mt),
                    total=length(sel_mt)*(length(sel_mt)-1)))
  
}) %>% do.call(rbind,.)

res_each_simi_3$high_cor_prob = res_each_simi_3$right_cor/res_each_simi_3$total





load("~/help_for_others/DYQ/output/2_response_mt_enricher//csv/df_rad_norm.rdata")
kk=1
res_each_simi_4 = lapply(1:length(clone_ls),function(kk){
  sel_mt = final %>% dplyr::filter(clone.id==clone_ls[kk]) %>% pull(mt)
  tmp = rad %>% dplyr::select(any_of(sel_mt))
  cor_tb = cor(tmp)
  return(data.frame(clone_No=clone_ls[kk],right_cor=table(cor_tb>=0.7)["TRUE"]-length(sel_mt),
                    total=length(sel_mt)*(length(sel_mt)-1)))
  
}) %>% do.call(rbind,.)

res_each_simi_4$high_cor_prob = res_each_simi_4$right_cor/res_each_simi_4$total


table(res_each_simi_1$high_cor_prob==1)
table(res_each_simi_2$high_cor_prob==1)
table(res_each_simi_3$high_cor_prob==1)
table(res_each_simi_4$high_cor_prob==1)


res1=res_each_simi_1 %>% dplyr::select(clone_No,high_cor_prob) %>% dplyr::rename(val=high_cor_prob)
res2=res_each_simi_2 %>% dplyr::select(clone_No,high_cor_prob) %>% dplyr::rename(chemo=high_cor_prob)
res3=res_each_simi_3 %>% dplyr::select(clone_No,high_cor_prob) %>% dplyr::rename(surg=high_cor_prob)
res4=res_each_simi_4 %>% dplyr::select(clone_No,high_cor_prob) %>% dplyr::rename(rad=high_cor_prob)
colnames(final)[1]="clone_No"
merge = final %>% inner_join(res1,by = "clone_No") %>% inner_join(res2,by = "clone_No") %>% 
  inner_join(res3,by = "clone_No") %>% inner_join(res4,by = "clone_No")


merge_plot = merge %>% dplyr::select(1,3:6) %>% pivot_longer(cols = !clone_No,names_to = "cohort",values_to = "high_cor_prob")
merge_plot$cohort = factor(merge_plot$cohort,levels = rev(c("val","chemo","surg","rad")),ordered = TRUE)
library(ggridges)
library(wesanderson)
pdf("./mt_block/mt_block_validation.pdf",width = 5,height = 4)
ggplot(merge_plot)+
  geom_density_ridges(aes(x=high_cor_prob,y=cohort,fill=cohort,color = cohort))+
  scale_color_manual(values = c("#8dd2c5","#bfbcda","#f47f72","#7fb2d5"))+
  scale_fill_manual(values = c("#8dd2c5","#bfbcda","#f47f72","#7fb2d5"))+
  theme_bw()
dev.off()
