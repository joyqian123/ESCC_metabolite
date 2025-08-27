rm(list = ls())
load("~/help_for_others/DYQ/output/1_preclean/data/data_R_preclean.rdata")

colnames(clinic_in_R)
clinic_need = clinic_in_R %>% dplyr::select(1,2,3,5,7:21,25,26)
table(clinic_need$CNSR1)
clinic_need$CNSR1 = ifelse(clinic_need$CNSR1==0,1,0)
clinic_need$is_long_PFS = ifelse(clinic_need$AVAL1<=14,ifelse(clinic_need$CNSR1==0,"unknown","No"),"Yes")
table(clinic_need$is_long_PFS)
clinic_need = clinic_need %>% dplyr::filter(!is_long_PFS=="unknown")




#########
data_C1D1_long = data_C1D1 %>% pivot_longer(cols = !c(1:2),names_to = "mt",values_to = "C1D1") %>% 
  dplyr::mutate(sample=str_remove(Sample,"C1D1"))
data_C6D1_long = data_C6D1 %>% pivot_longer(cols = !c(1:2),names_to = "mt",values_to = "C6D1") %>% 
  dplyr::mutate(sample=str_remove(Sample,"C6D1"))
sample_sel = intersect(data_C1D1_long$sample,data_C6D1_long$sample)
sample_sel_1 = c(str_c(sample_sel,"C1D1"),str_c(sample_sel,"C6D1"))
data_C1C6 = data_in_R %>% dplyr::select(1,any_of(sample_sel_1)) %>% column_to_rownames("Compounds") %>% 
  t() %>% data.frame(.,check.names = FALSE) %>% rownames_to_column("Sample") %>% 
  dplyr::inner_join(clinic_need[,c("Sample","is_long_PFS")],by = "Sample") %>% 
  dplyr::select(Sample,is_long_PFS,everything())
write.csv(data_C1C6,file = "~/help_for_others/DYQ/output/5_delta_metabo_response_14m_discovering/csv/C1C6_pair_df.csv",row.names = FALSE)

save(data_C1C6,clinic_need,file = "~/help_for_others/DYQ/output/5_delta_metabo_response_14m_discovering/csv/14m_paired_sample.rdata")



