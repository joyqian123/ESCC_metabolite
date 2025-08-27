rm(list = ls())
setwd("~/help_for_others/DYQ/")  #######
library(MetaboAnalystR)

library(tidyverse)



load("./output/1_preclean/data/data_preclean.rdata")
mt_info = data_in[,c(1:3)]

mt_class = table(mt_info$`Class I`) %>% data.frame() %>% dplyr::filter(Freq>0)
mt_class$Var1 = as.character(mt_class$Var1)
mt_class = mt_class %>% dplyr::mutate(class_final = c(rep("others",2),mt_class$Var1[3],rep("others",2),
                                                      mt_class$Var1[6],rep("others",1),rep("lipid",3),
                                                      rep("others",2),mt_class$Var1[13:14],rep("others",2),
                                                      rep("lipid",2),"others"
                                                      ))
mt_class_final = mt_class %>% group_by(class_final) %>% dplyr::summarize(sum=sum(Freq))
library(ggplot2)
mt_class_final$percentage <- round(mt_class_final$sum / sum(mt_class_final$sum), 2)
mt_class_final$label <- paste0(mt_class_final$class_final, "\n", mt_class_final$sum, " (", mt_class_final$percentage, "%)")

mt_class_final$class_final = factor(mt_class_final$class_final,levels = mt_class_final$class_final,ordered = TRUE)
mt_class_final <- mt_class_final %>%
  arrange(desc(class_final)) %>%
  mutate(
    cumulative = cumsum(percentage),
    midpoint = (cumulative - 0.5*sum/sum(mt_class_final$sum))* sum(mt_class_final$sum)    # 
  )


# 
pdf("./output/2_response_mt_enricher/plot/mt_class.pdf",width = 10,height = 5)
ggplot(mt_class_final, aes(x = "", y = sum, fill = class_final)) +
  geom_bar(stat = "identity", width = 1) +  
  coord_polar(theta = "y") +  
  theme_void() +  
  labs(title = "1384 metabolites") +
  # geom_segment(aes(
  #   x = 1.2, xend = 1.4,
  #   y = cumulative - sum, yend = midpoint
  # ), color = "black", size = 0.5) +
  geom_text(aes(x=1.5,y=midpoint,
    label = label
  ), hjust = 0, size = 3)+
  scale_fill_manual(values = c("#cdacb4","#e9d2c8", "#f4efcb","#d6e1d3","#b3cbe2","#8fabda"))+
  theme(legend.position = "none")
dev.off()







##############################################################
##############################################################
setwd("./output/2_response_mt_enricher/data_enr/")
# Create mSetObj
mSet<-InitDataObjects("conc", "msetora", FALSE)
mSet<-Read.TextData(mSet, "../../1_preclean/data/C1_df.csv", "rowu", "disc");

mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet)  ####
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "S10T0", ratio=FALSE, ratioNum=20)

mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###


norm = mSet$dataSet$norm
group = mSet$dataSet$cls


hvg = apply(norm,2,var) %>% data.frame() %>% dplyr::rename(var=colnames(.)[1]) %>% 
  arrange(desc(var))

target_cf = c(100,200,300,400,500)

for (i in 1:length(target_cf)){
  
  hvg_mt = hvg %>% rownames(.) %>% .[1:target_cf[i]]
  norm_use = norm %>% dplyr::select(any_of(hvg_mt)) 
  norm_use = cbind(data.frame(response=group),norm_use) %>% rownames_to_column("Sample")
  
  
  write.csv(norm_use,file = str_c("~/help_for_others/DYQ/output/2_response_mt_enricher/csv/hvg_df",target_cf[i],".csv",sep = ""),
            row.names = FALSE)
  
  
  setwd("../")
  dir.create(str_c("hvg",target_cf[i],sep = ""))
  setwd(str_c("./hvg",target_cf[i],sep = ""))
  # Create mSetObj
  mSet<-InitDataObjects("conc", "msetora", FALSE)
  mSet<-Read.TextData(mSet, str_c("~/help_for_others/DYQ/output/2_response_mt_enricher/csv/hvg_df",target_cf[i],".csv",sep = ""), 
                      "rowu", "disc");
  
  mSet<-SanityCheckData(mSet)
  
  mSet<-ReplaceMin(mSet)  ####
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "NULL", "NULL", "NULL", "S10T0", ratio=FALSE, ratioNum=20)
  
  mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###
  mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###
  
  
  #####################fold change analysis
  # Perform fold-change analysis on uploaded data, unpaired
  mSet<-FC.Anal(mSet, 2.0, 0, FALSE)  

  # Plot fold-change analysis
  mSet<-PlotFC(mSet, "fc_0_", "pdf", 72, width=NA)
  
  # To view fold-change 
  mSet$analSet$fc$fc.log
  
  FC_tb = mSet$analSet$fc$fc.log %>% data.frame() %>% rownames_to_column("mt")
  colnames(FC_tb)[2]="FC"
  
  
  #####################T-Test
  # Perform T-test (parametric)
  mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, pvalType = "fdr", TRUE)
  t_tb = mSet$analSet$tt$p.value%>% data.frame() %>% rownames_to_column("mt")
  FC_tb =  FC_tb %>% inner_join(t_tb,by = "mt")
  
  colnames(FC_tb)[3]="p"
  FC_tb$p.adj = p.adjust(FC_tb$p,method = "fdr")
  
  
  #####################T-Test
  # Perform T-test (parametric)
  mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, pvalType = "fdr", TRUE)
  t_tb = mSet$analSet$tt$p.value%>% data.frame() %>% rownames_to_column("mt")
  FC_tb =  FC_tb %>% inner_join(t_tb,by = "mt")
  
  # Perform oPLS-DA analysis
  mSet<-OPLSR.Anal(mSet, reg=TRUE)
  
  # Create a 2D oPLS-DA score plot
  mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", format = "pdf", dpi=72, width=NA, inx1 = 1,inx2 = 2,
                        reg = 0.95,show = 1,grey.scale = 0)
  mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", plotType = "all", format = "pdf", dpi = 72, width=NA)
  mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", format = "pdf", dpi = 72, width=NA, type = "vip", feat.nm = "tscore", 
                     feat.num = 15,color.BW = FALSE)
  mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", format = "pdf", dpi=72, width=NA)
  mSet<-OPLSDA.Permut(mSet, 100)
  mSet<-PlotOPLS.Permutation(mSet, "opls_perm_2_", format = "pdf", dpi=72, width=NA)

}






