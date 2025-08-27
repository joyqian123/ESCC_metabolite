rm(list = ls())
setwd("~/help_for_others/DYQ/output/5_delta_metabo_response_14m_discovering/metabo_process/") 
dir.create("./C1C6_delta")
setwd("./C1C6_delta/")
# .libPaths()
library(MetaboAnalystR)
library(tidyverse)

# .libPaths()
library(MetaboAnalystR)
library(tidyverse)
mSet<-InitDataObjects("conc", "stat", FALSE)


mSet<-Read.TextData(mSet, "~/help_for_others/DYQ/output/5_delta_metabo_response_14m_discovering/csv/C1C6_delta_df.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)

mSet$dataSet$filt <- mSet$dataSet$edit <- NULL
preproc <- qs::qread("preproc.qs")
int.mat <- preproc
mSet$dataSet$proc.feat.num <- ncol(int.mat)
qs::qsave(as.data.frame(int.mat), file = "data_proc.qs")
mSet$msgSet$replace.msg <- paste("Zero or missing values were replaced by 1/5 of the min positive value for each variable.")
invisible(gc())



mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL","CrNorm","MeanCenter", ratio=FALSE,ratioNum = 20)
mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)  



#####################fold change analysis
# Perform fold-change analysis on uploaded data, unpaired
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)  
mSet<-PlotFC(mSet, "fc_0_", "pdf", 72, width=NA)

# To view fold-change 
mSet$analSet$fc$fc.log

FC_tb = mSet$analSet$fc$fc.log %>% data.frame() %>% rownames_to_column("mt")
colnames(FC_tb)[2]="FC"

#####################T-Test
# Perform T-test (parametric)
mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, pvalType = "fdr", TRUE)
t_tb = data.frame(t_score=mSet$analSet$tt$t.score,p=mSet$analSet$tt$p.value)%>% data.frame() %>% rownames_to_column("mt")
FC_tb =  FC_tb %>% inner_join(t_tb,by = "mt")



###################Orthogonal Partial Least Squares - Discriminant Analysis (orthoPLS-DA)
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






save(FC_tb,file = "../../C1C6_delta/C1C6_delta_metabo_FC_tb.rdata")




########################################################################
#######################################################################
#######################################################################
sel_bl = FC_tb %>% dplyr::filter(p<0.3) %>%
  pull(mt)
norm = mSet$dataSet$norm %>% dplyr::select(any_of(sel_bl))
block_dt_train_sel = data.frame(response=mSet$dataSet$cls,norm) %>% rownames_to_column("Sample")

write.csv(block_dt_train_sel,file = "../../csv/C1C6_delta_df_after_selection.csv",row.names = FALSE)



setwd("~/help_for_others/DYQ/output/5_delta_metabo_response_14m_discovering/metabo_process/")  
dir.create("./C1C6_delta_after_selection")
setwd("./C1C6_delta_after_selection//")
library(MetaboAnalystR)
library(tidyverse)
mSet<-InitDataObjects("conc", "stat", FALSE)


mSet<-Read.TextData(mSet, "~/help_for_others/DYQ/output/5_delta_metabo_response_14m_discovering/csv/C1C6_delta_df_after_selection.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)


mSet$dataSet$filt <- mSet$dataSet$edit <- NULL
preproc <- qs::qread("preproc.qs")
int.mat <- preproc
mSet$dataSet$proc.feat.num <- ncol(int.mat)
qs::qsave(as.data.frame(int.mat), file = "data_proc.qs")
mSet$msgSet$replace.msg <- paste("Zero or missing values were replaced by 1/5 of the min positive value for each variable.")
invisible(gc())



mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL","CrNorm","MeanCenter", ratio=FALSE,ratioNum = 20)
mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)  





#####################fold change analysis
# Perform fold-change analysis on uploaded data, unpaired
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)  
mSet<-PlotFC(mSet, "fc_0_", "pdf", 72, width=NA)

# To view fold-change 
mSet$analSet$fc$fc.log

FC_tb = mSet$analSet$fc$fc.log %>% data.frame() %>% rownames_to_column("mt")
colnames(FC_tb)[2]="FC"

#####################T-Test
# Perform T-test (parametric)
mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, pvalType = "fdr", TRUE)
####参数：
t_tb = data.frame(t_score=mSet$analSet$tt$t.score,p=mSet$analSet$tt$p.value)%>% data.frame() %>% rownames_to_column("mt")
FC_tb =  FC_tb %>% inner_join(t_tb,by = "mt")

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