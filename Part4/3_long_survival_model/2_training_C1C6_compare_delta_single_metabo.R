# .libPaths()
rm(list = ls())
setwd("~/help_for_others/DYQ/output/6_long_survival_model/")  #######
dir.create("./training_C1C6_delta_metabo")
setwd("./training_C1C6_delta_metabo/")

library(MetaboAnalystR)
library(tidyverse)

mSet<-InitDataObjects("conc", "stat", FALSE)
mSet<-Read.TextData(mSet, "~/help_for_others/DYQ/output/6_long_survival_model/csv/training_delta_df.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)

# a0 <- mSet$dataSet$norm



mSet$dataSet$filt <- mSet$dataSet$edit <- NULL
preproc <- qs::qread("preproc.qs")
int.mat <- preproc
mSet$dataSet$proc.feat.num <- ncol(int.mat)
qs::qsave(as.data.frame(int.mat), file = "data_proc.qs")
mSet$msgSet$replace.msg <- paste("Zero or missing values were replaced by 1/5 of the min positive value for each variable.")
invisible(gc())


###
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL","CrNorm","MeanCenter", ratio=FALSE,ratioNum = 20)
mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###




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

# #####################Volcano Plot

###################Orthogonal Partial Least Squares - Discriminant Analysis (orthoPLS-DA)
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



VIP_tb <- mSet$analSet$oplsda$vip.mat %>% as.data.frame(check.names=FALSE) %>%
  rownames_to_column("mt") %>% dplyr::select(-V2)


sel <- FC_tb %>% inner_join(VIP_tb,by = "mt")
ggplot(sel)+
  geom_point(aes(x=V1, y = FC))



save(VIP_tb,sel,file = "./training_C1C6_delta_single_metabo_FC_tb.rdata")


