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
library(tidyverse)
library(stringr)

dir.create("~/help_for_others/DYQ/output/3_baseline_model/training_single_mt/")  #######
setwd("~/help_for_others/DYQ/output/3_baseline_model/training_single_mt/")
library(MetaboAnalystR)
library(tidyverse)
mSet<-InitDataObjects("conc", "stat", FALSE)


mSet<-Read.TextData(mSet, "~/help_for_others/DYQ/output/3_baseline_model/csv/C1D1_response_df_training.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet)  ####
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "S10T0", ratio=FALSE, ratioNum=20)

mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###

training = cbind(mSet$dataSet$cls,mSet$dataSet$norm) %>% dplyr::rename(response=colnames(.)[1])
save(training,file = "../csv/df_training_norm.rdata")

#########################
#########################
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

# Plot of the T-test results
mSet<-PlotTT(mSet, imgName = "tt_0_", format = "pdf", dpi = 72, width=NA)
t_tb = mSet$analSet$tt$p.value %>% data.frame() %>% rownames_to_column("mt")
colnames(t_tb)[2]="p"


FC_tb = FC_tb %>% inner_join(t_tb,by = "mt")


# #####################Volcano Plot
# Perform the volcano analysis
mSet<-Volcano.Anal(mSet, paired=FALSE, fcthresh = 2.0, cmpType = 0, nonpar = F, threshp = 0.1,
                   equal.var = TRUE, pval.type = "fdr")
# Create the volcano plot
# mSet<-PlotVolcano(mSet, "volcano_0_", 1, 0, format ="pdf", dpi=72, width=NA)

#
#
####################Orthogonal Partial Least Squares - Discriminant Analysis (orthoPLS-DA)
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


VIP_tb <- mSet$analSet$oplsda$vip.mat %>% as.data.frame(check.names=FALSE) %>% rownames_to_column("mt") %>% dplyr::select(-V2)



sel <- FC_tb %>% inner_join(VIP_tb,by = "mt")
ggplot(sel)+
  geom_point(aes(x=V1, y = FC))

save(sel,file = "./training_sel.rdata")





#############################################################validation
#############################################################validation
#############################################################validation
dir.create("~/help_for_others/DYQ/output/3_baseline_model/validation_single_mt/")  #######
setwd("~/help_for_others/DYQ/output/3_baseline_model/validation_single_mt/")
library(MetaboAnalystR)
library(tidyverse)
mSet<-InitDataObjects("conc", "stat", FALSE)


mSet<-Read.TextData(mSet, "~/help_for_others/DYQ/output/3_baseline_model/csv/C1D1_response_df_validation.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet)  ####
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "S10T0", ratio=FALSE, ratioNum=20)

mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###

validation = cbind(mSet$dataSet$cls,mSet$dataSet$norm) %>% dplyr::rename(response=colnames(.)[1])
save(validation,file = "../csv/df_validation_norm.rdata")


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

# Plot of the T-test results
mSet<-PlotTT(mSet, imgName = "tt_0_", format = "pdf", dpi = 72, width=NA)
t_tb = mSet$analSet$tt$p.value %>% data.frame() %>% rownames_to_column("mt")
colnames(t_tb)[2]="p"


FC_tb = FC_tb %>% inner_join(t_tb,by = "mt")


# #####################Volcano Plot
# Perform the volcano analysis
mSet<-Volcano.Anal(mSet, paired=FALSE, fcthresh = 2.0, cmpType = 0, nonpar = F, threshp = 0.1,
                   equal.var = TRUE, pval.type = "fdr")
# Create the volcano plot
# mSet<-PlotVolcano(mSet, "volcano_0_", 1, 0, format ="pdf", dpi=72, width=NA)

#
#
####################Orthogonal Partial Least Squares - Discriminant Analysis (orthoPLS-DA)
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


VIP_tb <- mSet$analSet$oplsda$vip.mat %>% as.data.frame(check.names=FALSE) %>% rownames_to_column("mt") %>% dplyr::select(-V2)







