rm(list = ls())
setwd("~/help_for_others/DYQ/output/6_long_survival_model/")  #######
# .libPaths()
library(MetaboAnalystR)
library(tidyverse)

####################################################
####################################################
load("./training_C1C6_delta_metabo/training_C1C6_delta_single_metabo_FC_tb.rdata")
delta_single_FC = sel
sel_mt = delta_single_FC %>% dplyr::filter(p<0.3) %>%
  pull(mt)

final_sel_mt = sel_mt


load("./csv/training_validation_df.rdata")
training = training %>% dplyr::select(1,2,any_of(final_sel_mt))
validation = validation %>% dplyr::select(1,2,any_of(final_sel_mt))
all=rbind(training,validation)

write.csv(training,file = "~/help_for_others/DYQ/output/6_long_survival_model/csv/training_C1C6_delta_after_select.csv",row.names = FALSE)
write.csv(validation,file = "~/help_for_others/DYQ/output/6_long_survival_model/csv/validation_C1C6_delta_after_select.csv",row.names = FALSE)
write.csv(all,file = "~/help_for_others/DYQ/output/6_long_survival_model/csv/all_C1C6_delta_after_select.csv",row.names = FALSE)

save(training,validation,file = "./csv/after_select_metabo_training_validation.rdata")










rm(list = ls())
setwd("~/help_for_others/DYQ/output/6_long_survival_model/")  #######在最开始设置工作目录
dir.create("./training_C1C6_delta_after_selection")
setwd("./training_C1C6_delta_after_selection/")
# .libPaths()
library(MetaboAnalystR)
library(tidyverse)

# .libPaths()
# library(MetaboAnalystR)
# library(tidyverse)
mSet<-InitDataObjects("conc", "stat", FALSE)


mSet<-Read.TextData(mSet, "~/help_for_others/DYQ/output/6_long_survival_model/csv/training_C1C6_delta_after_select.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)

# a0 <- mSet$dataSet$norm


####运行但不进行replacemin
# mSet<-ReplaceMin(mSet)  ####替换0值或缺失值
mSet$dataSet$filt <- mSet$dataSet$edit <- NULL
preproc <- qs::qread("preproc.qs")
int.mat <- preproc
mSet$dataSet$proc.feat.num <- ncol(int.mat)
qs::qsave(as.data.frame(int.mat), file = "data_proc.qs")
mSet$msgSet$replace.msg <- paste("Zero or missing values were replaced by 1/5 of the min positive value for each variable.")
invisible(gc())


###运行
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL","CrNorm","MeanCenter", ratio=FALSE,ratioNum = 20)
mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###按代谢物作图
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###按样本作图


a <- mSet$dataSet$norm


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
####第二个：use a non-parametric test（是否使用非参检验）
####第三个：the adjusted p-value (FDR) cutoff
####第四个：配对
####最后：all_results	Logical, if TRUE, returns T-Test analysis results for all compounds
t_tb = data.frame(t_score=mSet$analSet$tt$t.score,p=mSet$analSet$tt$p.value)%>% data.frame() %>% rownames_to_column("mt")
FC_tb =  FC_tb %>% inner_join(t_tb,by = "mt")

# #####################Volcano Plot
# Perform the volcano analysis
# mSet<-Volcano.Anal(mSet, paired=FALSE, fcthresh = 2.0, cmpType = 0, nonpar = F, threshp = 0.1,
#                    equal.var = TRUE, pval.type = "fdr")
# Create the volcano plot
# mSet<-PlotVolcano(mSet, "volcano_0_", 1, 0, format ="pdf", dpi=72, width=NA)

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





VIP_tb <- mSet$analSet$oplsda$vip.mat %>% as.data.frame(check.names=FALSE) %>%
  rownames_to_column("mt") %>% dplyr::select(-V2)


sel <- FC_tb %>% inner_join(VIP_tb,by = "mt")
ggplot(sel)+
  geom_point(aes(x=V1, y = FC))

save(VIP_tb,sel,file = "./after_select_metabo_training_VIP_sel.rdata")


norm_data = mSet$dataSet$norm
training = cbind(response=mSet$dataSet$cls,norm_data)

save(training,file = "../csv/CrNorm_normalized_training_data.rdata")

#######################################feature_selection
load("./after_select_metabo_training_VIP_sel.rdata")
intersect_mt = VIP_tb$mt
load("../csv/after_select_metabo_training_validation.rdata")


load("../../8_compare_with_mt_change_in_R_and_NR/VD/mt_list_from_data.rdata")


sel = sel %>% dplyr::filter(abs(FC)>1,V1>1.5,p<0.05) 

pdf("~/help_for_others/DYQ/output/6_long_survival_model_2/try.pdf",width = 5,height = 5)
sel_plot = sel %>% arrange(desc(V1)) %>% dplyr::mutate(rank=seq(1,nrow(.),1))
ggplot(sel_plot,aes(x=rank,y=V1))+
  geom_point()+
  geom_line()
dev.off()



i=4
t=2
sig_ls <- lapply(2:10,function(i){
  print(i)
  for(t in 2:i){
    if(t==2){
      sel_sig_1 = sel %>% arrange(desc(V1)) %>% .[1:t,] %>% pull(mt)
      training_try <- training[,colnames(training) %in% c(sel_sig_1)]
      cor_sel=cor(training_try)
      if(table(cor_sel>0.9)["TRUE"]==ncol(cor_sel)){
        sel_sig = sel_sig_1
      }else{
        sel_sig = sel_sig_1[-t]
      }
    }else{
      sel_sig_add <- sel %>% arrange(desc(V1)) %>% .[t,] %>% pull(mt)
      sel_sig_1 <- c(sel_sig,sel_sig_add) %>% unlist()
      training_try <- training[,colnames(training) %in% c(sel_sig_1)]
      cor_sel=cor(training_try)
      print(table(cor_sel>0.9)["TRUE"]==ncol(cor_sel))
      if(table(cor_sel>0.9)["TRUE"]==ncol(cor_sel)){
        sel_sig = sel_sig_1
      }else{
        sel_sig = sel_sig
      }
    }
  }
  return(sel_sig)
})

save(sig_ls,file = "~/help_for_others/DYQ/output/6_long_survival_model_2/training_C1C6_delta_after_selection/training_feature_selection_based_on_p.rdata")





