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
setwd("~/help_for_others/DYQ/output/7_JM_3/")
dir.create("pre_metabo")
setwd("./pre_metabo/")

load("norm.rdata")


setwd("~/help_for_others/DYQ/output/7_JM_3//")
library(tidyverse)
load("./joint_res_1/joint_tb_and_res_1.rdata")
res_1 = res %>% dplyr::select(metabo,coef,p)
colnames(res_1)[2:3]=str_c(colnames(res_1)[2:3],"R1",sep = "_")
res_1$coef_R1 = log(res_1$coef_R1)
load("./joint_res_2/joint_tb_and_res_2.rdata")
res_2 = res %>% dplyr::select(metabo,coef,p)
colnames(res_2)[2:3]=str_c(colnames(res_2)[2:3],"R2",sep = "_")

res = res_1 %>% inner_join(res_2,by = "metabo") %>% na.omit()

cor.test(res$coef_R1,res$coef_R2)



final =res %>% dplyr::mutate(label=ifelse(coef_R1*coef_R2>0,1,0))
final$label = factor(final$label,levels = c(0,1))


load("../1_preclean/data/data_C1_preclean.rdata")
mt_info = data_in_C1 %>% dplyr::select(Compounds,`Class I`,`Class II`) %>% 
  dplyr::rename(mt=Compounds)
# VIP_FC = sel %>% inner_join(mt_info,by="mt")

class_ID_1 = mt_info %>% dplyr::select(mt,`Class I`) %>% dplyr::rename(class_1=colnames(.)[2])
class_ID_2 = mt_info %>% dplyr::select(mt,`Class II`) %>% dplyr::rename(class=colnames(.)[2])

final = final %>% dplyr::rename(mt=metabo) %>% inner_join(class_ID_2,by = "mt") %>% inner_join(class_ID_1,by = "mt")




#######################################################################
#######################################################################
Cer = final %>% dplyr::filter(class %in% c("Cer-NS","Cer-AS","Cer-NP","Cer-AP","HexCer-AP","Cer-NDS")) %>% 
  dplyr::filter(p_R1<0.1 | p_R2<0.1,coef_R1>0 & coef_R2>0) %>% pull(mt)
SM = final %>% dplyr::filter(class %in% c("SM"))  %>% dplyr::filter(p_R1<0.1 | p_R2<0.1,coef_R1>0 & coef_R2>0) %>% pull(mt)
#"Cer-AP","Cer-AS","Cer-NP", 
# SP = "LPI(20:3/0:0)"

SP_norm = norm %>% dplyr::select(1,any_of(Cer),any_of(SM))


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
RRA_SP = RRA %>% dplyr::rename(Sample=Name) %>% 
  inner_join(clinic_in[,c("Sample","SUBJID","cycle","AVAL1","CNSR1","R")],by = "Sample")


RRA_SP$CNSR1 = ifelse(RRA_SP$CNSR1==0,1,0)

load("../../data/blood_time.rdata")
blood_time = blood_time %>% dplyr::mutate(cycle=ifelse(Cycle=="筛选期(D-7~D-1)","C1D1",Cycle))

RRA_SP = RRA_SP %>% inner_join(blood_time[,c("SUBJID","cycle","Time")],by = c("SUBJID","cycle")) %>% 
  dplyr::mutate(measure_time=ifelse(Time<0,0,Time)) %>% dplyr::mutate(measure_time=measure_time/30.44)

RRA_SP = RRA_SP %>% dplyr::rename(SP=Score)







PL = final %>% dplyr::filter(class %in% c("LPE","LPI","LPA","PA","PI","LPC","PC","PE")) %>% 
  dplyr::filter(p_R1<0.1 | p_R2<0.1,coef_R1<0 & coef_R2<0) %>% 
  pull(mt)

# PL = "LPI(20:3/0:0)"

PL_norm = norm %>% dplyr::select(1,any_of(PL))


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


load("../1_preclean/data/data_preclean.rdata")
RRA_PL = RRA %>% dplyr::rename(PL=Score,Sample=Name)

head(RRA_PL)
head(RRA_SP)



RRA = RRA_PL %>% inner_join(RRA_SP,by = "Sample")

save(RRA,file = "./ratio_pred/data_RRA_for_better_2.rdata")
save(Cer,PL,SM,file = "./ratio_pred/all_mt_use_for_model.rdata")







#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
library(survival)
library(survminer)
library(parallel)
library(JMbayes2) #,lib.loc = "/hwdata/home/fengzq7/R/4.0/JMbayes2" 
a = RRA
a$measure_time <- as.numeric(a$measure_time)
# if(length(unique(a$metabo))>5){

# a$start_time <- a$measure_time
df=2

lmefit1 <- lme(SP ~ ns(measure_time, 3), data = a,
               random = ~ ns(measure_time, 3) | SUBJID,
               control = lmeControl(opt = 'optim'))
lmefit2 <- lme(PL ~ ns(measure_time, 3), data = a,
               random = ~ ns(measure_time, 3) | SUBJID,
               control = lmeControl(opt = 'optim'))


# lmefit
# s=summary(lmefit)

a.id <- a[!duplicated(a$SUBJID), ]
coxfit <- coxph(Surv(AVAL1,CNSR1) ~ 1,
                data = a.id, x = TRUE)

# cl <- makePSOCKcluster(1, manual=TRUE, outfile='log.txt')
jmod<-jm(coxfit,list(lmefit1,lmefit2),
         time_var ="measure_time")#,
# n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)
summary(jmod)


fForms1 <- list("SP" = ~ value(SP),"PL" = ~ value(PL))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s1 = data.frame(
  type_SP="value",type_PL="value",DIC=s$fit_stats$conditional$DIC,
  WAIC=s$fit_stats$conditional$WAIC,LPML=s$fit_stats$conditional$WAIC,
  coef_SP=s$Survival["value(SP)","Mean"],p_SP=s$Survival["value(SP)","P"],l_SP=s$Survival["value(SP)","2.5%"],h_SP=s$Survival["value(SP)","97.5%"],
  coef_PL=s$Survival["value(PL)","Mean"],p_PL=s$Survival["value(PL)","P"],l_PL=s$Survival["value(PL)","2.5%"],h_PL=s$Survival["value(PL)","97.5%"],
  long_ns_df=df
  # long_sigma=s$Outcome1["sigma","Mean"],p_long_sigma=s$Outcome1["sigma","P"]
)

model_1 = jmod1


fForms1 <- list("SP" = ~ value(SP),"PL" = ~ slope(PL))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s2 = data.frame(
  type_SP="value",type_PL="slope",DIC=s$fit_stats$conditional$DIC,
  WAIC=s$fit_stats$conditional$WAIC,LPML=s$fit_stats$conditional$WAIC,
  coef_SP=s$Survival["value(SP)","Mean"],p_SP=s$Survival["value(SP)","P"],l_SP=s$Survival["value(SP)","2.5%"],h_SP=s$Survival["value(SP)","97.5%"],
  coef_PL=s$Survival["slope(PL)","Mean"],p_PL=s$Survival["slope(PL)","P"],l_PL=s$Survival["slope(PL)","2.5%"],h_PL=s$Survival["slope(PL)","97.5%"],
  long_ns_df=df
  # long_sigma=s$Outcome1["sigma","Mean"],p_long_sigma=s$Outcome1["sigma","P"]
)
model_2 = jmod1


fForms1 <- list("SP" = ~ value(SP),"PL" = ~ area(PL))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s3 = data.frame(
  type_SP="value",type_PL="area",DIC=s$fit_stats$conditional$DIC,
  WAIC=s$fit_stats$conditional$WAIC,LPML=s$fit_stats$conditional$WAIC,
  coef_SP=s$Survival["value(SP)","Mean"],p_SP=s$Survival["value(SP)","P"],l_SP=s$Survival["value(SP)","2.5%"],h_SP=s$Survival["value(SP)","97.5%"],
  coef_PL=s$Survival["area(PL)","Mean"],p_PL=s$Survival["area(PL)","P"],l_PL=s$Survival["area(PL)","2.5%"],h_PL=s$Survival["area(PL)","97.5%"],
  long_ns_df=df
  # long_sigma=s$Outcome1["sigma","Mean"],p_long_sigma=s$Outcome1["sigma","P"]
)
model_3 = jmod1


fForms1 <- list("SP" = ~ area(SP),"PL" = ~ value(PL))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s4 = data.frame(
  type_SP="area",type_PL="value",DIC=s$fit_stats$conditional$DIC,
  WAIC=s$fit_stats$conditional$WAIC,LPML=s$fit_stats$conditional$WAIC,
  coef_SP=s$Survival["area(SP)","Mean"],p_SP=s$Survival["area(SP)","P"],l_SP=s$Survival["area(SP)","2.5%"],h_SP=s$Survival["area(SP)","97.5%"],
  coef_PL=s$Survival["value(PL)","Mean"],p_PL=s$Survival["value(PL)","P"],l_PL=s$Survival["value(PL)","2.5%"],h_PL=s$Survival["value(PL)","97.5%"],
  long_ns_df=df
  # long_sigma=s$Outcome1["sigma","Mean"],p_long_sigma=s$Outcome1["sigma","P"]
)
model_4 = jmod1


fForms1 <- list("SP" = ~ area(SP),"PL" = ~ slope(PL))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s5 = data.frame(
  type_SP="area",type_PL="slope",DIC=s$fit_stats$conditional$DIC,
  WAIC=s$fit_stats$conditional$WAIC,LPML=s$fit_stats$conditional$WAIC,
  coef_SP=s$Survival["area(SP)","Mean"],p_SP=s$Survival["area(SP)","P"],l_SP=s$Survival["area(SP)","2.5%"],h_SP=s$Survival["area(SP)","97.5%"],
  coef_PL=s$Survival["slope(PL)","Mean"],p_PL=s$Survival["slope(PL)","P"],l_PL=s$Survival["slope(PL)","2.5%"],h_PL=s$Survival["slope(PL)","97.5%"],
  long_ns_df=df
  # long_sigma=s$Outcome1["sigma","Mean"],p_long_sigma=s$Outcome1["sigma","P"]
)
model_5 = jmod1


fForms1 <- list("SP" = ~ area(SP),"PL" = ~ area(PL))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s6 = data.frame(
  type_SP="area",type_PL="area",DIC=s$fit_stats$conditional$DIC,
  WAIC=s$fit_stats$conditional$WAIC,LPML=s$fit_stats$conditional$WAIC,
  coef_SP=s$Survival["area(SP)","Mean"],p_SP=s$Survival["area(SP)","P"],l_SP=s$Survival["area(SP)","2.5%"],h_SP=s$Survival["area(SP)","97.5%"],
  coef_PL=s$Survival["area(PL)","Mean"],p_PL=s$Survival["area(PL)","P"],l_PL=s$Survival["area(PL)","2.5%"],h_PL=s$Survival["area(PL)","97.5%"],
  long_ns_df=df
  # long_sigma=s$Outcome1["sigma","Mean"],p_long_sigma=s$Outcome1["sigma","P"]
)
model_6 = jmod1




fForms1 <- list("SP" = ~ slope(SP),"PL" = ~ value(PL))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s7 = data.frame(
  type_SP="slope",type_PL="value",DIC=s$fit_stats$conditional$DIC,
  WAIC=s$fit_stats$conditional$WAIC,LPML=s$fit_stats$conditional$WAIC,
  coef_SP=s$Survival["slope(SP)","Mean"],p_SP=s$Survival["slope(SP)","P"],l_SP=s$Survival["slope(SP)","2.5%"],h_SP=s$Survival["slope(SP)","97.5%"],
  coef_PL=s$Survival["value(PL)","Mean"],p_PL=s$Survival["value(PL)","P"],l_PL=s$Survival["value(PL)","2.5%"],h_PL=s$Survival["value(PL)","97.5%"],
  long_ns_df=df
  # long_sigma=s$Outcome1["sigma","Mean"],p_long_sigma=s$Outcome1["sigma","P"]
)
model_7 = jmod1


fForms1 <- list("SP" = ~ slope(SP),"PL" = ~ slope(PL))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s8 = data.frame(
  type_SP="slope",type_PL="slope",DIC=s$fit_stats$conditional$DIC,
  WAIC=s$fit_stats$conditional$WAIC,LPML=s$fit_stats$conditional$WAIC,
  coef_SP=s$Survival["slope(SP)","Mean"],p_SP=s$Survival["slope(SP)","P"],l_SP=s$Survival["slope(SP)","2.5%"],h_SP=s$Survival["slope(SP)","97.5%"],
  coef_PL=s$Survival["slope(PL)","Mean"],p_PL=s$Survival["slope(PL)","P"],l_PL=s$Survival["slope(PL)","2.5%"],h_PL=s$Survival["slope(PL)","97.5%"],
  long_ns_df=df
  # long_sigma=s$Outcome1["sigma","Mean"],p_long_sigma=s$Outcome1["sigma","P"]
)
model_8 = jmod1




fForms1 <- list("SP" = ~ slope(SP),"PL" = ~ area(PL))
jmod1 <- update(jmod, functional_forms = fForms1)
s = summary(jmod1)

res_s9 = data.frame(
  type_SP="slope",type_PL="area",DIC=s$fit_stats$conditional$DIC,
  WAIC=s$fit_stats$conditional$WAIC,LPML=s$fit_stats$conditional$LPML,
  coef_SP=s$Survival["slope(SP)","Mean"],p_SP=s$Survival["slope(SP)","P"],l_SP=s$Survival["slope(SP)","2.5%"],h_SP=s$Survival["slope(SP)","97.5%"],
  coef_PL=s$Survival["area(PL)","Mean"],p_PL=s$Survival["area(PL)","P"],l_PL=s$Survival["area(PL)","2.5%"],h_PL=s$Survival["area(PL)","97.5%"],
  long_ns_df=df
  # long_sigma=s$Outcome1["sigma","Mean"],p_long_sigma=s$Outcome1["sigma","P"]
)
model_9 = jmod1





res_s = rbind(res_s1,res_s2,res_s3,res_s4,res_s5,res_s6,res_s7,res_s8,res_s9)
model_ls = list(model_1,model_2,model_3,model_4,model_5,model_6,model_7,model_8,model_9)

save(res_s,file = "./ratio_pred/area_slope_area_try_res_for_better.rdata")
save(model_ls,file = "./ratio_pred/area_slope_value_model_for_better.rdata")




######################################################
######################################################
load("./ratio_pred/final_JM_model.rdata")
ND = RRA %>% dplyr::filter(SUBJID=="05020",measure_time<9) %>% arrange(measure_time)
ND$AVAL1=9
ND$CNSR1=0
predSurv <- predict(jmod1, newdata = ND, process = "event",
                    times = seq(9, 35, length.out = 51),
                    return_newdata = TRUE)
predLong <- predict(jmod1, newdata = ND,
                    times = seq(9, 35, length.out = 51),
                    return_newdata = TRUE)

library(randomcoloR)
cols <- randomColor(2)
colnames(predLong$newdata2)

source("~/help_for_others/DYQ/code/7_JM_3/source_script//source_predict_plot.R")
pdf("./ratio_pred/JM_model.pdf",width = 5,height = 5)
plot_modified(x=predLong, x2=predSurv, outcomes = 1:2, subject = 1,
              # fun_long = list(exp, identity, identity),
              fun_event = function (x) 1 - x,
              ylab_event = "Survival Probabilities",
              ylab_long = colnames(predLong$newdata2)[2:3],
              col_points = cols, col_line_long = cols,
              fill_CI_long = c("#F25C7880", "#D973B580", "#F2832280"),
              fill_CI_event = "#F7F7FF80",
              pos_ylab_long = c(1.9, 1.9, 0.08))
dev.off()
