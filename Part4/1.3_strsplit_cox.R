rm(list = ls())
setwd("~/help_for_others/DYQ/output/1_preclean/data//")
library(tidyverse)
load("data_preclean.rdata")


clinic_3 = clinic_in %>% distinct(SUBJID,.keep_all = T)

library(survminer) # 加载包
library(survival) # 加载包
table(clinic_3$SUSMFREQ)
clinic_3$SUSMFREQ <- factor(clinic_3$SUSMFREQ,levels = c("从不吸烟","已经戒烟","正在吸烟"),ordered = T)
table(clinic_3$SUDRFREQ)
clinic_3$SUDRFREQ <- factor(clinic_3$SUDRFREQ,levels = c("从不饮酒","已经戒酒","正在饮酒"),ordered = T)

# PFS_tb <- clinic_trt %>% dplyr::select(AVAL1,CNSR1,AVAL4)
clinic_3$CNSR1 <- sapply(1:nrow(clinic_3),function(uu){
     if(clinic_3$CNSR1[uu]==0){
         return(1)
     }else{
           return(0)
     }
})





data_in = clinic_3 %>% dplyr::select(AVAL1,CNSR1)

fit<-survfit(Surv(AVAL1, CNSR1)~1, data = data_in)
ggsurvplot(fit, data = data_in)
ggsurvplot(fit, data = data_in,fun = "cumhaz")
# hazard <- -log(survfit(Surv(AVAL3, CNSR3) ~ 1, data = data_in)$surv)

times <- seq(1,40,1)
sur_tb = data.frame(time=summary(fit,times)$time,sur_prob=summary(fit,times)$surv)
sur_tb$cumHR = 1-sur_tb$sur_prob

ggplot(sur_tb)+
  geom_line(aes(x=time,y=cumHR))



#"weibull", "exponential", "gaussian", "logistic","lognormal" and "loglogistic"
coxregfit<-coxph(Surv(AVAL1, CNSR1) ~ 1,ties="breslow", data=data_in)
###模型检验（-2log(beita),AIC,BIC)  ####第一个数字是模型中参数数量 q的值
-2*coxregfit$loglik
extractAIC(coxregfit)
extractAIC(coxregfit,k=log(coxregfit$nevent))

fit1=coxregfit

setwd("~/help_for_others/DYQ/output/5_14m/")
dir.create("survival_regression")
pdf("./survival_regression/pfs_coxregfit.pdf",width = 5,height = 4)
plot(fit,conf.int=F,col = 4,lwd = 2,
     main="Kaplan-Meier estimate",
     xlab = "months", ylab = "survival function",
     cex.main=0.8)
text(x=5,y=0.2,label=paste("AIC=",round(extractAIC(coxregfit)[2],3),sep = ""))
dev.off()






weibullregfit <- survreg(Surv(AVAL1, CNSR1)~1, data = data_in,dist="weibull")
-2*weibullregfit $loglik
extractAIC(weibullregfit)
# extractAIC(weibullregfit,k=log(weibullregfit$))

fit1=weibullregfit
pct <- 1:99/100  
ptime <- predict(fit1,newdata = data.frame(1),type='quantile',p=pct, se=T)
temp0 <- data.frame(time=ptime$fit,
                    ymin=ptime$fit -1.96*ptime$se.fit,
                    ymax=ptime$fit +1.96*ptime$se.fit,
                    p=1-pct)

pdf("./survival_regression/pfs_weibullregfit.pdf",width = 5,height = 4)
plot(fit,conf.int=F,col = 4,lwd = 2,
     main="Kaplan-Meier estimate and weibull",
     xlab = "months", ylab = "survival function",
     cex.main=0.8)
lines(temp0$time, temp0$p, type = "l",lwd = 2,col=9)
lines(temp0$ymin, temp0$p, type = "l",lwd = 2,col=9,lty=3)
lines(temp0$ymax, temp0$p, type = "l",lwd = 2,col=9,lty=3)
legend("topright",pch = c(15,15),legend = c("KM-estimator","weibull"),
       col =c(4,9),bty="n")
text(x=5,y=0.2,label=paste("AIC=",round(extractAIC(weibullregfit)[2],3),sep = ""))
dev.off()





loglogiregfit <- survreg(Surv(AVAL1, CNSR1)~1, data = data_in,dist="loglogistic")
-2*loglogiregfit$loglik
extractAIC(loglogiregfit)

fit1=loglogiregfit
pct <- 1:99/100  ##定义死亡率
ptime <- predict(fit1,newdata = data.frame(1),type='quantile',p=pct, se=T)
temp0 <- data.frame(time=ptime$fit,
                    ymin=ptime$fit -1.96*ptime$se.fit,
                    ymax=ptime$fit +1.96*ptime$se.fit,
                    p=1-pct)

# dir.create("survival_regression")
pdf("./survival_regression/pfs_loglogiregfit.pdf",width = 5,height = 4)
plot(fit,conf.int=F,col = 4,lwd = 2,
     main="Kaplan-Meier estimate and loglogistic",
     xlab = "months", ylab = "survival function",
     cex.main=0.8)
lines(temp0$time, temp0$p, type = "l",lwd = 2,col=9)
lines(temp0$ymin, temp0$p, type = "l",lwd = 2,col=9,lty=3)
lines(temp0$ymax, temp0$p, type = "l",lwd = 2,col=9,lty=3)
legend("topright",pch = c(15,15),legend = c("KM-estimator","loglogistic"),
       col =c(4,9),bty="n")
text(x=5,y=0.2,label=paste("AIC=",round(extractAIC(loglogiregfit)[2],3),sep = ""))
dev.off()







lognormalregfit <- survreg(Surv(AVAL1, CNSR1)~1, data = data_in,dist="lognormal")
-2*lognormalregfit$loglik
extractAIC(lognormalregfit)

fit1=lognormalregfit
pct <- 1:99/100  ##定义死亡率
ptime <- predict(fit1,newdata = data.frame(1),type='quantile',p=pct, se=T)
temp0 <- data.frame(time=ptime$fit,
                    ymin=ptime$fit -1.96*ptime$se.fit,
                    ymax=ptime$fit +1.96*ptime$se.fit,
                    p=1-pct)

pdf("./survival_regression/pfs_lognormalregfit.pdf",width = 5,height = 4)
plot(fit,conf.int=F,col = 4,lwd = 2,
     main="Kaplan-Meier estimate and lognormal",
     xlab = "months", ylab = "survival function",
     cex.main=0.8)
lines(temp0$time, temp0$p, type = "l",lwd = 2,col=9)
lines(temp0$ymin, temp0$p, type = "l",lwd = 2,col=9,lty=3)
lines(temp0$ymax, temp0$p, type = "l",lwd = 2,col=9,lty=3)
legend("topright",pch = c(15,15),legend = c("KM-estimator","lognormal"),
       col =c(4,9),bty="n")
text(x=5,y=0.2,label=paste("AIC=",round(extractAIC(lognormalregfit)[2],3),sep = ""))
dev.off()






exponentialregfit <- survreg(Surv(AVAL1, CNSR1)~1, data = data_in,dist="exponential")
-2*exponentialregfit$loglik
extractAIC(exponentialregfit)


fit1=exponentialregfit
pct <- 1:99/100  ##定义死亡率
ptime <- predict(fit1,newdata = data.frame(1),type='quantile',p=pct, se=T)
temp0 <- data.frame(time=ptime$fit,
                    ymin=ptime$fit -1.96*ptime$se.fit,
                    ymax=ptime$fit +1.96*ptime$se.fit,
                    p=1-pct)

pdf("./survival_regression/pfs_exponentialregfit.pdf",width = 5,height = 4)
plot(fit,conf.int=F,col = 4,lwd = 2,
     main="Kaplan-Meier estimate and exponential",
     xlab = "months", ylab = "survival function",
     cex.main=0.8)
lines(temp0$time, temp0$p, type = "l",lwd = 2,col=9)
lines(temp0$ymin, temp0$p, type = "l",lwd = 2,col=9,lty=3)
lines(temp0$ymax, temp0$p, type = "l",lwd = 2,col=9,lty=3)
legend("topright",pch = c(15,15),legend = c("KM-estimator","exponential"),
       col =c(4,9),bty="n")
text(x=5,y=0.2,label=paste("AIC=",round(extractAIC(exponentialregfit)[2],3),sep = ""))
dev.off()









fit<-survfit(Surv(AVAL1, CNSR1)~1, data = data_in)
surv_probs <- summary(fit)$surv
time_points <- summary(fit)$time
rate_of_change <- diff(surv_probs) / diff(time_points)
dat_for_plot = data.frame(time=time_points[-1],rate_of_change=rate_of_change)
pdf("./survival_regression/rate_of_survival_change.pdf",width = 4,height = 3)
ggplot(dat_for_plot)+
  geom_line(aes(x=time,y=rate_of_change))+
  geom_point(aes(x=time,y=rate_of_change),shape=21)+
  geom_vline(xintercept=14,lty=2,col="red")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()
# plot(time_points[-1], rate_of_change, type = "b", xlab = "Time", ylab = "Rate of Change", 
#      main = "Rate of Change of Survival Probability")











fit1=loglogiregfit
summary(fit1)
intercept=fit1$coefficients
scale=fit1$scale

alpha=exp(-(intercept/scale))
gamma=1/scale
sur=function(t){
  return(1/(1+alpha*t^gamma))
}
sur(ptime$fit[1])
sur(ptime$fit[2])
sur(ptime$fit[50])
sur(1)
hazard=function(t){
  return((alpha*gamma*t^(gamma-1))/(1+alpha*t^gamma))
}
h_res = lapply(seq(0,50,0.1),function(t){
  return(hazard(t))
})
h_tb = data.frame(time=seq(0,50,0.1),hr=as.numeric(h_res))
pdf("./survival_regression/pfs_loglogistic_hazard.pdf",width = 4,height = 3)
ggplot(h_tb)+
  geom_line(aes(x=time,y=hr))+
  geom_vline(xintercept = c(8,14.2),lty=2,lwd=0.3,col="black")+
  geom_vline(xintercept = c(14),lty=2,col="red")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()

f=expression((alpha*gamma*t^(gamma-1))/(1+alpha*t^gamma))
# h_dao = deriv3(f,"t",function.arg = TRUE)
print(D(f,"t"))
h_dao = function(t){
  alpha * gamma * (t^((gamma - 1) - 1) * (gamma - 1))/(1 + alpha * t^gamma) - (alpha * gamma * t^(gamma - 1)) * (alpha * (t^(gamma - 1) * gamma))/(1 + alpha * t^gamma)^2
}
h_dao_res = lapply(seq(0,50,0.1),function(t){
  return(h_dao(t))
})
h_dao_tb = data.frame(time=seq(0,50,0.1),hr_d=as.numeric(h_dao_res))
pdf("./survival_regression/pfs_loglogistic_hazard_D.pdf",width = 4,height = 3)
ggplot(h_dao_tb)+
  geom_line(aes(x=time,y=hr_d))+
  geom_vline(xintercept = c(8,14.2),lty=2,lwd=0.3,col="black")+
  geom_vline(xintercept = c(14),lty=2,col="red")+
  geom_hline(yintercept = c(0),lty=2,col="black")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()


