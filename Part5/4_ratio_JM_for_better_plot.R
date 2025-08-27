rm(list = ls())
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



#########################
load("../1_preclean/data/data_C1_preclean.rdata")
mt_info = data_in_C1 %>% dplyr::select(Compounds,`Class I`,`Class II`) %>% 
  dplyr::rename(mt=Compounds)
# VIP_FC = sel %>% inner_join(mt_info,by="mt")

class_ID_1 = mt_info %>% dplyr::select(mt,`Class I`) %>% dplyr::rename(class_1=colnames(.)[2])
class_ID_2 = mt_info %>% dplyr::select(mt,`Class II`) %>% dplyr::rename(class=colnames(.)[2])

final = final %>% dplyr::rename(mt=metabo) %>% inner_join(class_ID_2,by = "mt") %>% inner_join(class_ID_1,by = "mt")

load("./ratio_pred/all_mt_use_for_model.rdata")

SP = c(Cer,SM)




####################################################
####################################################
label_data1 = final %>% dplyr::filter(mt %in% PL) %>% dplyr::mutate(kind="PL")
label_data2 = final %>% dplyr::filter(mt %in% SP) %>% dplyr::mutate(kind="SP")

label_data = rbind(label_data1,label_data2)

label_data = label_data %>% arrange(kind,coef_R1) %>% dplyr::mutate(X=seq(1,nrow(.),1))
label_data$mt = factor(label_data$mt,levels = unique(label_data$mt),ordered = T)

# label_data <- 
number_of_bar <- nrow(label_data)

angle <-  90 - 360 * (label_data$X-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)
# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
#查看标签数据
head(label_data)


library(ggnewscale)
p <- ggplot(filter(label_data,kind=="PL"), aes(x=as.factor(mt), y=abs(coef_R1),fill=abs(coef_R1))) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  # This add the bars with a blue color
  geom_bar(stat="identity") +
  scale_fill_gradient2(low = "white",mid = "#E6E6FA",high = "#483D8B")+
  new_scale_fill()+
  # Note that id is a factor. If x is numeric, there is some space between the first bar
  # This add the bars with a blue color
  geom_bar(inherit.aes = FALSE,data=filter(label_data,kind=="SP"), aes(x=as.factor(mt), y=abs(coef_R1),fill=abs(coef_R1)),stat="identity") +
  scale_fill_gradient2(low = "white",mid = "#DDA0DD",high = "#8B008B")+
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-0.3,2) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
    # plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=as.factor(mt), y=abs(coef_R1)+0.01,label=mt, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=2, 
            angle= label_data$angle, inherit.aes = FALSE ) #+
# geom_text(data=label_data, aes(x=as.factor(mt), y=coef_R1-0.18,label=round(coef_R1,2), hjust=hjust), 
#           color="black", fontface="bold",alpha=0.6, size=2.5, 
#           angle= label_data$angle, inherit.aes = FALSE ) 

pdf("./ratio_pred/cor_barplot_all_coef_R1.pdf",width = 14,height = 10)
p
dev.off()



####################################################
####################################################
label_data1 = final %>% dplyr::filter(mt %in% PL) %>% dplyr::mutate(kind="PL")
label_data2 = final %>% dplyr::filter(mt %in% SP) %>% dplyr::mutate(kind="SP")

label_data = rbind(label_data1,label_data2)

label_data = label_data %>% arrange(kind,coef_R2) %>% dplyr::mutate(X=seq(1,nrow(.),1))
label_data$mt = factor(label_data$mt,levels = unique(label_data$mt),ordered = T)

# label_data <- 
number_of_bar <- nrow(label_data)

angle <-  90 - 360 * (label_data$X-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)
# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
#查看标签数据
head(label_data)


# 绘制带标签的环状条形图
p <- ggplot(filter(label_data,kind=="PL"), aes(x=as.factor(mt), y=log(abs(coef_R2)+1),fill=log(abs(coef_R2)+1))) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  # This add the bars with a blue color
  geom_bar(stat="identity") +
  scale_fill_gradient2(low = "white",mid = "#E6E6FA",high = "#483D8B")+
  new_scale_fill()+
  # Note that id is a factor. If x is numeric, there is some space between the first bar
  # This add the bars with a blue color
  geom_bar(inherit.aes = FALSE,data=filter(label_data,kind=="SP"), aes(x=as.factor(mt), y=log(abs(coef_R2)+1),fill=log(abs(coef_R2)+1)),stat="identity") +
  scale_fill_gradient2(low = "white",mid = "#DDA0DD",high = "#8B008B")+
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-0.62,4.2) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
    # plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=as.factor(mt), y=log(abs(coef_R2)+1)+0.01,label=mt, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=2, 
            angle= label_data$angle, inherit.aes = FALSE ) #+
# geom_text(data=label_data, aes(x=as.factor(mt), y=coef_R1-0.18,label=round(coef_R1,2), hjust=hjust), 
#           color="black", fontface="bold",alpha=0.6, size=2.5, 
#           angle= label_data$angle, inherit.aes = FALSE ) 

pdf("./ratio_pred/cor_barplot_all_coef_R2.pdf",width = 16,height = 12)
p
dev.off()

