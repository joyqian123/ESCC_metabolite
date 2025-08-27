rm(list = ls())
setwd("~/help_for_others/DYQ/output/8_compare_with_mt_change_in_R_and_NR/")
library(tidyverse)
dir.create("metabo_process")
setwd("./metabo_process/")
dir.create("./all_C1C6")
setwd("./all_C1C6/")

load("~/help_for_others/DYQ/output/1_preclean/data/data_preclean.rdata")
table(clinic_in$cycle)



clinic_1 = clinic_in %>% dplyr::filter(cycle %in% c("C1D1","C6D1"))


# .libPaths()
library(MetaboAnalystR)
library(tidyverse)

mSet<-InitDataObjects("conc", "stat", FALSE)
mSet<-Read.TextData(mSet, "../../csv/data_all_C1C6.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet)  
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE,ratioNum = 20)

mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   

norm = mSet$dataSet$norm
info = mSet$dataSet$meta.info %>% unlist()




norm_1 = norm %>% rownames_to_column("Sample") %>% 
  inner_join(clinic_1[,c("Sample","SUBJID","cycle")],by = "Sample") %>% 
  dplyr::select(Sample,SUBJID,cycle,everything())
# norm_2 = norm_1 %>% pivot_longer(cols = !c(1:3),names_to = "mt",values_to="value")

clinic_1 = clinic_1 %>% dplyr::select(SUBJID,R) %>% distinct(.,.keep_all = T)
norm_2 = norm_1%>% inner_join(clinic_1[,c("SUBJID","R")],by = "SUBJID") %>% dplyr::select(Sample,SUBJID,R,cycle,everything())


hvg = apply(norm_2[,-c(1:4)],2,var) %>% data.frame() %>% dplyr::rename(v=colnames(.)[1]) %>% rownames_to_column("mt")
hvg_mt = hvg %>% arrange(v)%>% pull(mt)

norm_use = norm_2 %>% dplyr::select(any_of(hvg_mt))
library(factoextra)
library(FactoMineR)
# df.pca<- PCA(norm_use, graph = FALSE)
group <- str_c(norm_2$cycle,norm_2$R,sep = "_")

pca_res <- prcomp(norm_use, scale. = TRUE)  
explained_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

pc1_var <- explained_var[1]  
pc2_var <- explained_var[2]  

pca_data <- as.data.frame(pca_res$x)
pca_data$Group <- group  
centroids <- pca_data %>%
  group_by(Group) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))
sample_points = separate(centroids,col = Group,into = c("cycle","response"),sep = "_") %>% 
  pivot_wider(id_cols = response,names_from = cycle,values_from = c("PC1","PC2"))



library(wesanderson)
pdf("../../plot/pca_cycle_change.pdf",width = 5,height = 3.5)
ggplot() +
  # geom_point(alpha = 0.6) +
  geom_point(data = centroids, aes(x = PC1, y = PC2,col=Group), 
             shape = 20, size = 5,stroke = 1.5) +
  geom_text(data = centroids, aes(x = PC1, y = PC2, label = Group), 
            vjust = -1, color = "black") +
  scale_color_manual(values = rep(wes_palette("GrandBudapest1",2),2))+
  geom_vline(xintercept = 0,lty=3,lwd=0.5,col="#4B0082")+
  geom_hline(yintercept = 0,lty=3,lwd=0.5,col="#4B0082")+
  geom_segment(data = sample_points, 
               aes(x = PC1_C1D1 + 0.05, xend = PC1_C6D1 - 0.05, 
                   y = PC2_C1D1 , yend = PC2_C6D1 ),
               arrow = arrow(length = unit(0.1, "inches")), lwd=0.5,lty=2,
               color = wes_palette("GrandBudapest1",3)[3])+
  xlim(-6,6)+
  ylim(-1,0.75)+
  theme_bw()
dev.off()
