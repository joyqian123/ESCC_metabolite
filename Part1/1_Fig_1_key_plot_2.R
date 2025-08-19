rm(list = ls())
setwd("~/help_for_others/DYQ/")  #######在最开始设置工作目录
library(MetaboAnalystR)

library(tidyverse)

setwd("./output/2_response_mt_enricher/data_enr/")
# Create mSetObj
mSet<-InitDataObjects("conc", "msetora", FALSE)
mSet<-Read.TextData(mSet, "../../1_preclean/data/C1_df.csv", "rowu", "disc");

mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet)  ####替换0值或缺失值
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "S10T0", ratio=FALSE, ratioNum=20)

mSet<-PlotNormSummary(mSet, "norm_0_", format ="pdf", dpi=72, width=NA)   ###按代谢物作图
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "pdf", dpi=72, width=NA)   ###按样本作图

norm = mSet$dataSet$norm
group = mSet$dataSet$cls


#########################广义线性回归
kk=1
res = lapply(1:ncol(norm),function(kk){
  norm_use = cbind(data.frame(response=group),norm[,kk]) %>% dplyr::rename(mt=colnames(.)[2])
  model <- glm(response ~ ., data = norm_use, family = binomial)
  s=summary(model)
  return(data.frame(mt=colnames(norm)[kk],z=s$coefficients["mt",3],p=s$coefficients["mt",4]))
}) %>% do.call(rbind,.)

res$p.adj = p.adjust(res$p,method = "BH")


res = res %>% arrange(z) %>% dplyr::mutate(rank=seq(1,nrow(.),1))


save(res,file = "../glm/res.rdata")



load("../glm/res.rdata")
load("../../1_preclean/data/data_C1_preclean.rdata")
mt_info = data_in_C1 %>% dplyr::select(Compounds,`Class I`,`Class II`) %>% 
  dplyr::rename(mt=Compounds)
res = res %>% inner_join(mt_info,by="mt")


class_ID_1 = res %>% dplyr::select(mt,`Class I`) %>% dplyr::rename(class=colnames(.)[2])
class_ID_2 = res %>% dplyr::select(mt,`Class II`) %>% dplyr::rename(class=colnames(.)[2])
# class_ID_2$class = str_remove(class_ID_2$class,"-O")
# class_ID_2$class = str_remove(class_ID_2$class,"-P")
class_ID_2$class[grepl("Cer",class_ID_2$class)] = "Cer"
# class_ID_2$class[grepl("LP",class_ID_2$class)] = "LP"
# class_ID_2$class[grepl("LPE",class_ID_2$class)] = "LPC/LPE"

table(class_ID_2$class)
# class_ID_2_hex = class_ID_2 %>% dplyr::filter()


class_ID = rbind(class_ID_1,class_ID_2)
unique(class_ID_1$class)
unique(class_ID_2$class)

geneList= res$z #%>% rank()
# geneList = (geneList-min(geneList))/(max(geneList)-min(geneList))-0.5
names(geneList) = res$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=class_ID_1[,c(2,1)], 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_class_I = egmt@result

egmt <- GSEA(geneList, TERM2GENE=class_ID_2[,c(2,1)], 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_class_II = egmt@result

egmt_res_class_II = separate(egmt_res_class_II,col = "leading_edge",sep = ", ",into = c("tag","list","signal"))
egmt_res_class_II$signal = str_extract(egmt_res_class_II$signal,"[0-9]+") %>% as.numeric()
egmt_res_class_II$core_enrichment_N = sapply(1:nrow(egmt_res_class_II),function(kk){
  tmp = egmt_res_class_II[kk,] %>% pull(core_enrichment) %>% strsplit("/") %>% unlist() %>% length()
  return(tmp)
})
egmt_res_class_II = egmt_res_class_II %>% dplyr::filter(signal>30,core_enrichment_N>1)
  
save(egmt_res_class_I,egmt_res_class_II,file = "../glm/egmt_res.rdata")

# load("egmt_res.rdata")
# options(scipen = 200)


myPalette <- colorRampPalette(c("blue","green","yellow","red"))
gradientColors <- myPalette(1000)

library(ggrepel)

# egmt_res_class_II_R = egmt_res_class_II %>% dplyr::filter(NES>0)

pdf("../plot/egmt.pdf",width = 5,height = 4)
ggplot(egmt_res_class_II,aes(x=NES,y=-log(pvalue)))+
  geom_point(aes(size=-log(pvalue),fill=abs(NES)),shape=21)+
  scale_fill_gradientn(colours = gradientColors)+
  geom_text_repel(aes(label=ifelse(pvalue<0.05,ID,"")),min.segment.length = 0.1)+
  geom_hline(yintercept = -log(0.05),lty=2)+
  geom_vline(xintercept = 0,lty=2)+
  ggtitle('metabolites gsea')+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"), # 添加 x 和 y 轴线 
        axis.ticks = element_line(color = "black"),
        # legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()



# "#FF69B4","#6A5ACD"
pdf("../plot/glm_res.pdf",width = 4,height = 4)
ggplot(res,aes(x=rank,y=z))+
  geom_point(aes(col=ifelse(p<0.05,ifelse(z>0,"#e5965a","#746fb1"),"lightgrey")),size=1)+
  scale_color_identity()+
  geom_text_repel(aes(label=ifelse(p<0.05,mt,"")),size=3,min.segment.length = 0.1)+
  geom_hline(yintercept = 0,lwd=0.5,lty=2)+
  theme_bw() + # 使用主题基础 
  theme( panel.background = element_blank(), # 移除面板背景 
         panel.grid = element_blank(), # 移除网格线 
         axis.title = element_text(size = 12), # 调整轴标题大小 
         axis.line = element_line(color = "black"), # 添加 x 和 y 轴线 
         axis.ticks = element_line(color = "black"), # 保留刻度线 
         legend.position = "none") # 移除图例 )
dev.off()






# load("../../../reference/kegg_ID_sum.rdata")
# up_R = res %>% dplyr::filter(p<0.1,z>0) %>% pull(mt)
# res_R = enricher(up_R,TERM2GENE = class_ID_2[,c(2,1)])
# res_R = res_R@result
# 
# plot = separate(data = res_R,col = "GeneRatio",sep = "/",into = c("A","B"),remove = FALSE) %>% 
#   dplyr::mutate(A=as.numeric(A),B=as.numeric(B)) %>% dplyr::mutate(r=A/B)
# myPalette <- colorRampPalette(c("blue","green","yellow","red"))
# gradientColors <- myPalette(1000)
# 
# pdf("../plot/enricher_R.pdf",width = 5,height = 4)
# ggplot(plot,aes(x=r,y=-log(pvalue)))+
#   geom_point(aes(size=-log(pvalue),fill=-log(pvalue)),shape=21)+
#   scale_fill_gradientn(colours = gradientColors)+
#   geom_text_repel(aes(label=ifelse(pvalue<0.1,ID,"")),min.segment.length = 0.1)+
#   geom_hline(yintercept = -log(0.1),lty=2)+
#   geom_vline(xintercept = 0,lty=2)+
#   ggtitle('metabolites gsea')+
#   theme_bw()+
#   theme(panel.grid.major=element_line(colour=NA),
#         panel.border = element_blank(),
#         axis.line = element_line(color = "black"), # 添加 x 和 y 轴线 
#         axis.ticks = element_line(color = "black"),
#         # legend.position = "none",
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank())
# dev.off()
# 
# 
# 
# up_NR = res %>% dplyr::filter(p<0.1,z<0) %>% pull(mt)
# res_NR = enricher(up_NR,TERM2GENE = class_ID_2[,c(2,1)])
# res_NR = res_NR@result
# 
# plot = separate(data = res_NR,col = "GeneRatio",sep = "/",into = c("A","B"),remove = FALSE) %>% 
#   dplyr::mutate(A=as.numeric(A),B=as.numeric(B)) %>% dplyr::mutate(r=A/B)
# myPalette <- colorRampPalette(c("blue","green","yellow","red"))
# gradientColors <- myPalette(1000)
# 
# pdf("../plot/enricher_R.pdf",width = 5,height = 4)
# ggplot(plot,aes(x=r,y=-log(pvalue)))+
#   geom_point(aes(size=-log(pvalue),fill=-log(pvalue)),shape=21)+
#   scale_fill_gradientn(colours = gradientColors)+
#   geom_text_repel(aes(label=ifelse(pvalue<0.1,ID,"")),min.segment.length = 0.1)+
#   geom_hline(yintercept = -log(0.1),lty=2)+
#   geom_vline(xintercept = 0,lty=2)+
#   ggtitle('metabolites gsea')+
#   theme_bw()+
#   theme(panel.grid.major=element_line(colour=NA),
#         panel.border = element_blank(),
#         axis.line = element_line(color = "black"), # 添加 x 和 y 轴线 
#         axis.ticks = element_line(color = "black"),
#         # legend.position = "none",
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank())
# dev.off()




# class = egmt_res_class_II %>% dplyr::filter(pvalue<0.05) %>% pull(ID)
# label_id = egmt_res_class_II %>% dplyr::filter(pvalue<0.05) %>% pull(core_enrichment) 

library(ComplexHeatmap)
library(circlize)
kk=3



mt_ls_1 = up_R = res %>% dplyr::filter(p<0.05,z>0) %>% pull(mt)
mt_ls_2 = up_NR = res %>% dplyr::filter(p<0.05,z<0) %>% pull(mt)

norm_sel = norm %>% dplyr::select(any_of(mt_ls_1),any_of(mt_ls_2)) #%>% cbind(response=group,.) 
scale_tb = lapply(1:ncol(norm_sel),function(pp){
  a = norm_sel[,pp]
  tmp = data.frame(mt=(a-min(a))/(max(a)-min(a)))
}) %>% do.call(cbind,.)
mt_ls = c(mt_ls_1,mt_ls_2)
colnames(scale_tb)=mt_ls
scale_tb = scale_tb %>% t() %>% data.frame(.,check.names = FALSE) %>% 
  distinct(.,.keep_all = T)

mt_ls_1 = mt_ls_1 %>% intersect(rownames(scale_tb))
mt_ls_2 = mt_ls_2 %>% intersect(rownames(scale_tb))

mt_ls_class_1 = class_ID_1 %>% dplyr::filter(mt %in% mt_ls_1) %>% arrange(class)
mt_ls_class_2 = class_ID_1 %>% dplyr::filter(mt %in% mt_ls_2) %>% arrange(class)

scale_tb = scale_tb[c(mt_ls_class_1$mt,mt_ls_class_2$mt),]
scale_tb_R = scale_tb[c(mt_ls_class_1$mt),]
scale_tb_NR = scale_tb[c(mt_ls_class_2$mt),]
scale_tb = rbind(scale_tb_NR,scale_tb_R) %>% rownames_to_column("mt") %>% 
  dplyr::mutate(mt=factor(mt,levels = mt,ordered = T)) %>% column_to_rownames("mt")

anno_mt = c("(R)-(-)-1-Amino-2-propanol","3-Methoxyphenylacetic acid",
            "TXB2","Tetraethylene-glycol","Hexaethylene-glycol","Heptethylene-glycol",
            "TG(18:2_20:5_22:6)","TG(18:2_22:5_22:6)",
            "S-Allyl-L-cysteine","N-Cinnamylglycine","4-Phenylbutyric acid",
            "1,5-Anhydro-D-Glucitol","PI(18:0_22:4)","PI(16:0_16:1)","PI(15:0_17:1)",
            "PA(18:0_20:3)","PA(18:0_22:4)","PC(18:0_20:2)","PC(16:0_20:2)","PC(18:0_22:1)",
            "TG(15:0_18:0_18:1)","TG(16:0_18:1_18:1)","TG(14:0_18:0_18:2)",
            "Lithocholic acid","O-Phospho-L-tyrosine")

library(colorspace)
swatchplot(qualitative_hcl(12, palette = "Harmonic"))
# col_use = qualitative_hcl(12, palette = "Set 3")
col_use = c("#478bbe","#e6966a","#7c6aa4","#83c5db","#aedacb","#aeb0d8","#d75821","#f7bfa4",
            "#e6b119","#b5c574","#b65741","#4eb69e")
names(col_use) = unique(c(mt_ls_class_2$class,mt_ls_class_1$class))

col_ha = HeatmapAnnotation(response=group,col = list(response=c("NR"="#b5a1e3","R"="#f0c2a2")))
row_ha = rowAnnotation(mt=c(rep("R",length(mt_ls_1)),rep("NR",length(mt_ls_2))),
                       type = c(mt_ls_class_2$class,mt_ls_class_1$class),
                       col = list(mt=c("NR"="#b5a1e3","R"="#f0c2a2"),
                                  type = col_use))

# col_fun = colorRamp2(c(0,0.01,1), c("navy", "white", "red"))
# Heatmap(t(scale_tb),top_annotation = col_ha,cluster_columns = FALSE,
#         column_names_gp = gpar(fontsize = 0))  
color = colorRampPalette(c("#9E423A","#E17B64","#DFD0A9","#B3D2C7","#72BFE6","#0065AA"))(100) %>% rev()
# ,)(100)
color = colorRampPalette(c("navy", "white", "red"))(100) 

pdf("../plot/heatmap.pdf",width = 10,height = 8)
Heatmap(scale_tb,top_annotation = col_ha,right_annotation = row_ha,
        cluster_columns = TRUE,
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 0),
        column_split = group,
        row_split = c(rep("NR",length(mt_ls_2)),rep("R",length(mt_ls_1))),
        col = color,
        # row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 0))   + 
  rowAnnotation(link = anno_mark(at = which(rownames(scale_tb) %in% anno_mt), 
                                 labels = rownames(scale_tb)[which(rownames(scale_tb) %in% anno_mt)], 
                                 labels_gp = gpar(fontsize = 10)))
dev.off()



############################################################################################
############################################################################################
label_id = egmt_res_class_II %>% dplyr::filter(pvalue<0.05) %>% pull(core_enrichment)


# LPC = label_id[[1]] %>% strsplit(.,"/LPC") %>% unlist()
# LPC[-1] = str_c("LPC",LPC[-1],sep = "")
# LPC = rev(LPC)
# library(ggpubr)
# library(ggpmisc)
# norm_sel = norm %>% dplyr::select(any_of(LPC[c(1:7)])) %>% rownames_to_column("Sample")
# response = data.frame(Sample=rownames(norm),response=group)
# library(RobustRankAggreg)
# method_choose="RRA"
# glist<-list()
# i=2
# for (i in 2:ncol(norm_sel)){
#   data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
#   glist[[i-1]]<-data$Sample
#   names(glist)[i-1]<-colnames(norm_sel)[i]
# }
# library(RobustRankAggreg)
# r<-RobustRankAggreg::rankMatrix(glist,full = TRUE)
# # RRA<-aggregateRanks(glist,rmat=r,
# #                         method=method_choose,full=TRUE,exact=F,topCutoff=NA)
# 
# RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
# # RRA$Name = rownames(r)
# # RRA = data.frame(Name=norm_sel$Sample,Score=apply(norm_sel[,-1],1,mean))
# 
# RRA = RRA %>% dplyr::rename(Sample=Name) %>% inner_join(response,by = "Sample")
# RRA_LPC = RRA %>% dplyr::mutate(mt="LPC")
# t.test(RRA$Score~RRA$response)
# 
# ggplot(RRA,aes(x=response,y=Score,col=response))+
#   geom_boxplot()+
#   geom_jitter()+
#   # geom_segment(aes(x=response-0.5,xend=response+0.5,y=rank,yend = rank,col=response))+
#   scale_color_manual(values = c("#FF69B4","#6A5ACD"))+
#   stat_compare_means(comparisons = list(c("R","NR")))+
#   # facet_grid(~mt)+
#   theme_bw()#+
#   # theme(strip.background = element_blank())



oxlipid = label_id[[2]] %>% strsplit(.,"/") %>% unlist()
oxlipid = oxlipid %>% rev() #%>% .[1:4]
# LPC[-1] = str_c("LPC",LPC[-1],sep = "")

library(ggpubr)
library(ggpmisc)
norm_sel = norm %>% dplyr::select(any_of(oxlipid)) %>% rownames_to_column("Sample")
response = data.frame(Sample=rownames(norm),response=group)
library(RobustRankAggreg)
method_choose="mean"
glist<-list()
i=2
for (i in 2:ncol(norm_sel)){
  data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
  glist[[i-1]]<-data$Sample
  names(glist)[i-1]<-colnames(norm_sel)[i]
}
library(RobustRankAggreg)
r<-RobustRankAggreg::rankMatrix(glist,full = TRUE)
# RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
RRA<-aggregateRanks(glist,rmat=r,
                        method=method_choose,full=TRUE,exact=F,topCutoff=NA)

RRA = RRA %>% dplyr::rename(Sample=Name) %>% inner_join(response,by = "Sample")
RRA_ox = RRA %>% dplyr::mutate(mt="oxlipid")
t.test(RRA_ox$Score~RRA_ox$response)
# p1=ggplot(RRA,aes(x=response,y=Score,col=response))+
#   geom_boxplot()+
#   geom_jitter()+
#   ggtitle("oxlipid")+
#   scale_color_manual(values = c("#FF69B4","#6A5ACD"))+
#   stat_compare_means(comparisons = list(c("R","NR")))+
#   theme_bw()+
#   theme(panel.grid.major=element_line(colour=NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         legend.position = "none")




PC = label_id[[5]] %>% strsplit(.,"/") %>% unlist()
# LPC[-1] = str_c("LPC",LPC[-1],sep = "")

library(ggpubr)
library(ggpmisc)
norm_sel = norm %>% dplyr::select(any_of(PC)) %>% rownames_to_column("Sample")
response = data.frame(Sample=rownames(norm),response=group)
library(RobustRankAggreg)
method_choose="mean"
glist<-list()
i=2
for (i in 2:ncol(norm_sel)){
  data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
  glist[[i-1]]<-data$Sample
  names(glist)[i-1]<-colnames(norm_sel)[i]
}
library(RobustRankAggreg)
r<-rankMatrix(glist,full = TRUE)
RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
# RRA<-aggregateRanks(glist,rmat=r,
#                     method=method_choose,full=TRUE,exact=F,topCutoff=NA)

RRA = RRA %>% dplyr::rename(Sample=Name) %>% inner_join(response,by = "Sample")

RRA_PC = RRA %>% dplyr::mutate(mt="PC")
t.test(RRA_PC$Score~RRA_PC$response)
# p2=ggplot(RRA,aes(x=response,y=Score,col=response))+
#   geom_boxplot()+
#   geom_jitter()+
#   ggtitle("PC")+
#   scale_color_manual(values = c("#FF69B4","#6A5ACD"))+
#   stat_compare_means(comparisons = list(c("R","NR")))+
#   theme_bw()+
#   theme(panel.grid.major=element_line(colour=NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         legend.position = "none")





PI = label_id[[4]] %>% strsplit(.,"/") %>% unlist()

library(ggpubr)
library(ggpmisc)
norm_sel = norm %>% dplyr::select(any_of(PI)) %>% rownames_to_column("Sample")
response = data.frame(Sample=rownames(norm),response=group)
library(RobustRankAggreg)
method_choose="mean"
glist<-list()
i=2
for (i in 2:ncol(norm_sel)){
  data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
  glist[[i-1]]<-data$Sample
  names(glist)[i-1]<-colnames(norm_sel)[i]
}
library(RobustRankAggreg)
r<-rankMatrix(glist,full = TRUE)
RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
# RRA<-aggregateRanks(glist,rmat=r,
#                     method=method_choose,full=TRUE,exact=F,topCutoff=NA)

RRA = RRA %>% dplyr::rename(Sample=Name) %>% inner_join(response,by = "Sample")
RRA_PI = RRA %>% dplyr::mutate(mt="PI")
t.test(RRA_PI$Score~RRA_PI$response)

# p3=ggplot(RRA,aes(x=response,y=Score,col=response))+
#   geom_boxplot()+
#   geom_jitter()+
#   ggtitle("PI")+
#   scale_color_manual(values = c("#FF69B4","#6A5ACD"))+
#   stat_compare_means(comparisons = list(c("R","NR")))+
#   theme_bw()+
#   theme(panel.grid.major=element_line(colour=NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         legend.position = "none")






sugar_alchol = label_id[[1]] %>% strsplit(.,"/") %>% unlist()

library(ggpubr)
library(ggpmisc)
norm_sel = norm %>% dplyr::select(any_of(sugar_alchol)) %>% rownames_to_column("Sample")
response = data.frame(Sample=rownames(norm),response=group)
library(RobustRankAggreg)
method_choose="mean"
glist<-list()
i=2
for (i in 2:ncol(norm_sel)){
  data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
  glist[[i-1]]<-data$Sample
  names(glist)[i-1]<-colnames(norm_sel)[i]
}
library(RobustRankAggreg)
r<-rankMatrix(glist,full = TRUE)
RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])

RRA = RRA %>% dplyr::rename(Sample=Name) %>% inner_join(response,by = "Sample")

RRA_sugaralcohol = RRA %>% dplyr::mutate(mt="sugar_alcohol")
t.test(RRA_sugaralcohol$Score~RRA_sugaralcohol$response)
# p4=ggplot(RRA,aes(x=response,y=Score,col=response))+
#   geom_boxplot()+
#   geom_jitter()+
#   ggtitle("sugar_alcohol")+
#   scale_color_manual(values = c("#FF69B4","#6A5ACD"))+
#   stat_compare_means(comparisons = list(c("R","NR")))+
#   theme_bw()+
#   theme(panel.grid.major=element_line(colour=NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         legend.position = "none")






alchol = label_id[[3]] %>% strsplit(.,"/") %>% unlist()

library(ggpubr)
library(ggpmisc)
norm_sel = norm %>% dplyr::select(any_of(alchol)) %>% rownames_to_column("Sample")
response = data.frame(Sample=rownames(norm),response=group)
library(RobustRankAggreg)
method_choose="mean"
glist<-list()
i=2
for (i in 2:ncol(norm_sel)){
  data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
  glist[[i-1]]<-data$Sample
  names(glist)[i-1]<-colnames(norm_sel)[i]
}
library(RobustRankAggreg)
r<-RobustRankAggreg::rankMatrix(glist,full = TRUE)
RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
RRA<-aggregateRanks(glist,rmat=r,
                    method=method_choose,full=TRUE,exact=F,topCutoff=NA)

RRA = RRA %>% dplyr::rename(Sample=Name) %>% inner_join(response,by = "Sample")

RRA_alcohol = RRA %>% dplyr::mutate(mt="alcohol")
t.test(RRA_alcohol$Score~RRA_alcohol$response)


# p5=ggplot(RRA,aes(x=response,y=Score,col=response))+
#   geom_boxplot()+
#   geom_jitter()+
#   ggtitle("alchol")+
#   scale_color_manual(values = c("#FF69B4","#6A5ACD"))+
#   stat_compare_means(comparisons = list(c("R","NR")))+
#   theme_bw()+
#   theme(panel.grid.major=element_line(colour=NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         legend.position = "none")





# hex = label_id[[5]] %>% strsplit(.,"/HexCer") %>% unlist()
# hex[-1] = str_c("HexCer",hex[-1],sep = "")
# 
# library(ggpubr)
# library(ggpmisc)
# norm_sel = norm %>% dplyr::select(any_of(hex)) %>% rownames_to_column("Sample")
# response = data.frame(Sample=rownames(norm),response=group)
# library(RobustRankAggreg)
# method_choose="mean"
# glist<-list()
# i=2
# for (i in 2:ncol(norm_sel)){
#   data<-norm_sel %>% dplyr::select(1,any_of(i)) %>% dplyr::rename(mt=colnames(.)[2]) %>% arrange(mt)
#   glist[[i-1]]<-data$Sample
#   names(glist)[i-1]<-colnames(norm_sel)[i]
# }
# library(RobustRankAggreg)
# r<-RobustRankAggreg::rankMatrix(glist,full = TRUE)
# RRA = rowMeans(r) %>% data.frame() %>% rownames_to_column("Name") %>% dplyr::rename(Score=colnames(.)[2])
# 
# RRA = RRA %>% dplyr::rename(Sample=Name) %>% inner_join(response,by = "Sample")

# ggplot(RRA,aes(x=response,y=Score,col=response))+
#   geom_boxplot()+
#   geom_jitter()+
#   ggtitle("alchol")+
#   scale_color_manual(values = c("#FF69B4","#6A5ACD"))+
#   stat_compare_means(comparisons = list(c("R","NR")))+
#   theme_bw()+
#   theme(legend.position = "none")


RRA = rbind(RRA_alcohol,RRA_ox,RRA_PC,RRA_PI,RRA_sugaralcohol)


library(ggpubr)
library(ggpmisc)
pdf("../plot/5_class_RRA_predict_new.pdf",width = 12,height = 4)
library(wesanderson)
# median_tb <- RRA %>%
#   group_by(mt,response) %>%
#   dplyr::summarise(median_expr = median(Score, na.rm = TRUE)) %>%
#   ungroup()

ggplot(RRA,aes(x=response,y=Score))+
  # geom_violin(aes(fill=response),alpha=1,scale = "width",trim = TRUE)+
  geom_boxplot(aes(color=response),alpha=1,outliers=FALSE)+
  # geom_jitter()+
  stat_compare_means(comparisons = list(c("NR","R")),method = "t.test")+
  # geom_segment(data = median_tb,  # 使用提前计算的中位数数据
  #              aes(x= as.numeric(as.factor(response)) - 0.15,  # 调整x的值来控制线段的起点
  #                  xend = as.numeric(as.factor(response)) + 0.15,  # 调整xend的值来控制线段的终点
  #                  y = median_expr,
  #                  yend = median_expr),
  #              linewidth=1,lty=1) +
  scale_fill_manual(values = c("#b5a1e3","#f0c2a2"))+
  scale_color_manual(values = c("#b5a1e3","#f0c2a2"))+
  xlab("abundance")+
  facet_wrap(~mt,ncol=5,scale="free_y")+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        strip.background = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()

# 
# pdf("../plot/5_class_RRA_predict.pdf",width = 12,height = 4)
# cowplot::plot_grid(p1,p5,p2,p3,p4,nrow = 1)
# dev.off()
# 


# 
# RRA_LPC = RRA_LPC %>% dplyr::rename(LPC=colnames(.)[2])
# RRA_PC = RRA_PC %>% dplyr::rename(PC=colnames(.)[2]) %>% inner_join(RRA_LPC,by = c("Sample","response")) %>% 
#   dplyr::mutate(rel=PC/LPC)
# 
# pdf("../plot/PC_LPC_ratio_class_RRA_predict.pdf",width = 3,height = 4)
# ggplot(RRA_PC,aes(x=response,y=log(rel),col=response))+
#   geom_boxplot()+
#   geom_jitter()+
#   ggtitle("PC/LPC")+
#   scale_color_manual(values = c("#FF69B4","#6A5ACD"))+
#   stat_compare_means(comparisons = list(c("R","NR")))+
#   theme_bw()+
#   theme(legend.position = "none")
# dev.off()





######################################
r = pROC::roc(response$response,norm$`S-Allyl-L-cysteine`,levels=c("NR","R"),direction="<")
auc = r %>% pROC::auc()
auc_ci <- pROC::ci(r)
tb1 = data.frame(mt="S-Allyl-L-cysteine",auc=auc,auc_h=auc_ci[3],auc_l=auc_ci[1],col="#e5965a")

r = pROC::roc(response$response,norm$`1,5-Anhydro-D-Glucitol`)
auc = r %>% pROC::auc()
auc_ci <- pROC::ci(r)
tb2 = data.frame(mt="1,5-Anhydro-D-Glucitol",auc=auc,auc_h=auc_ci[3],auc_l=auc_ci[1],col="#e5965a")

r = pROC::roc(response$response,norm$`N-Cinnamylglycine`)
auc = r %>% pROC::auc()
auc_ci <- pROC::ci(r)
tb3 = data.frame(mt="N-Cinnamylglycine",auc=auc,auc_h=auc_ci[3],auc_l=auc_ci[1],col="#e5965a")

r = pROC::roc(response$response,norm$`(R)-(-)-1-Amino-2-propanol`)
auc = r %>% pROC::auc()
auc_ci <- pROC::ci(r)
tb4 = data.frame(mt="(R)-(-)-1-Amino-2-propanol",auc=auc,auc_h=auc_ci[3],auc_l=auc_ci[1],col="#746fb1")

r = pROC::roc(response$response,norm$`Rhamnose`)
auc = r %>% pROC::auc()
auc_ci <- pROC::ci(r)
tb5 = data.frame(mt="Rhamnose",auc=auc,auc_h=auc_ci[3],auc_l=auc_ci[1],col="#e5965a")

r = pROC::roc(response$response,norm$`TG(17:1_17:1_24:6)`)
auc = r %>% pROC::auc()
auc_ci <- pROC::ci(r)
tb6 = data.frame(mt="TG(17:1_17:1_24:6)",auc=auc,auc_h=auc_ci[3],auc_l=auc_ci[1],col="#e5965a")


tb = rbind(tb1,tb2,tb3,tb4,tb5,tb6)
tb = tb %>% arrange(auc)
tb$mt = factor(tb$mt,levels = tb$mt,ordered = T)


pdf("../plot/top_sig_mt_AUC.pdf",width = 5,height = 3)
ggplot(tb)+
  geom_segment(aes(x=auc_l,xend=auc_h,y=mt,yend=mt,col=col))+
  geom_segment(aes(x=auc_l,xend=auc_l,y=as.numeric(mt)-0.1,yend=as.numeric(mt)+0.1,col=col))+
  geom_segment(aes(x=auc_h,xend=auc_h,y=as.numeric(mt)-0.1,yend=as.numeric(mt)+0.1,col=col))+
  geom_point(aes(x=auc,y=mt,col=col))+
  scale_color_identity()+
  theme_bw()
dev.off()
