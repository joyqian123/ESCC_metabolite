rm(list = ls())
setwd("~/help_for_others/DYQ/output/7_JM_1//")
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

class_ID_1 = mt_info %>% dplyr::select(mt,`Class I`) %>% dplyr::rename(class_1=colnames(.)[2])
class_ID_2 = mt_info %>% dplyr::select(mt,`Class II`) %>% dplyr::rename(class=colnames(.)[2]) %>% 
  dplyr::filter(!class=="Hormones and hormone related compounds") %>% 
  dplyr::filter(!class=="Organic acid and Its derivatives")

final = final %>% dplyr::rename(mt=metabo) %>% inner_join(class_ID_2,by = "mt") %>% inner_join(class_ID_1,by = "mt")


library(ggpubr)
pdf("./plot/plot_R1_R2.pdf",width = 5,height = 4)
ggplot(final,aes(x=coef_R1,y=coef_R2))+
  geom_point(aes(col=label,size=label))+
  xlim(-0.6,0.6)+
  ylim(-2,2)+
  geom_hline(yintercept = 0,lty=2)+
  geom_vline(xintercept = 0,lty=2)+
  geom_smooth(method = "lm",lwd=1,lty=2,col="#C71585")+
  stat_cor(method = "pearson",
           label.x.npc = 0,
           label.y.npc = 1) +
  scale_color_manual(values = c("grey","#6A5ACD"))+
  scale_size_manual(values = c(0.2,1))+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()

final_need =final %>% dplyr::filter(label==1)


final_need_risk = final_need %>% dplyr::filter(coef_R1>0.1,coef_R2>0.1) %>% pull(mt)
final_need_protect = final_need %>% dplyr::filter(coef_R1<(-0.1),coef_R2<(-0.1)) %>% pull(mt)


#########################
load("~/help_for_others/DYQ//reference/kegg_ID_sum.rdata")
library(clusterProfiler)
risk_up_res = enricher(final_need_risk,TERM2GENE = ID_sum[,c(2,1)])
risk_up = risk_up_res@result %>% dplyr::filter(pvalue<0.05) %>% dplyr::mutate(group="risk")

protect_up_res = enricher(final_need_protect,TERM2GENE = ID_sum[,c(2,1)])
protect_up = protect_up_res@result %>% dplyr::filter(pvalue<0.05) %>% dplyr::mutate(group="protect")


library(ComplexHeatmap)
plot=final_need %>% dplyr::filter(mt %in% c(final_need_protect,final_need_risk)) %>% arrange(class,coef_R1+coef_R2) %>% 
  column_to_rownames("mt") %>% dplyr::select(coef_R1,coef_R2) 
plot_anno = final_need %>% dplyr::filter(mt %in% c(final_need_protect,final_need_risk)) %>%  arrange(class,coef_R1+coef_R2) %>% 
  column_to_rownames("mt") %>% dplyr::select(class_1)


plot_scale = scale(plot) %>% data.frame()
row_ha = rowAnnotation(class=plot_anno$class_1)



final_enrich = rbind(risk_up,protect_up) %>% separate(.,col="GeneRatio",into=c("A","B"),sep="/")
final_enrich$A = as.numeric(final_enrich$A)
final_enrich$B = as.numeric(final_enrich$B)  
final_enrich = final_enrich %>% dplyr::mutate(GeneRatio=A/B)
final_enrich = final_enrich %>% arrange(group,GeneRatio)
sel_ID = final_enrich %>% distinct(geneID,.keep_all = TRUE) %>% pull(ID)
# final_enrich = final_enrich %>% dplyr::filter(ID %in% sel_ID)
final_enrich$ID = factor(final_enrich$ID,levels = final_enrich$ID,ordered = TRUE)

pdf("./plot/enricher_plot.pdf",width = 8,height = 5)
ggplot(final_enrich)+
  geom_point(aes(x=GeneRatio,y=ID,size=-log(p.adjust),col=-log(p.adjust)))+
  facet_grid(~group)+
  scale_color_gradient2(low = "#E6E6FA",mid = "#6495ED",high = "#9932CC")+
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_blank())
dev.off()



################################################################
################################################################
en_res = enricher(final_need_risk,TERM2GENE=class_ID_2[,c(2,1)])
en_res_class_II = en_res@result

plot = separate(data = en_res_class_II,col = "GeneRatio",sep = "/",into = c("A","B"),remove = FALSE) %>% 
  dplyr::mutate(A=as.numeric(A),B=as.numeric(B)) %>% dplyr::mutate(r=A/B)

myPalette <- colorRampPalette(c("blue","green","yellow","red"))
gradientColors <- myPalette(1000)

library(ggrepel)

pdf("./plot/risk_mt_enrichment.pdf",width = 4,height = 3)
ggplot(plot,aes(x=r,y=-log(pvalue)))+
  geom_point(aes(size=-log(pvalue),fill=-log(pvalue)),shape=21)+
  scale_fill_gradientn(colours = gradientColors)+
  geom_text_repel(aes(label=ifelse(pvalue<0.1,ID,"")))+
  geom_hline(yintercept = -log(0.1),lty=2)+
  ggtitle('risk_mt enrichment')+
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



en_res = enricher(final_need_protect,TERM2GENE=class_ID_2[,c(2,1)])
en_res_class_II = en_res@result

plot = separate(data = en_res_class_II,col = "GeneRatio",sep = "/",into = c("A","B"),remove = FALSE) %>% 
  dplyr::mutate(A=as.numeric(A),B=as.numeric(B)) %>% dplyr::mutate(r=A/B)

myPalette <- colorRampPalette(c("blue","green","yellow","red"))
gradientColors <- myPalette(1000)

library(ggrepel)

pdf("./plot/protect_mt_enrichment.pdf",width = 4,height = 3)
ggplot(plot,aes(x=r,y=-log(pvalue)))+
  geom_point(aes(size=-log(pvalue),fill=-log(pvalue)),shape=21)+
  scale_fill_gradientn(colours = gradientColors)+
  geom_text_repel(aes(label=ifelse(pvalue<0.1,ID,"")))+
  geom_hline(yintercept = -log(0.1),lty=2)+
  ggtitle('protect_mt enrichment')+
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





##############################################
library(ComplexHeatmap)
library(circlize)
plot=final_need %>% dplyr::filter(mt %in% c(final_need_protect,final_need_risk)) %>% 
  dplyr::filter(class %in% c("SM","TG","Cer-NS")) %>% 
  arrange(class,coef_R1+coef_R2) %>% 
  column_to_rownames("mt") %>% dplyr::select(coef_R1,coef_R2) 
plot_anno = final_need %>% dplyr::filter(mt %in% c(final_need_protect,final_need_risk)) %>% 
  dplyr::filter(class %in% c("SM","TG","Cer-NS")) %>% 
  arrange(class,coef_R1+coef_R2) %>% 
  column_to_rownames("mt") %>% dplyr::select(class)

# plot_scale = scale(plot) %>% data.frame()
row_ha = rowAnnotation(class=plot_anno$class,
                       col = list(class = c("SM" = "#FF1493", "TG" = "#8B008B", "Cer-NS" = "#DDA0DD")))
col_fun = colorRamp2(breaks = c(-1,0,1),colors = c("#6A5ACD","white","#C71585"))
pdf("./plot//heatmap_R1_R2_risk.pdf",width = 3.5,height = 7)
Heatmap(plot,cluster_rows = FALSE,cluster_columns = FALSE,
        left_annotation = row_ha,col = col_fun,row_names_gp = gpar(fontsize=7),
        row_names_side = "right",border = "black")
dev.off()





library(ComplexHeatmap)
plot=final_need %>% dplyr::filter(mt %in% c(final_need_protect,final_need_risk)) %>% 
  dplyr::filter(class %in% c("PI","LPC","LPE")) %>% 
  arrange(class,coef_R1+coef_R2) %>% 
  column_to_rownames("mt") %>% dplyr::select(coef_R1,coef_R2) 
plot_anno = final_need %>% dplyr::filter(mt %in% c(final_need_protect,final_need_risk)) %>% 
  dplyr::filter(class %in% c("PI","LPC","LPE")) %>% 
  arrange(class,coef_R1+coef_R2) %>% 
  column_to_rownames("mt") %>% dplyr::select(class)

# plot_scale = scale(plot) %>% data.frame()
row_ha = rowAnnotation(class=plot_anno$class,
                       col = list(class = c("PI" = "#000080", "PC" = "#0000CD", "PA" = "#4169E1",
                                            "LPE" = "#6495ED", "LPC-O" = "#87CEFA", "LPC" = "#00BFFF")))
col_fun = colorRamp2(breaks = c(-1,0,1),colors = c("#6A5ACD","white","#C71585"))
pdf("./plot//heatmap_R1_R2_protect.pdf",width = 3.5,height = 7)
Heatmap(plot,cluster_rows = FALSE,cluster_columns = FALSE,
        left_annotation = row_ha,col = col_fun,
        row_names_gp = gpar(fontsize=7),
        row_names_side = "right",border = "black")
dev.off()




##########################################################
##########################################################
plot = final %>% dplyr::select(mt,coef_R1,coef_R2) %>% column_to_rownames("mt")
library(RobustRankAggreg)
method_choose="mean"
glist<-list()
i=2
for (i in 1:ncol(plot)){
  data<-plot %>% dplyr::select(any_of(i)) %>% dplyr::rename(mt=colnames(.)[1]) %>% arrange(mt)
  glist[[i]]<-rownames(data)
  names(glist)[i]<-colnames(plot)[i]
}
library(RobustRankAggreg)
r<-rankMatrix(glist,full = TRUE)
RRA<-aggregateRanks(glist,rmat=r,
                    method=method_choose,full=TRUE,exact=F,topCutoff=NA)
# rho_neg<-data.frame(apply(r_neg,1,rhoScores),check.names = FALSE)
RRA$p.adj<-p.adjust(RRA$Score,method = "BH")
RRA$rank<-rank(RRA$Score)


target_ID = class_ID_1 %>% dplyr::filter(class_1 %in% c("SP","GP")) %>% dplyr::rename(class=class_1)
target_ID_2 = class_ID_2 %>% dplyr::filter(grepl("Cer",class)) %>% dplyr::mutate(class="Cer")
target_ID_3 = class_ID_2 %>% dplyr::filter(grepl("Cer",class))
target_ID_4 = class_ID_2 %>% dplyr::filter(class %in% c("SM","SPH"))
target_ID_5 = class_ID_2 %>% dplyr::filter(class %in% c("LPC","LPE","LPI","LPS","LPA","LPG","LPC-O","LPE-P",
                                                        "PC","PE","PI","PS","PA","PG","PC-O","PE-P","PE-O"))
all_target_1 = rbind(target_ID,target_ID_2,target_ID_4)
all_target_2 = rbind(target_ID_3,target_ID_4,target_ID_5)

geneList= final$coef_R1
names(geneList) = final$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
library(GseaVis)
library(GSEABase)
egmt <- GSEA(geneList, TERM2GENE=all_target_1[,c(2,1)], 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_coef_R1 = egmt@result %>% dplyr::mutate(type="R1")
p1=gseaNb(object = egmt,
          geneSetID = c("SP","GP","SM","Cer"),
          newGsea = T,
          addPval = T,
          pvalX = 0.5,pvalY = 0.75,
          pCol = 'black',
          pHjust = 1,
          subPlot = 2)

geneList= final$coef_R2
names(geneList) = final$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=all_target_1[,c(2,1)], 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_coef_R2 = egmt@result %>% dplyr::mutate(type="R2")
p2=gseaNb(object = egmt,
          geneSetID = c("SP","GP","SM","Cer"),
          newGsea = T,
          addPval = T,
          pvalX = 0.5,pvalY = 0.75,
          pCol = 'black',
          pHjust = 1,
          subPlot = 2)

geneList= (RRA$rank-median(RRA$rank))/(max(RRA$rank)-min(RRA$rank))
names(geneList) = RRA$Name
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=all_target_1[,c(2,1)], 
             minGSSize = 0,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_RRA = egmt@result %>% dplyr::mutate(type="RRA")
p3=gseaNb(object = egmt,
          geneSetID = c("SP","GP","SM","Cer"),
          newGsea = T,
          addPval = T,
          pvalX = 0.5,pvalY = 0.75,
          pCol = 'black',
          pHjust = 1,
          subPlot = 2)


pdf("./plot/GSEA_SP_GP_SM_Cer.pdf",width = 12,height = 6)
cowplot::plot_grid(p1,p2,p3,nrow = 1)
dev.off()



#######################################
set.seed(123)
geneList= final$coef_R1
names(geneList) = final$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=all_target_2[,c(2,1)], 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_coef_R1 = egmt@result %>% dplyr::mutate(type="R1")

geneList= final$coef_R2
names(geneList) = final$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=all_target_2[,c(2,1)], 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_coef_R2 = egmt@result %>% dplyr::mutate(type="R2")

geneList= (RRA$rank-median(RRA$rank))/(max(RRA$rank)-min(RRA$rank))
names(geneList) = RRA$Name
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=all_target_2[,c(2,1)], 
             minGSSize = 0,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_RRA = egmt@result %>% dplyr::mutate(type="RRA")


final_egmt = rbind(egmt_res_coef_R1,egmt_res_coef_R2,egmt_res_RRA)


final_plot = final_egmt %>% dplyr::select(ID,NES,pvalue,type)
final_NES = final_plot %>% pivot_wider(id_cols = ID,names_from = type,values_from = NES) %>% 
  dplyr::mutate(R1=as.numeric(R1),R2=as.numeric(R2),RRA=as.numeric(RRA)) %>% arrange(desc(RRA)) %>%  
  column_to_rownames("ID") %>% as.matrix()
rownames(final_NES) = factor(rownames(final_NES),levels = rownames(final_NES),ordered = T)
final_p = final_plot %>% pivot_wider(id_cols = ID,names_from = type,values_from = pvalue) %>% 
  dplyr::mutate(R1=as.numeric(R1),R2=as.numeric(R2),RRA=as.numeric(RRA)) %>% 
  dplyr::mutate(ID=factor(ID,levels = rownames(final_NES),ordered = T)) %>% 
  arrange(ID) %>% 
  dplyr::mutate(R1=ifelse(R1<0.001,"***",ifelse(R1<0.01,"**",ifelse(R1<0.05,"*","")))) %>% 
  dplyr::mutate(R2=ifelse(R2<0.001,"***",ifelse(R2<0.01,"**",ifelse(R2<0.05,"*","")))) %>% 
  dplyr::mutate(RRA=ifelse(RRA<0.001,"***",ifelse(RRA<0.01,"**",ifelse(RRA<0.05,"*","")))) %>% 
  column_to_rownames("ID") %>% as.matrix()


library(ComplexHeatmap)
pdf("./plot/GSEA_subtype.pdf",width = 4,height = 6)
Heatmap(final_NES,cluster_columns = FALSE,cluster_rows = FALSE,border = TRUE,
        rect_gp = gpar(col = "white"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(final_p[i, j], x, y, gp = gpar(fontsize = 10))
        })
dev.off()



###############################################
hit_color = c("SM" = "#FF1493", "TG" = "#8B008B", 
              "Cer-NS" = "#FF00FF","Cer-NP" = "#FF00FF",
              "Cer-AS" = "#FF00FF","Cer-AP" = "#FF00FF","HexCer-AP" = "#FF00FF",
              "PI" = "#000080","LPI" = "#000080", 
              "PC" = "#0000CD", "LPC-O" = "#0000CD", "LPC" = "#0000CD", 
              "PA" = "#4169E1","LPA" = "#4169E1",
              "PE" = "#6495ED","LPE" = "#6495ED","PE-O" = "#6495ED","PE-P" = "#6495ED","LPE-P" = "#6495ED") %>% data.frame() %>% 
  rownames_to_column("type") %>% dplyr::rename(col=colnames(.)[2])

plot_ls = lapply(1:nrow(hit_color),function(kk){
  
  hit = final %>% dplyr::filter(class %in% hit_color$type[kk])
  p=ggplot(final,aes(x=coef_R1,y=coef_R2))+
    geom_point(col="grey")+
    xlim(-0.6,0.6)+
    ylim(-2,2)+
    geom_hline(yintercept = 0,lty=2)+
    geom_vline(xintercept = 0,lty=2)+
    geom_density_2d(data = hit, color = hit_color$col[kk]) +
    # ggrepel::geom_text_repel(aes(x=coef_R1,y=coef_R2,label=ifelse(label==1,mt,NA)),min.segment.length = 0.5)+
    # scale_color_manual(values = c("grey","#6A5ACD"))+
    scale_size_manual(values = c(0.2,1))+
    ggtitle(hit_color$type[kk])+
    theme_bw() +
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank())

})


pdf("./plot/density/density_one_type.pdf",width = 18,height = 9)
wrap_plots(plot_ls,nrow = 3)
dev.off()

all_target_1=all_target_1 %>% dplyr::filter(!class == "SPH")
hit_color=c("#000080","#FF1493","#FF00FF","#FF1493")
kk=1
plot_ls = lapply(1:length(unique(all_target_1$class)),function(kk){
  
  mt_sel = all_target_1 %>% dplyr::filter(class %in% unique(all_target_1$class)[kk]) %>% pull(mt)
  hit = final %>% dplyr::filter(mt %in% mt_sel)
  p=ggplot(final,aes(x=coef_R1,y=coef_R2))+
    geom_point(col="grey")+
    xlim(-0.6,0.6)+
    ylim(-2,2)+
    geom_hline(yintercept = 0,lty=2)+
    geom_vline(xintercept = 0,lty=2)+
    geom_density_2d(data = hit, color = hit_color[kk]) +
    # ggrepel::geom_text_repel(aes(x=coef_R1,y=coef_R2,label=ifelse(label==1,mt,NA)),min.segment.length = 0.5)+
    # scale_color_manual(values = c("grey","#6A5ACD"))+
    scale_size_manual(values = c(0.2,1))+
    ggtitle(unique(all_target_1$class)[kk])+
    theme_bw() +
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank())
  
})
library(patchwork)
pdf("./plot/density/density_merge_type.pdf",width = 12,height = 3)
wrap_plots(plot_ls,nrow = 1)
dev.off()







##################################################################################
##################################################################################
all_rank = RRA %>% dplyr::rename(mt=Name) %>% inner_join(final[,c("mt","coef_R1","coef_R2","class")],by = "mt") %>% 
  dplyr::mutate(RRA=(RRA$rank-median(RRA$rank))/(max(RRA$rank)-min(RRA$rank)))
all_rank_TG = all_rank %>% dplyr::filter(class=="TG") %>% 
  dplyr::mutate(mt1=str_remove(mt,"^TG\\(")) %>% dplyr::mutate(mt1=str_remove(mt1,"\\)")) %>% 
  separate(.,col=mt1,into=c("A","B","C"),sep="_")%>% 
  separate(.,col=A,into=c("A1","A2"),remove = FALSE,sep=":")%>% 
  separate(.,col=B,into=c("B1","B2"),remove = FALSE,sep=":")%>% 
  separate(.,col=C,into=c("C1","C2"),remove = FALSE,sep=":") %>% 
  dplyr::mutate(A1=as.numeric(A1),A2=as.numeric(A2),B1=as.numeric(B1),B2=as.numeric(B2),C1=as.numeric(C1),C2=as.numeric(C2)) %>% 
  dplyr::mutate(total_C = A1+B1+C1,total_unsature=A2+B2+C2)

cor.test(all_rank_TG$total_unsature,all_rank_TG$RRA)


facet_means <- all_rank_TG %>%
  group_by(total_C, total_unsature) %>%
  summarise(mean_coef = mean(coef_R1, na.rm = TRUE)) %>%
  ungroup()
pdf("./plot/TG/coef_R1_facet.pdf",width = 10,height = 6)
ggplot(all_rank_TG)+
  geom_rect(data = facet_means, 
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = mean_coef), 
            alpha = 0.7) +  # 透明度控制背景色
  geom_boxplot(aes(x=coef_R1,y=1))+
  scale_fill_gradient2(breaks=c(-0.2,0,0.75),low = "#6A5ACD",mid = "white", high = "#C71585") + 
  facet_grid(total_C~total_unsature)+
  scale_x_continuous(breaks = c(-0.5,0,0.5))+
  theme_bw()+
  theme(strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
dev.off()

facet_means <- all_rank_TG %>%
  group_by(total_C, total_unsature) %>%
  summarise(mean_coef = mean(coef_R2, na.rm = TRUE)) %>%
  ungroup()
pdf("./plot/TG/coef_R2_facet.pdf",width = 10,height = 6)
ggplot(all_rank_TG)+
  geom_rect(data = facet_means, 
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = mean_coef), 
            alpha = 0.7) +  # 透明度控制背景色
  geom_boxplot(aes(x=coef_R2,y=1))+
  scale_fill_gradient2(breaks=c(-0.75,0,1),low = "#6A5ACD",mid = "white", high = "#C71585") + 
  facet_grid(total_C~total_unsature)+
  scale_x_continuous(breaks = c(-0.5,0,0.5))+
  theme_bw()+
  theme(strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
dev.off()

facet_means <- all_rank_TG %>%
  group_by(total_C, total_unsature) %>%
  summarise(mean_coef = mean(RRA, na.rm = TRUE)) %>%
  ungroup()
pdf("./plot/TG/coef_RRA_facet.pdf",width = 10,height = 6)
ggplot(all_rank_TG)+
  geom_rect(data = facet_means, 
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = mean_coef), 
            alpha = 0.7) +  # 透明度控制背景色
  geom_boxplot(aes(x=RRA,y=1))+
  scale_fill_gradient2(breaks=c(-0.5,0,0.5),low = "#6A5ACD",mid = "white", high = "#C71585") + 
  facet_grid(total_C~total_unsature)+
  scale_x_continuous(breaks = c(-0.5,0,0.5))+
  theme_bw()+
  theme(strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
dev.off()



# need_list = c("22:6","20:4","20:5","22:5","18:0","18:1","18:2","18:3","16:0","16:1","16:2")
need_list = unique(c(all_rank_TG$A,all_rank_TG$B,all_rank_TG$C))
kk=1
TG_spe_ls = lapply(1:length(need_list),function(nn){
  tb = all_rank_TG 
  tb$tmp = sapply(1:nrow(tb),function(kk){
    tmp = (need_list[nn] %in% tb$A[kk]) | (need_list[nn] %in% tb$B[kk])| (need_list[nn] %in% tb$C[kk])
    if(tmp==TRUE){
      return(need_list[nn])
    }else{
      return("others")
    }
  })
  need = tb %>% dplyr::filter(tmp==need_list[nn]) %>% dplyr::select(tmp,mt)
}) %>% do.call(rbind,.)




geneList= all_rank$coef_R1
names(geneList) = all_rank$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=TG_spe_ls, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_coef_R1 = egmt@result %>% dplyr::mutate(type="R1")

geneList= all_rank$coef_R2
names(geneList) = all_rank$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=TG_spe_ls, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_coef_R2 = egmt@result %>% dplyr::mutate(type="R2")

geneList= all_rank$RRA
names(geneList) = all_rank$mt
geneList=sort(geneList,decreasing = T)
head(geneList)
library(clusterProfiler)
egmt <- GSEA(geneList, TERM2GENE=TG_spe_ls, 
             minGSSize = 0,
             pvalueCutoff = 0.99,
             verbose=FALSE)
egmt_res_RRA = egmt@result %>% dplyr::mutate(type="RRA")



# ggplot(egmt_res_RRA)+
#   geom_point(aes(x=NES,y=-log(pvalue),size=-log(pvalue)))
myPalette <- colorRampPalette(c("blue","green","yellow","red"))
gradientColors <- myPalette(1000)

library(ggrepel)

# egmt_res_class_II_R = egmt_res_class_II %>% dplyr::filter(NES>0)

pdf("./plot/TG/chain_egmt.pdf",width = 5,height = 4)
ggplot(egmt_res_RRA,aes(x=NES,y=-log(pvalue)))+
  geom_point(aes(size=-log(pvalue),fill=abs(NES)),shape=21)+
  scale_fill_gradientn(colours = gradientColors)+
  geom_text_repel(aes(label=ifelse(pvalue<0.05,ID,"")),min.segment.length = 0.1)+
  geom_hline(yintercept = -log(0.05),lty=2)+
  geom_vline(xintercept = 0,lty=2)+
  ggtitle('TG_chain gsea')+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),  
        axis.ticks = element_line(color = "black"),
        # legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()
