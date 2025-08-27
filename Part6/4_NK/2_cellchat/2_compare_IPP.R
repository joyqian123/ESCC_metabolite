rm(list = ls())
library(copykat)
#Main Celltype Annotation----
# rm(list = ls())
library(remotes)
# memory.limit()
# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(future)
setwd("~/help_for_others/DYQ/scRNA//")

library(harmony)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(copykat)
library(Rtsne)
library(CellChat)
library(ggsci)

load("./15_NK/cellchat/IPP/cellchat_list.rdata")

cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use = pal_d3("category10")(10))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",
                           color.use = pal_d3("category10")(10))
pdf("./15_NK/cellchat/IPP/compareInteraction.pdf",width = 6,height = 3)
gg1 + gg2
dev.off()

order = cellchat@idents[[1]] %>% unique() %>% as.character() %>% sort()
cellchat@idents[[1]] <- factor(cellchat@idents[[1]],
                               levels = order,ordered = TRUE)
cellchat@idents[[2]] <- factor(cellchat@idents[[2]],
                               levels = order,ordered = TRUE)
cellchat@idents[[3]] <- factor(cellchat@idents[[3]],
                               levels = order,ordered = TRUE)

col_use = c("#bdc3d2","#88c4e8","#0074b3",
            "#e5ce81","#4eb69e","#aedacb","#7c6aa4",
            "#f47720","#f7bfa4","#e6966a",
            "#982b2d","#de5a69","#ee9d9f","#fccfc9")

pdf("./15_NK/cellchat/IPP/netVisual_diffInteraction.pdf",width = 8,height = 4)
par(mfrow = c(1,2), xpd=TRUE)
p1 = netVisual_diffInteraction(cellchat, weight.scale = T,color.use = col_use,
                               vertex.size.max = 9,vertex.label.cex = 0.8)
p2 = netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",color.use = col_use,
                               vertex.size.max = 9,vertex.label.cex = 0.8)
# p1 + p2
dev.off()

pdf("./15_NK/cellchat/IPP/netVisual_diffInteraction_to_NK.pdf",width = 8,height = 4)
par(mfrow = c(1,2), xpd=TRUE)
p1 = netVisual_diffInteraction(cellchat, weight.scale = T,color.use = col_use,targets.use = order[c(11:14)],
                               vertex.size.max = 9,vertex.label.cex = 0.8)
p2 = netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",color.use = col_use,targets.use = order[c(11:14)],
                               vertex.size.max = 9,vertex.label.cex = 0.8)
# p1 + p2
dev.off()




gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
pdf("./15_NK/cellchat/IPP/netVisual_heatmap.pdf",width = 10,height = 6)
gg1 + gg2
dev.off()


i=1

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", 
               # sources.use = c("Mono_macro"),
               targets.use = c("NK_Cd11b+Cd27-_Klrc2"), 
               stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", 
               # sources.use = c("Mono_macro"),
               targets.use = c("NK_Cd11b+Cd27-_Klrc2"), 
               stacked = F, do.stat = TRUE)
pdf("./15_NK/cellchat/IPP/all_to_NK_Cd11b+Cd27-_Klrc2.pdf",width = 10,height = 6)
gg1 + gg2
dev.off()


pathways.show <- c("CCL") 

# weight.max <- getMaxWeight(cellchat_list, slot.name = c("netP"), attribute = pathways.show)
pdf("./15_NK/cellchat/IPP/CCL.pdf",width = 15,height =15)
par(mfrow = c(1,2))
netVisual_aggregate(cellchat_list[[1]], signaling = pathways.show, layout = "chord",
                    # sources.use = order[c(1:10)],
                    targets.use = c("NK_Cd11b+Cd27-_Klrc2"),
                    cell.order = order,color.use = col_use,
                    signaling.name = paste(pathways.show, names(cellchat_list)[1]))
netVisual_aggregate(cellchat_list[[2]], signaling = pathways.show, layout = "chord",
                    # sources.use = order[c(1:10)],
                    targets.use = c("NK_Cd11b+Cd27-_Klrc2"),
                    cell.order = order,color.use = col_use,
                    signaling.name = paste(pathways.show, names(cellchat_list)[2]))

# }
dev.off()


pdf("./15_NK/cellchat/IPP/CCL_detail.pdf",width = 7,height = 4)
signaling = data.frame(pathway_name="CCL")
gg1 <- netVisual_bubble(cellchat, signaling= signaling,
                        sources.use = order,
                        targets.use = c("NK_Cd11b+Cd27-_Klrc2"),
                        comparison = c(1,2),
                        angle.x = 45, remove.isolate = T)
gg1
dev.off()
