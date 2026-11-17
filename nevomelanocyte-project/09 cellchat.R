####cellchat
#加载
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ROGUE)
library(dplyr)
library(tidyverse)
library(CellChat)
library(circlize)
library(future)

table(subset@meta.data$sub_celltype)
table(subset@meta.data$orig.ident)


expr <- GetAssayData(object = subset@assays$RNA, layer = "data")


cellchat <- createCellChat(expr,
                           meta=subset@meta.data,
                           group.by="sub_celltype")


CellChatDB <- CellChatDB.human


showDatabaseCategory(CellChatDB)


cellchat@DB
#CellChatDB <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB


cellchat <- subsetData(cellchat)


future::plan("multisession", workers = 15)


plan(sequential)


cellchat <- identifyOverExpressedGenes(cellchat)


cellchat <- identifyOverExpressedInteractions(cellchat)


cellchat <- projectData(cellchat, PPI.human)


cellchat <- computeCommunProb(cellchat, raw.use = T) 
cellchat <- filterCommunication(cellchat, min.cells = 10)


cellchat <- computeCommunProbPathway(cellchat)


df.net <- subsetCommunication(cellchat)
write.csv(df.net, "1Gene.csv", row.names = F)


df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "2Pathway.csv", row.names = F)


cellchat <- aggregateNet(cellchat)


qsave(cellchat, file = "cellchat.qs")


groupSize <- as.numeric(table(cellchat@idents))
png("NetVisual_overview_all_number.png", width = 8, height = 6, units = "in", res = 300)
par(xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
dev.off()
png("NetVisual_overview_all_weights.png", width = 8, height = 6, units = "in", res = 300)
par(xpd = TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()


dir_name <- "NetVisual_separate"
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}

mat <- cellchat@net$weight
groupSize <- as.numeric(table(cellchat@idents))

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  
  filename <- file.path(dir_name, paste0("NetVisual_overview_split_", rownames(mat)[i], ".png"))
  
  
  png(filename, width = 6, height = 5, units = "in", res = 300)
  par(xpd = TRUE)
  
 
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
  
  
  dev.off()
}





mypathways <- cellchat@netP$pathways
mypathways <- c("LAMININ","COLLAGEN", "SPP1", "FN1", "CD99","PTN","THBS","TENASCIN","TGFb")
mypathways <- c("LAMININ","COLLAGEN","THBS")
mypathways


for (pathways.show in mypathways) {
  png(paste0("NetVisual_pathways_circle_", pathways.show, ".png"), 
      width = 2000, height = 1600, res = 300)  # 高分辨率 PNG
  par(xpd = TRUE)  
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
}




for(pathways.show in mypathways){
  png(paste0("NetVisual_pathways_heatmap_", pathways.show, ".png"),
      width = 2000, height = 1600, res = 300, units = "px")  # 高分辨率 PNG
  p <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  print(p)  # 
  dev.off()
}


#信号通路内的配体-受体
# 创建保存目录
dir.create("Pathways", showWarnings = FALSE)

for(pathways.show in mypathways){
  
  # 1️⃣ 配体-受体贡献度图
  png(paste0("Pathways/", pathways.show, "_contribution.png"),
      width = 2400, height = 1800, res = 300, units = "px")
  print(netAnalysis_contribution(cellchat, signaling = pathways.show))
  dev.off()
  
  # 提取富集的配体-受体对
  pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)$interaction_name
  
  # 2️⃣ 网络图(circle)
  for(LR.show in pairLR){
    png(paste0("Pathways/", pathways.show, "_", LR.show, "_circle.png"),
        width = 2400, height = 1800, res = 300, units = "px")
    netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
    dev.off()
  }
  
  # # 3️⃣ 网络图(chord)
  # for(LR.show in pairLR){
  #   png(paste0("Pathways/", pathways.show, "_", LR.show, "_chord.png"),
  #       width = 2400, height = 1800, res = 300, units = "px")
  #   netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  #   dev.off()
  # }
  # 
}

levels(cellchat@idents)

cellchat2 <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

for(pathways.show in mypathways){
  png(paste0("SignalingRole_", pathways.show, ".png"),
      width = 2400, height = 1000, res = 300, units = "px")
  
  netAnalysis_signalingRole_network(cellchat2,
                                    signaling = pathways.show,
                                    width = 8, height = 2.5, font.size = 10)
  dev.off()
}



#### CCI----
p <- netVisual_bubble(cellchat, 
                      sources.use = 1:length(levels(cellchat@idents)), 
                      targets.use = 1:length(levels(cellchat@idents)), 
                      remove.isolate = FALSE)
ggsave("10CCI_all.png", plot = p, 
       width = 8, height = 30, units = "in", dpi = 300, limitsize = FALSE)




send <- "TGFβ_responsive_Melanocytes"    
p <- netVisual_bubble(cellchat, 
                      sources.use = send, 
                      targets.use = 1:5, 
                      signaling = c("LAMININ","COLLAGEN", "FN1", "CD99","PTN","THBS","TENASCIN"), #可以按照通路筛选
                      #pairLR.use =  pairLR.use, #指定配体-受体
                      remove.isolate = FALSE)

num_LR_pairs <- length(unique(p$data$interaction_name_2))
dynamic_height <- 2 + (num_LR_pairs * 0.25)
if(dynamic_height < 4) { dynamic_height <- 4 }
if(dynamic_height > 20) { dynamic_height <- 20 }

file_name <- paste0("CCI_subcell_", send, ".png")
ggsave(file_name, plot = p, 
       width = 8, height = dynamic_height, units = "in", dpi = 300, limitsize = FALSE)





p <- plotGeneExpression(cellchat, signaling = "COLLAGEN")
ggsave("14GeneExpression_violin_sig.png", p, width = 10, height = 9, units = "in", dpi = 300, limitsize = F)
p <- plotGeneExpression(cellchat, signaling = "COLLAGEN", enriched.only = FALSE)
ggsave("15GeneExpression_violin_all.png", p, width = 10, height = 9, units = "in", dpi = 300, limitsize = F)

####
# 
fibro_subset <- subset(seurat, idents = "Fibroblasts")


fibro_subset$cellchat_anno <- "Fibroblasts"


subset$cellchat_anno <- subset$sub_celltype 


combined_chat <- merge(fibro_subset, y = subset, 
                       add.cell.ids = c("Fibro", "Mel"), 
                       project = "Nevus_Microenvironment")


table(combined_chat$cellchat_anno)



combined_chat <- JoinLayers(combined_chat)


table(combined_chat@meta.data$cellchat_anno)

expr <- GetAssayData(object = combined_chat, assay = "RNA", layer = "data")

cellchat <- createCellChat(expr,
                           meta = combined_chat@meta.data,
                           group.by = "cellchat_anno")


CellChatDB <- CellChatDB.human


showDatabaseCategory(CellChatDB)


cellchat@DB

cellchat@DB <- CellChatDB


cellchat <- subsetData(cellchat)


future::plan("multisession", workers = 15)


plan(sequential)


cellchat <- identifyOverExpressedGenes(cellchat)


cellchat <- identifyOverExpressedInteractions(cellchat)




cellchat <- computeCommunProb(cellchat, raw.use = T) 
cellchat <- filterCommunication(cellchat, min.cells = 10)


cellchat <- computeCommunProbPathway(cellchat)


df.net <- subsetCommunication(cellchat)
write.csv(df.net, "1Gene.csv", row.names = F)


df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "2Pathway.csv", row.names = F)


cellchat <- aggregateNet(cellchat)


qsave(cellchat, file = "cellchat.qs")


groupSize <- as.numeric(table(cellchat@idents))
png("NetVisual_overview_all_number.png", width = 8, height = 6, units = "in", res = 300)
par(xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
dev.off()
png("NetVisual_overview_all_weights.png", width = 8, height = 6, units = "in", res = 300)
par(xpd = TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()


mat <- cellchat@net$weight
groupSize <- as.numeric(table(cellchat@idents))

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  # 定义文件名
  filename <- paste0("NetVisual_overview_split_", rownames(mat)[i], ".png")
  
  # 打开 PNG 设备
  png(filename, width = 6, height = 5, units = "in", res = 300)
  par(xpd = TRUE)
  
  # 绘图
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
  
  # 关闭设备
  dev.off()
}





mypathways <- cellchat@netP$pathways
mypathways <- c("LAMININ","COLLAGEN", "SPP1", "FN1", "CD99","PTN","THBS","TENASCIN")
mypathways <- c("COLLAGEN")


for (pathways.show in mypathways) {
  png(paste0("NetVisual_pathways_circle_", pathways.show, ".png"), 
      width = 2000, height = 1600, res = 300)  # 高分辨率 PNG
  par(xpd = TRUE)  # 允许绘图元素超出边界
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  dev.off()
}




for(pathways.show in mypathways){
  png(paste0("NetVisual_pathways_heatmap_", pathways.show, ".png"),
      width = 2000, height = 1600, res = 300, units = "px")  # 高分辨率 PNG
  p <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  print(p)  # 
  dev.off()
}



dir.create("Pathways", showWarnings = FALSE)

for(pathways.show in mypathways){
  
  # 1️⃣ 配体-受体贡献度图
  png(paste0("Pathways/", pathways.show, "_contribution.png"),
      width = 2400, height = 1800, res = 300, units = "px")
  print(netAnalysis_contribution(cellchat, signaling = pathways.show))
  dev.off()
  
  # 提取富集的配体-受体对
  pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)$interaction_name
  
  # 2️⃣ 网络图(circle)
  for(LR.show in pairLR){
    png(paste0("Pathways/", pathways.show, "_", LR.show, "_circle.png"),
        width = 2400, height = 1800, res = 300, units = "px")
    netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
    dev.off()
  }
  
  # # 3️⃣ 网络图(chord)
  # for(LR.show in pairLR){
  #   png(paste0("Pathways/", pathways.show, "_", LR.show, "_chord.png"),
  #       width = 2400, height = 1800, res = 300, units = "px")
  #   netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  #   dev.off()
  # }
  # 
}

levels(cellchat@idents)


p <- netVisual_bubble(cellchat, 
                      sources.use = 1:length(levels(cellchat@idents)), 
                      targets.use = 1:length(levels(cellchat@idents)), 
                      remove.isolate = FALSE)
ggsave("10CCI_all.png", plot = p, 
       width = 8, height = 30, units = "in", dpi = 300, limitsize = FALSE)



p <- netVisual_bubble(cellchat, 
                      sources.use = "Proliferative_Melanocytes", 
                      targets.use = 1:5, 
                      signaling = c("LAMININ","COLLAGEN", "FN1", "CD99","PTN","THBS","TENASCIN"), #可以按照通路筛选
                      #pairLR.use =  pairLR.use, 
                      remove.isolate = FALSE)
ggsave("CCI_subcell_Proliferative_Melanocytes.png", plot = p, 
       width = 8, height = 12, units = "in", dpi = 300, limitsize = FALSE)





p <- plotGeneExpression(cellchat, signaling = "COLLAGEN")
ggsave("14GeneExpression_violin_sig.png", p, width = 10, height = 9, units = "in", dpi = 300, limitsize = F)
p <- plotGeneExpression(cellchat, signaling = "COLLAGEN", enriched.only = FALSE)
ggsave("15GeneExpression_violin_all.png", p, width = 10, height = 9, units = "in", dpi = 300, limitsize = F)


cellchat2 <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

for(pathways.show in mypathways){
  png(paste0("SignalingRole_", pathways.show, ".png"),
      width = 2400, height = 1000, res = 300, units = "px")
  
  netAnalysis_signalingRole_network(cellchat2,
                                    signaling = pathways.show,
                                    width = 8, height = 2.5, font.size = 10)
  dev.off()
}






my_colors <- c(
  "Proliferative_Melanocytes" = "#BF812D",
  "TGFβ_responsive_Melanocytes" = "#DFC27D",
  "MX1+_MX2+_Melanocytes" = "#80CDC1",
  "ECM_Melanocytes" = "#01665E",
  "Fibroblasts" = "#d73027"
   
)


groupSize <- as.numeric(table(cellchat@idents))

pdf("NetVisual_overview_all_number.pdf", width = 7, height = 7) 

par(xpd = TRUE, mar = c(1, 1, 3, 1)) 
netVisual_circle(
  cellchat@net$count, 
  vertex.weight = groupSize, 
  weight.scale = TRUE, 
  label.edge = FALSE, 
  color.use = my_colors,      # my_colors_matched
  edge.width.max = 8,         
  title.name = "Number of interactions"
)
dev.off()



pdf("NetVisual_overview_all_weights.pdf", width = 7, height = 7)
par(xpd = TRUE, mar = c(1, 1, 3, 1))
netVisual_circle(
  cellchat@net$weight, 
  vertex.weight = groupSize, 
  weight.scale = TRUE, 
  label.edge = FALSE, 
  color.use = my_colors,      # my_colors_matched
  edge.width.max = 8,
  title.name = "Interaction weights/strength"
)
dev.off()


dir_name <- "NetVisual_separate_PDF"
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}

mat <- cellchat@net$weight
groupSize <- as.numeric(table(cellchat@idents))


weight.max <- max(mat) 


for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  
  sender_name <- rownames(mat)[i]
  
  
  filename <- file.path(dir_name, paste0("NetVisual_split_", sender_name, ".pdf"))
  pdf(filename, width = 6, height = 6)
  
  
  par(xpd = TRUE, mar = c(1, 1, 3, 1))
  
  
  netVisual_circle(
    mat2, 
    vertex.weight = groupSize, 
    weight.scale = TRUE, 
    edge.weight.max = weight.max, 
    title.name = paste("Sender:", sender_name), 
    color.use = my_colors                       
  )
  
  dev.off()
}


mypathways <- cellchat@netP$pathways

mypathways <- c("LAMININ", "COLLAGEN","FN1" ,"THBS", "TENASCIN")


dir.create("Pathway_Specific_PDF", showWarnings = FALSE)

for (pathways.show in mypathways) {
  
  
  pdf(file.path("Pathway_Specific_PDF", paste0("NetVisual_circle_", pathways.show, ".pdf")), 
      width = 8, height = 8)
  
  
  circlize::circos.clear()
  par(xpd = TRUE, mar = c(2, 2, 3, 2)) 
  
  
  netVisual_aggregate(
    cellchat, 
    signaling = pathways.show, 
    layout = "circle", 
    color.use = my_colors,           
  )
  
  dev.off()
}





library(NMF)
library(ggalluvial)


cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 3)

pattern_colors <- c(
  "Pattern 1" = "#A6CEE3",  # 浅蓝
  "Pattern 2" = "#1F78B4",  # 深蓝
  "Pattern 3" = "#B2DF8A"   # 浅绿
)

signaling <- c("LAMININ","COLLAGEN","FN1","THBS","TENASCIN","HSPG")


pdf("CellChat_incoming_RiverPlot.pdf", width = 8, height = 10)
p_river <- netAnalysis_river(
  cellchat, 
  pattern = "incoming",   #"outgoing" or "incoming"
  # signaling = signaling,
  color.use = my_colors,
  color.use.pattern = pattern_colors,
  color.use.signaling = "grey80",
  do.order = T
)
print(p_river)
dev.off()




library(ComplexHeatmap)

pdf("CellChat_Global_Role_Heatmap.pdf", width = 15, height = 8)

signaling <- c("LAMININ","COLLAGEN","FN1","THBS","TENASCIN","HSPG")


ht_out <- netAnalysis_signalingRole_heatmap(
  cellchat2, 
  signaling = signaling,
  pattern = "outgoing", 
  color.use = my_colors,
  title = "Outgoing Signaling (Sender)"
)


ht_in <- netAnalysis_signalingRole_heatmap(
  cellchat2, 
  signaling = signaling,
  pattern = "incoming", 
  color.use = my_colors,
  title = "Incoming Signaling (Receiver)"
)


draw(ht_out + ht_in, ht_gap = unit(0.5, "cm"))

dev.off()



library(ggplot2)
library(viridis)

send <- "Fibroblasts"

p_bubble <- netVisual_bubble(
  cellchat, 
  sources.use = send, 
  targets.use = 1:5, 
  signaling = c("LAMININ","COLLAGEN", "FN1", "CD99","PTN","THBS","TENASCIN","HSPG"), 
  remove.isolate = TRUE
)


p_bubble_advanced <- p_bubble + 
  
  scale_color_viridis_c(option = "magma", direction = -1, na.value = NA) + 
  
  theme_bw() + 
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3), # 极细的网格线
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 11),
    axis.text.y = element_text(color = "black", face = "italic", size = 10), # 🌟 基因名必须斜体！
    strip.background = element_rect(fill = "#2C3E50"), # 顶部标签的背景色
    strip.text = element_text(color = "white", face = "bold")
  )


num_LR_pairs <- length(unique(p_bubble$data$interaction_name_2))
dynamic_height <- 2 + (num_LR_pairs * 0.25)
dynamic_height <- max(4, min(dynamic_height, 20)) 


file_name <- paste0("CCI_subcell_", send, "_Advanced.pdf")
ggsave(file_name, plot = p_bubble_advanced, 
       width = 7, height = dynamic_height, dpi = 300, limitsize = FALSE)



dir.create("Pathways", showWarnings = FALSE)

mypathways <- c("COLLAGEN")

for(pathways.show in mypathways){
  
  # 1️⃣ 配体-受体贡献度图
  pdf(paste0("Pathways/", pathways.show, "_contribution.pdf"),
      width = 6, height = 10)
  print(netAnalysis_contribution(cellchat, signaling = pathways.show))
  dev.off()
  
  # 提取富集的配体-受体对
  pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)$interaction_name
  
  # 2️⃣ 网络图(circle)
  for(LR.show in pairLR){
    pdf(paste0("Pathways/", pathways.show, "_", LR.show, "_circle.pdf"),
        width = 7, height = 7)
    netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle", color.use = my_colors)
    dev.off()
  }
  
}
