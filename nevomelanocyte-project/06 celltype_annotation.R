library(qs)
library(scop)
library(scales)
library(scCustomize)
library(scRNAtoolVis)
library(Seurat)
library(dplyr)  


seurat <- qread("seurat_singlet.qs")


colnames(seurat@meta.data)


DimPlot(seurat,group.by = "seurat_clusters",label = T)




######## FindAllMarker----
table(Idents(seurat))
Idents(seurat) <- "seurat_clusters"
f = "Allcelltype_markers.csv"
if(!file.exists(f)){
  allmarkers <- FindAllMarkers(seurat, 
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25 
  )
  write.csv(allmarkers,file='Allcelltype_markers.csv')
} else {
  allmarkers = read.csv('Allcelltype_markers.csv',row.names = 1)
}
head(allmarkers)


library(tidyverse)
g = allmarkers %>% 
  group_by(cluster) %>% 
  top_n(5,wt = avg_log2FC) %>%  
  pull(gene) %>% 
  unique()

DotPlot(seurat,
        features = g,cols = "RdYlBu")+
  RotatedAxis()
ggsave("balloonplot.png",width = 24,height = 6)





genes_to_check = c("MITF","MLANA","TYR","PMEL","TYRP1","DCT", #Melanocytes
                   "KRT5","KRT14","TP63", #Keratinocytes
                   "DCN","COL1A1","COL3A1","PDGFRA", #Fibroblasts
                   "PECAM1","VWF","CLDN5","CDH5" #Endothelial
                   )  

library(stringr)  

genes_to_check=str_to_upper(genes_to_check) 
genes_to_check


for (gene in genes_to_check) {
  p <- FeaturePlot(seurat, features = gene) + ggtitle(gene)
  
  
  ggsave(
    filename = paste0(gene, ".png"),
    plot = p,
    width = 6, height = 5, units = "in", dpi = 300
  )
}


for (gene in genes_to_check) {
  p <- FeatureDimPlot(
    srt = seurat,
    features = gene,
    reduction = "UMAP",
    theme_use = "theme_blank"
  ) + ggtitle(gene)
  
  
  ggsave(
    filename = paste0(gene, "_scop.png"),
    plot = p,
    width = 6, height = 5, units = "in", dpi = 300
  )
}


Idents(seurat) <- "seurat_clusters"
table(Idents(seurat))
Idents(seurat) <- factor(Idents(seurat), 
                         levels = sort(as.numeric(levels(Idents(seurat)))))

p = DotPlot(seurat, 
            features = unique(genes_to_check),
            assay='RNA',cols = "RdYlBu"  )  + 
  coord_flip()+
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14,angle = 45, hjust = 1))
p
ggsave('DotPlot.png',plot = p, width = 8,height = 9)


markers_list <- list(
  Melanocytes = c("MITF","MLANA","TYR","PMEL","TYRP1","DCT"),
  Keratinocytes  = c("KRT5","KRT14","TP63"),
  Fibroblasts = c("DCN","COL1A1","COL3A1","PDGFRA"),
  Endothelial = c("PECAM1","VWF","CLDN5","CDH5")
)


genes_to_check <- lapply(markers_list, toupper)

dup <- names(table(unlist(genes_to_check)))[table(unlist(genes_to_check)) > 1]

genes_to_check <- lapply(genes_to_check, function(x) x[!x %in% dup])


Idents(seurat) <- factor(Idents(seurat),
                        levels = sort(as.numeric(levels(Idents(seurat)))))


p <- DotPlot(
  seurat,
  features = genes_to_check,   
  assay = "RNA",
  cols = "RdYlBu"
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1));p
w <- length(unique(unlist(genes_to_check))) / 5 + 6
ggsave("check_for_Nevi.png", plot = p, width = w)


cluster_ids <- sort(unique(seurat@meta.data$seurat_clusters))


celltype <- data.frame(
  ClusterID = cluster_ids,
  celltype = rep("Unknown", length(cluster_ids))
)


celltype[celltype$ClusterID %in% c(0,1,2,3,4,5,6,7,8,10,11,12,13,14,16,17 ),2]='Melanocytes'
celltype[celltype$ClusterID %in% c(19 ),2]='Keratinocytes'
celltype[celltype$ClusterID %in% c(9,18 ),2]='Fibroblasts'
celltype[celltype$ClusterID %in% c(15 ),2]='Endothelial'



table(celltype$celltype)

table(seurat@meta.data$seurat_clusters)

seurat@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  seurat@meta.data[which(seurat@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
#
table(seurat@meta.data$celltype)



CellDimPlot(seurat, group.by = "celltype", reduction = "UMAP")
ggsave("CellDimPlot.celltype.umap.pdf",width =9,height = 7)

CellDimPlot(seurat, group.by = "celltype", reduction = "TSNE")
ggsave("CellDimPlot.celltype.tsne.pdf",width =9,height = 7)


CellDimPlot(seurat, group.by = "orig.ident", reduction = "UMAP")
ggsave("orig.ident.pdf",width =9,height = 7)




Idents(seurat) <- "celltype"
genes_to_check = c("MITF","MLANA","TYR","PMEL","TYRP1","DCT", #Melanocytes
                   "KRT5","KRT14","TP63", #Keratinocytes
                   "DCN","COL1A1","COL3A1","PDGFRA", #Fibroblasts
                   "PECAM1","VWF","CLDN5","CDH5" #Endothelial
                   # "PMAIP1","GADD45B","TP53INP1","FOS","JUN","HSPA5","ATF3",
                   # "ATF4","DDIT3","PPP1R15A"#Apoptosis
)  


for (gene in genes_to_check) {
  p <- VlnPlot(seurat,features = gene) + ggtitle(gene)
  
  
  ggsave(
    filename = paste0(gene, ".png"),
    plot = p,
    width = 6, height = 5, units = "in", dpi = 300
  )
}


for (gene in genes_to_check) {
  p <- FeatureStatPlot(
    srt = seurat,
    group.by = "celltype",
    #bg.by = "celltype",
    stat.by = gene,
    add_box = TRUE
    # comparisons = list(
    #  c("T_NK", "B_Plasma"),
    #  c("Mast", "Myeloid"),
    #  c("Fibroblasts", "Endothelial"))
  ) + ggtitle(gene)
  

  ggsave(
    filename = paste0(gene, ".png"),
    plot = p,
    width = 6, height = 5, units = "in", dpi = 300
  )
}





DotPlot(seurat, features = genes_to_check) + RotatedAxis()
ggsave("DotPlot.celltype.pdf",width =9,height = 7)


DoHeatmap(seurat, features = genes_to_check, group.by = "celltype")
ggsave("DoHeatmap.celltype.pdf",width =9,height = 7)



Idents(seurat) <- "celltype"
table(Idents(seurat))

seurat_sub <- subset(seurat, downsample = 300)


p <- DoHeatmap(seurat_sub, 
               features = genes_to_check, 
               group.by = "celltype",
               size = 4,           # 调大上方分组文字的大小
               angle = 0,          # 调整文字角度，0度更易读
               label = TRUE,       # 显示标签
               draw.lines = TRUE)  # 在组间画分割线


ggsave("DoHeatmap_balanced.pdf", plot = p, width = 9, height = 7)


ht <- GroupHeatmap(
  srt = seurat,
  features = genes_to_check,
  group.by = c("celltype"),
  heatmap_palette = "YlOrRd",
  # cell_annotation = c("Phase", "G2M_score", "Cdh2"),
  # cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  # show_row_names = TRUE,
  row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)
pdf("GroupHeatmap_res.pdf", width = 10, height = 8)
print(ht$plot)
dev.off()


seurat$celltype <- factor(seurat$celltype,
                         levels = c("Melanocytes", "Keratinocytes",
                                    "Fibroblasts", "Endothelial"))
Idents(seurat) <- "celltype"

pdf("average_heatmap.pdf", width = 8, height = 10)
averageHeatmap(object = seurat,
               markerGene = genes_to_check,
               group.by = "celltype",
               cluster_rows = FALSE)
dev.off()

# decontX
library(decontX)
set.seed(123)
counts <- GetAssayData(object = seurat, layer = "counts")
decontX_res <- decontX(counts)
seurat$contamination <- decontX_res$contamination
table(seurat$contamination)

seuratnew = seurat[,seurat$contamination < 0.1]

DimPlot(seuratnew, label = T,group.by = "celltype")+DimPlot(seurat, label = T,group.by = "celltype")
ggsave("decontX.png",width = 14,height = 7)
seuratnew
seurat

# check
seurat$RNA_com <- ifelse(seurat$contamination < 0.1,
                        "RNA_None","RNA_com")


CellDimPlot(seurat, 
            group.by = "celltype", 
            reduction = "UMAP")
ggsave("final_celltype.umap.pdf",width =9,height = 7)

CellDimPlot(seurat, 
            group.by = "celltype", 
            reduction = "TSNE")
ggsave("final_celltype.tsne.pdf",width =9,height = 7)

qsave(seurat, "seurat_celltype.qs")

